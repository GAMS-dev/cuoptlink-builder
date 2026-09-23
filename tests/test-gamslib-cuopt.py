"""Regression test: solve gamslib models with cuOpt and compare against CPLEX.

Each model listed in baseline.txt (name, model type, CPLEX objective) is
extracted with `gamslib`, solved with `solver=cuopt`, and the objective value
reported in the GAMS trace file is compared against the CPLEX reference.

Run with:
    python3 tests/test-gamslib-cuopt.py [-g GAMS_DIR] [-j JOBS] [model ...]

Exits with status 1 if any model mismatches or fails to produce a solution.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

_DEFAULT_BASELINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "baseline.txt")

# Column indices in a GAMS trace record written with traceopt=3.
_TRC_MODEL_STATUS = 13
_TRC_SOLVER_STATUS = 14
_TRC_OBJ = 15
_TRC_SOLVER_TIME = 17


@dataclass
class Result:
    name: str
    model_type: str
    ref: float
    obj: float | None = None
    model_status: str = "NA"
    solver_status: str = "NA"
    solver_time: str = "NA"
    error: str | None = None

    def rel_diff(self) -> float | None:
        if self.obj is None:
            return None
        return abs(self.obj - self.ref) / max(1.0, abs(self.ref))


def read_baseline(path: str) -> list[tuple[str, str, float]]:
    entries = []
    with open(path) as f:
        for line in f:
            fields = line.split()
            if len(fields) == 3:
                entries.append((fields[0], fields[1], float(fields[2])))
    return entries


def find_gams_dir(gams_dir: str | None) -> str:
    if gams_dir is None:
        gams = shutil.which("gams")
        if gams is None:
            sys.exit("error: gams not found in PATH, use --gams-dir")
        gams_dir = os.path.dirname(os.path.realpath(gams))
    if not os.path.isfile(os.path.join(gams_dir, "gams")):
        sys.exit(f"error: no gams executable in {gams_dir}")
    return gams_dir


def run_model(gams_dir: str, work_dir: str, name: str, model_type: str, ref: float,
              reslim: int, timeout: int) -> Result:
    res = Result(name, model_type, ref)
    model_dir = os.path.join(work_dir, name)
    os.makedirs(model_dir, exist_ok=True)
    try:
        subprocess.run([os.path.join(gams_dir, "gamslib"), "-q", name], cwd=model_dir,
                       check=True, capture_output=True)
        proc = subprocess.run(
            [os.path.join(gams_dir, "gams"), name, "solver=cuopt", "optcr=0", "optca=0",
             f"reslim={reslim}", "trace=trc.txt", "traceopt=3", "lo=2"],
            cwd=model_dir, capture_output=True, timeout=timeout)
    except subprocess.CalledProcessError:
        res.error = "gamslib failed"
        return res
    except subprocess.TimeoutExpired:
        res.error = f"timeout after {timeout}s"
        return res

    trc_path = os.path.join(model_dir, "trc.txt")
    records = []
    if os.path.isfile(trc_path):
        with open(trc_path) as f:
            records = [l.strip().split(",") for l in f if l.strip() and not l.startswith("*")]
    if not records:
        res.error = f"no solve recorded (gams rc={proc.returncode})"
        return res

    rec = records[-1]
    res.model_status = rec[_TRC_MODEL_STATUS]
    res.solver_status = rec[_TRC_SOLVER_STATUS]
    res.solver_time = rec[_TRC_SOLVER_TIME]
    try:
        res.obj = float(rec[_TRC_OBJ])
    except ValueError:
        res.error = "no objective value"
    if proc.returncode != 0 and res.error is None:
        res.error = f"gams rc={proc.returncode}"
    return res


def main() -> int:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("models", nargs="*", help="subset of models to run (default: all)")
    parser.add_argument("-g", "--gams-dir", help="GAMS system directory (default: from PATH)")
    parser.add_argument("-b", "--baseline", default=_DEFAULT_BASELINE, help="reference file")
    parser.add_argument("-r", "--reslim", type=int, default=120, help="time limit per solve in s")
    parser.add_argument("-j", "--jobs", type=int, default=1, help="models to run in parallel")
    parser.add_argument("-t", "--rtol", type=float, default=1e-6,
                        help="tolerance for |obj-ref|/max(1,|ref|)")
    parser.add_argument("-k", "--keep", action="store_true", help="keep the work directory")
    args = parser.parse_args()

    gams_dir = find_gams_dir(args.gams_dir)
    entries = read_baseline(args.baseline)
    if args.models:
        unknown = set(args.models) - {e[0] for e in entries}
        if unknown:
            sys.exit(f"error: not in baseline: {' '.join(sorted(unknown))}")
        entries = [e for e in entries if e[0] in args.models]

    work_dir = tempfile.mkdtemp(prefix="gamslib-cuopt-")
    # Some models set their own resLim, so leave generous headroom beyond --reslim.
    timeout = 2 * args.reslim + 60

    row = "{:6} {:12} {:6} {:>9} {:>9} {:>9} {:>17} {:>17} {}"
    print(f"cuOpt vs. CPLEX on {len(entries)} gamslib models "
          f"(gams={gams_dir}, reslim={args.reslim}, rtol={args.rtol:g})")
    header = row.format("Result", "Model", "Type", "ModelStat", "SolveStat", "Time",
                        "cuOpt Obj", "CPLEX Obj", "RelDiff")
    print(header)
    print("-" * len(header), flush=True)

    def job(entry):
        name, model_type, ref = entry
        res = run_model(gams_dir, work_dir, name, model_type, ref, args.reslim, timeout)
        rel = res.rel_diff()
        ok = res.error is None and rel is not None and rel <= args.rtol
        obj = "NA" if res.obj is None else f"{res.obj:.10g}"
        detail = res.error or f"{rel:.2e}"
        print(row.format("OK" if ok else "FAIL", res.name, res.model_type, res.model_status,
                         res.solver_status, res.solver_time, obj, f"{res.ref:.10g}", detail),
              flush=True)
        return res, ok

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(job, entries))

    failed = [r for r, ok in results if not ok]
    print(f"\n{len(results) - len(failed)}/{len(results)} models match the CPLEX reference")
    if failed:
        print("Mismatches / failures:")
        for r in failed:
            rel = r.rel_diff()
            detail = r.error or f"cuopt={r.obj:.10g} cplex={r.ref:.10g} rel={rel:.2e}"
            print(f"  {r.name:12} ms={r.model_status:3} {detail}")

    if args.keep or failed:
        print(f"\nWork directory: {work_dir}")
    else:
        shutil.rmtree(work_dir, ignore_errors=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
