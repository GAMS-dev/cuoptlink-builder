#!/usr/bin/env python3
"""Generate the option table of the GAMS solver manual from a GAMS option definition file.

Reads e.g. assets/optcuopt.def and writes optcuopt_s.md in the format that the GAMS
documentation includes into the solver manual page (same layout as the tables written by
GamsOptionsSpecWriter for other solver links): one subsection per option group, one table row
per visible option, with doxygen anchors <SOLVER><option name> to link to options.

Usage: python3 optdef2md.py [assets/optcuopt.def] [-o optcuopt_s.md] [--source "<repo> <tag>"]
"""

import argparse
import re
import sys
from pathlib import Path

INFINITE = {"maxint": "&infin;", "maxdouble": "&infin;", "minint": "-&infin;", "mindouble": "-&infin;"}

# Options taking their default from a GAMS parameter say so in their description,
# e.g. "(default GAMS ResLim)" or "(default GAMS OptCR, at most 0.1)".
GAMS_DEFAULT = re.compile(r"\s*\(default GAMS (\w+)[^)]*\)")

TOKEN = re.compile(r'\s*("[^"]*"|\S+)')

# number of fields between type and description, per option type
NFIELDS = {"integer": 6, "double": 6, "boolean": 4, "string": 4, "enumint": 4, "enumstr": 4}


def split_fields(line, n):
    """Split off n whitespace-separated (possibly double-quoted) fields, return them and the rest."""
    fields = []
    pos = 0
    for _ in range(n):
        m = TOKEN.match(line, pos)
        if not m:
            raise ValueError(f"expected {n} fields in line: {line}")
        fields.append(m.group(1))
        pos = m.end()
    return fields, line[pos:].strip()


def escape(s):
    """Make text safe for a doxygen markdown table cell (as makeValidMarkdownString in GamsOptionsSpecWriter)."""
    for old, new in (("<=", "&le;"), (">=", "&ge;"), ("<", "&lt;"), (">", "&gt;"), ("|", "\\|"), ("$", "\\f$")):
        s = s.replace(old, new)
    return s


def number(s):
    """Format a numeric field of the def file like the C++ writer does (%g), with infinities as symbols."""
    if s in INFINITE:
        return INFINITE[s]
    v = float(s)
    return str(int(v)) if v.is_integer() and abs(v) < 1e15 and "e" not in s.lower() else "%g" % v


def parse_def(path):
    groups = {}  # group number -> (name, title)
    options = []
    current = None  # last enum option, receives the indented enum value lines
    for lineno, raw in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw.rstrip()
        if not line.strip() or line.startswith("*"):
            continue
        try:
            if line[0].isspace():
                if current is None:
                    raise ValueError("enum value without enum option")
                (value, visible), text = split_fields(line, 2)
                if visible != "0":
                    current["enum"].append((value.strip('"'), text))
                continue
            name, kind = line.split()[:2]
            current = None
            if kind == "group":
                (_, _, number_, _), title = split_fields(line, 4)
                groups[int(number_)] = (name, title)
                continue
            if kind not in NFIELDS:
                raise ValueError(f"unsupported option type '{kind}'")
            fields, descr = split_fields(line, 2 + NFIELDS[kind])
            fields = fields[2:]
            opt = {"name": name, "type": kind, "default": fields[1].strip('"'),
                   "visible": fields[-2] != "0", "group": int(fields[-1]), "descr": descr, "enum": []}
            if kind in ("integer", "double"):
                opt["low"], opt["high"] = fields[2], fields[3]
            if kind.startswith("enum"):
                current = opt
            options.append(opt)
        except ValueError as e:
            sys.exit(f"{path}:{lineno}: {e}")
    return groups, options


def table_row(solver, opt):
    descr = opt["descr"]
    default = opt["default"]
    m = GAMS_DEFAULT.search(descr)
    if m:
        descr = descr[:m.start()] + descr[m.end():]
        param = m.group(1).lower()
        default = f'\\ref GAMSAO{param} "GAMS {param}"'
    elif opt["type"] in ("integer", "double"):
        default = number(default)

    cell = escape(descr.strip())
    kind = opt["type"]
    if opt["enum"]:
        for value, text in opt["enum"]:
            cell += f"<br/>{escape(value)}: {escape(text)}" if text else f"<br/>{escape(value)}"
    elif kind == "integer":
        cell += f"<br/>Range: {{{number(opt['low'])}, ..., {number(opt['high'])}}}"
    elif kind == "double":
        low, high = number(opt["low"]), number(opt["high"])
        cell += "<br/>Range: real" if (low, high) == ("-&infin;", "&infin;") else f"<br/>Range: [{low}, {high}]"
    elif kind == "boolean":
        cell += "<br/>Range: boolean"
    elif kind == "string":
        cell += "<br/>Range: string"

    return f"| \\anchor {solver.upper()}{opt['name']} {opt['name']} | {cell} | {escape(default)}|"


def write_md(solver, deffile, groups, options, source=None):
    origin = f"{Path(deffile).name} of {source}" if source else Path(deffile).name
    lines = [f"<!-- This file is autogenerated by optdef2md.py from {origin}. -->"]
    for number_, (name, title) in sorted(groups.items()):
        opts = sorted((o for o in options if o["group"] == number_ and o["visible"]), key=lambda o: o["name"].lower())
        if not opts:
            continue
        lines.append("")
        if len(groups) > 1:
            gid = re.sub(r"\W", "_", name)
            lines += [f"\\subsection {solver.upper()}_gr_{gid} {title}", ""]
        lines += ["| Option | Description | Default |", "|:-------|:------------| ------: |"]
        lines += [table_row(solver, o) for o in opts]
    unknown = sorted({o["group"] for o in options} - set(groups))
    if unknown:
        sys.exit(f"{deffile}: options refer to undefined group(s) {unknown}")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("deffile", nargs="?", default=str(Path(__file__).parent / "assets" / "optcuopt.def"),
                        help="option definition file opt<solver>.def (default: %(default)s)")
    parser.add_argument("-o", "--output", help="output file (default: opt<solver>_s.md in the current directory)")
    parser.add_argument("--source", help="origin of the def file recorded in the output header, "
                                         "e.g. \"GAMS-dev/cuoptlink-builder v0.0.9b\"")
    args = parser.parse_args()
    if args.source and "--" in args.source:
        sys.exit("--source must not contain '--' (it goes into an HTML comment)")

    stem = Path(args.deffile).stem
    if not stem.startswith("opt"):
        sys.exit(f"{args.deffile}: expected a file named opt<solver>.def")
    solver = stem[3:]
    groups, options = parse_def(args.deffile)
    output = Path(args.output or f"opt{solver}_s.md")
    output.write_text(write_md(solver, args.deffile, groups, options, args.source))
    print(f"Wrote {output} ({sum(o['visible'] for o in options)} options)")


if __name__ == "__main__":
    main()
