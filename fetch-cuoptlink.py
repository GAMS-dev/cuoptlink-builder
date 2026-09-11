from __future__ import annotations

import argparse
import ctypes
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
import zipfile
from typing import Optional, Tuple

SOLVER_NAME = "cuopt"
REPOSITORY = "GAMS-dev/cuoptlink-builder"
BASE_RELEASE_URL = f"https://api.github.com/repos/{REPOSITORY}/releases"
CUDA_VERSIONS = ("12", "13")
DEFAULT_CUDA_VERSION = "13"
CONFIG_FILE = "gamsconfig.yaml"
CUOPT_CONFIG_FILE = "gamsconfig_cuopt.yaml"
BACKUP_FILE = "gamsconfig.yaml.cuopt_backup"
MANIFEST_FILE = ".cuopt_installed_files.txt"
SOLVER_CONFIG_SECTION = "solverConfig"
MARKER_BEGIN = f"# begin {SOLVER_NAME} solver configuration of cuoptlink"
MARKER_END = f"# end {SOLVER_NAME} solver configuration of cuoptlink"


def _get_manifest_path(gams_dir: str) -> str:
    return os.path.join(gams_dir, MANIFEST_FILE)


def get_installed_files(gams_dir: str) -> list[str]:
    try:
        with open(_get_manifest_path(gams_dir), encoding="utf-8") as file:
            return [line for line in file.read().splitlines() if line]
    except FileNotFoundError:
        return []


def _set_installed_files(gams_dir: str, files: list[str]) -> None:
    with open(_get_manifest_path(gams_dir), "w", encoding="utf-8") as file:
        file.write("\n".join(sorted(files)) + "\n")


def _get_architecture() -> str:
    if platform.system() != "Linux":
        print(f"`{SOLVER_NAME}` link installer only supports Linux systems.")
        raise SystemExit(1)

    machine = platform.machine().lower()
    if machine in ("x86_64", "amd64"):
        return "x86_64"
    if machine in ("aarch64", "arm64"):
        return "arm64"

    print(f"`{SOLVER_NAME}` supports x86_64 and arm64, found {machine}.")
    raise SystemExit(1)


def _detect_cuda_version() -> str | None:
    for version in reversed(CUDA_VERSIONS):
        for library in (f"libcudart.so.{version}", f"libcublas.so.{version}"):
            try:
                _ = ctypes.CDLL(library)
                return version
            except OSError:
                continue
    return None


def _get_default_gams_dir() -> Optional[str]:
    gams_bin = shutil.which("gams")
    if gams_bin:
        resolved_bin = os.path.realpath(gams_bin)
        return os.path.dirname(resolved_bin)
    return None


def _prompt(text: str, default: Optional[str] = None) -> str:
    suffix = f" [{default}]" if default is not None else ""
    while True:
        value = input(f"{text}{suffix}: ").strip()
        if value:
            return value
        if default is not None:
            return default


def _confirm(text: str, default: bool) -> bool:
    choices = "Y/n" if default else "y/N"
    while True:
        value = input(f"{text} [{choices}]: ").strip().lower()
        if not value:
            return default
        if value in ("y", "yes"):
            return True
        if value in ("n", "no"):
            return False
        print("Please answer 'y' or 'n'.")


def _prompt_gams_dir() -> str:
    default_gams_dir = _get_default_gams_dir()
    while True:
        gams_dir_input = _prompt("GAMS system directory path", default=default_gams_dir)
        gams_dir = os.path.abspath(os.path.expanduser(gams_dir_input))

        if os.path.isdir(gams_dir) and os.access(gams_dir, os.W_OK):
            return gams_dir
        print(f"Error: Directory '{gams_dir}' does not exist or is not writable. Try again.\n")


def _prompt_cuda_version() -> Tuple[str, bool]:
    detected_cuda = _detect_cuda_version()
    if detected_cuda:
        print(f"Detected CUDA runtime on system: CUDA {detected_cuda}")
        default_cuda = detected_cuda
        has_detected_cuda = True
    else:
        print("No CUDA runtime automatically detected on system.")
        default_cuda = DEFAULT_CUDA_VERSION
        has_detected_cuda = False

    versions_str = ", ".join(CUDA_VERSIONS)
    cuda_version = _prompt(
        f"CUDA version (supported: {versions_str})",
        default=default_cuda,
    )

    # If CUDA was detected on the host system, default runtime download to False.
    # Otherwise, default runtime download to True.
    default_runtime_download = not has_detected_cuda

    cuda_runtime = _confirm(
        "Download and install bundled CUDA runtime libraries?",
        default=default_runtime_download,
    )

    return cuda_version, cuda_runtime


def _get_asset_urls(names: list[str], release_tag: Optional[str] = None) -> list[str]:
    if release_tag and release_tag.lower() != "latest":
        url = f"{BASE_RELEASE_URL}/tags/{release_tag}"
    else:
        url = f"{BASE_RELEASE_URL}/latest"

    request = urllib.request.Request(url, headers={"Accept": "application/vnd.github+json"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            body = response.read()
    except urllib.error.HTTPError as e:
        print(f"Failed to fetch release info ({e.code}): {e.read().decode(errors='replace')}")
        raise SystemExit(1) from e
    except urllib.error.URLError as e:
        print(f"Could not reach GitHub API ({url}): {e}")
        raise SystemExit(1) from e

    release = json.loads(body)
    assets = {
        asset["name"]: asset["browser_download_url"] for asset in release["assets"]
    }

    urls = []
    for name in names:
        if name not in assets:
            print(
                f"Release `{release['tag_name']}` does not contain `{name}`. "
                f"Available assets: {sorted(assets)}"
            )
            raise SystemExit(1)
        urls.append(assets[name])

    print(f"Installing `{SOLVER_NAME}` from release `{release['tag_name']}`...")
    return urls


def _format_size(num_bytes: int) -> str:
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024:
            return f"{size:.1f}{unit}"
        size /= 1024
    return f"{size:.1f}TB"


def _download(url: str, path: str) -> None:
    name = os.path.basename(path)
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            total = int(response.headers.get("Content-Length", 0))

            with open(path, "wb") as file:
                downloaded = 0
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    _ = file.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = downloaded * 100 // total
                        sys.stdout.write(
                            f"\r{name}: {pct:3d}% ({_format_size(downloaded)}/{_format_size(total)})"
                        )
                    else:
                        sys.stdout.write(f"\r{name}: {_format_size(downloaded)}")
                    sys.stdout.flush()
                sys.stdout.write("\n")
    except (urllib.error.HTTPError, urllib.error.URLError) as e:
        print(f"Could not download {url}: {e}")
        raise SystemExit(1) from e


def _extract(path: str, directory: str) -> list[str]:
    names = []
    with zipfile.ZipFile(path) as archive:
        for info in archive.infolist():
            if info.is_dir():
                continue

            name = os.path.basename(info.filename)
            if not name:
                continue

            target_path = os.path.join(directory, name)
            with archive.open(info) as source, open(target_path, "wb") as target:
                shutil.copyfileobj(source, target)

            mode = (info.external_attr >> 16) & 0o777
            if mode:
                os.chmod(target_path, mode)

            names.append(name)
    return names


def _backup_config(gams_dir: str, installed_files: list[str]) -> None:
    config_path = os.path.join(gams_dir, CONFIG_FILE)
    if not os.path.isfile(config_path) or CONFIG_FILE in installed_files:
        return

    backup_path = os.path.join(gams_dir, BACKUP_FILE)
    shutil.copy2(config_path, backup_path)
    print(f"Backed up original `{config_path}` to `{backup_path}`.")


_TOP_LEVEL_KEY_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_.\-]*)[ \t]*:([ \t]|$)")
_ENTRY_KEY_RE = re.compile(r"[ \t]*([A-Za-z_][A-Za-z0-9_.\-]*)[ \t]*:([ \t]|$)")


def _get_top_level_keys(lines: list[str], path: str) -> list[str]:
    keys: list[str] = []
    for line in lines[: _get_document_end(lines)]:
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or re.match(r"^---([ \t]|$)", line):
            continue
        if line[0].isspace() or stripped == "-" or stripped.startswith("- "):
            continue

        match = _TOP_LEVEL_KEY_RE.match(line)
        if not match:
            print(
                f"`{path}` must map configuration sections such as "
                f"`{SOLVER_CONFIG_SECTION}` to their entries. Could not parse line: "
                f"`{line}`."
            )
            raise SystemExit(1)
        keys.append(match.group(1))

    return keys


def _load_top_level_keys(path: str) -> list[str]:
    with open(path, encoding="utf-8") as file:
        lines = file.read().splitlines()
    return _get_top_level_keys(lines, path)


def _solver_config_has_entry(lines: list[str], start: int, end: int, name: str) -> bool:
    entries = lines[start:end]
    dash_indent = _get_entry_indent(entries)

    if dash_indent is not None:
        for entry in entries:
            if not entry.startswith(dash_indent) or not entry[len(dash_indent) :].startswith("-"):
                continue
            match = _ENTRY_KEY_RE.match(entry[len(dash_indent) + 1 :])
            if match and match.group(1) == name:
                return True
        return False

    base_indent: Optional[str] = None
    for entry in entries:
        stripped = entry.strip()
        if not stripped or stripped.startswith("#"):
            continue
        indent = entry[: len(entry) - len(entry.lstrip(" "))]
        if base_indent is None:
            base_indent = indent
        if indent != base_indent:
            continue
        match = _ENTRY_KEY_RE.match(entry[len(base_indent) :])
        if match and match.group(1) == name:
            return True
    return False


def _find_section(lines: list[str], section: str, path: str) -> int | None:
    for index, line in enumerate(lines):
        if not re.match(rf"{re.escape(section)}[ \t]*:", line):
            continue

        if not re.match(rf"{re.escape(section)}[ \t]*:[ \t]*(#.*)?$", line):
            print(
                f"`{section}` of `{path}` is not written as a block of entries, "
                f"hence we cannot add `{SOLVER_NAME}` to it. Please add the "
                f"entry by hand."
            )
            raise SystemExit(1)

        return index

    return None


def _get_section_end(lines: list[str], start: int) -> int:
    end = start
    for index in range(start, len(lines)):
        line = lines[index]
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        is_entry = stripped == "-" or stripped.startswith("- ")
        if not line[0].isspace() and not is_entry:
            break

        end = index + 1

    return end


def _get_document_end(lines: list[str]) -> int:
    for index, line in enumerate(lines):
        if re.match(r"\.\.\.([ \t]|$)", line):
            return index

    return len(lines)


def _get_entry_indent(entries: list[str]) -> str | None:
    for entry in entries:
        match = re.match(r"([ \t]*)-([ \t]|$)", entry)
        if match is not None:
            return match.group(1)

    return None


def _reindent(entries: list[str], indent: str) -> list[str]:
    base = _get_entry_indent(entries) or ""
    return [
        indent + entry[len(base) :] if entry.startswith(base) else entry
        for entry in entries
    ]


def _get_solver_config(gams_dir: str) -> list[str]:
    path = os.path.join(gams_dir, CUOPT_CONFIG_FILE)
    if _load_top_level_keys(path) != [SOLVER_CONFIG_SECTION]:
        print(
            f"`{path}` of the release archive holds more than a "
            f"`{SOLVER_CONFIG_SECTION}` section, which cannot be merged into "
            f"`{CONFIG_FILE}`. Please add its content to `{CONFIG_FILE}` by hand."
        )
        raise SystemExit(1)

    with open(path, encoding="utf-8") as file:
        lines = file.read().splitlines()

    index = _find_section(lines, SOLVER_CONFIG_SECTION, path)
    if index is None:
        print(f"`{path}` of the release archive does not register a solver.")
        raise SystemExit(1)

    return lines[index + 1 : _get_section_end(lines, index + 1)]


def _merge_solver_config(text: str, entries: list[str], path: str) -> str:
    lines = text.splitlines()
    index = _find_section(lines, SOLVER_CONFIG_SECTION, path)
    if index is None:
        block = [f"{SOLVER_CONFIG_SECTION}:", *entries]
        end = _get_document_end(lines)
    else:
        end = _get_section_end(lines, index + 1)
        indent = _get_entry_indent(lines[index + 1 : end])
        block = entries if indent is None else _reindent(entries, indent)

    merged = [*lines[:end], MARKER_BEGIN, *block, MARKER_END, *lines[end:]]
    return "\n".join(merged) + "\n"


def _strip_solver_config(text: str) -> str:
    remaining = []
    is_contributed = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped == MARKER_BEGIN:
            is_contributed = True
        elif stripped == MARKER_END:
            is_contributed = False
        elif not is_contributed:
            remaining.append(line)

    return "\n".join(remaining) + "\n" if remaining else ""


def _merge_cuopt_config(gams_dir: str) -> None:
    cuopt_cfg = os.path.join(gams_dir, CUOPT_CONFIG_FILE)
    if not os.path.isfile(cuopt_cfg):
        return

    try:
        config_path = os.path.join(gams_dir, CONFIG_FILE)

        if os.path.isfile(config_path):
            has_cuopt = False
            if SOLVER_CONFIG_SECTION in _load_top_level_keys(config_path):
                with open(config_path, encoding="utf-8") as file:
                    existing_lines = file.read().splitlines()
                index = _find_section(existing_lines, SOLVER_CONFIG_SECTION, config_path)
                if index is not None:
                    end = _get_section_end(existing_lines, index + 1)
                    has_cuopt = _solver_config_has_entry(existing_lines, index + 1, end, SOLVER_NAME)

            if has_cuopt:
                print(f"`{SOLVER_NAME}` entry already exists in `{CONFIG_FILE}`. Skipping merge.")
                return

        entries = _get_solver_config(gams_dir)

        text = ""
        try:
            with open(config_path, encoding="utf-8") as file:
                text = _strip_solver_config(file.read())
        except FileNotFoundError:
            pass

        with open(config_path, "w", encoding="utf-8") as file:
            file.write(_merge_solver_config(text, entries, config_path))

        print(f"Merged `{CUOPT_CONFIG_FILE}` into `{CONFIG_FILE}`.")
    except Exception as e:
        print(f"Failed to merge config:\n{e}")
        raise SystemExit(1) from e
    finally:
        if os.path.exists(cuopt_cfg):
            os.unlink(cuopt_cfg)


def _remove_cuopt_from_config(gams_dir: str) -> None:
    config_path = os.path.join(gams_dir, CONFIG_FILE)
    try:
        with open(config_path, encoding="utf-8") as file:
            remaining = _strip_solver_config(file.read())
    except FileNotFoundError:
        return

    try:
        if not remaining.strip():
            os.unlink(config_path)
        else:
            with open(config_path, "w", encoding="utf-8") as file:
                file.write(remaining)

        print(f"Removed `{SOLVER_NAME}` entry from `{CONFIG_FILE}`.")
    except Exception as e:
        print(f"Warning: Failed to remove `{SOLVER_NAME}` from `{CONFIG_FILE}`: {e}")


def _test_installation(gams_dir: str) -> None:
    print("Testing installation with GAMS trnsport model...")
    gamslib_bin = os.path.join(gams_dir, "gamslib")
    gams_bin = os.path.join(gams_dir, "gams")

    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Fetch model 1 (trnsport) using gamslib
            subprocess.run(
                [gamslib_bin, "-q", "1"],
                cwd=temp_dir,
                check=True,
            )

            # Execute gams with cuopt solver
            subprocess.run(
                [gams_bin, "trnsport", "solver=cuopt", "lo=0"],
                cwd=temp_dir,
                check=True,
            )
            print("Installation test completed successfully!")
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"Installation test failed: {e}")
            raise SystemExit(1) from e


def _run_interactive_install() -> None:
    gams_dir = _prompt_gams_dir()
    cuda_version, cuda_runtime = _prompt_cuda_version()

    release = _prompt(
        "cuoptlink-builder release version",
        default="latest",
    )

    _execute_install(gams_dir, cuda_version, cuda_runtime, release)


def _execute_install(
    gams_dir: str,
    cuda_version: Optional[str],
    cuda_runtime: bool,
    release: Optional[str],
) -> None:
    architecture = _get_architecture()

    if cuda_version is None:
        cuda_version = _detect_cuda_version()
        if cuda_version is None:
            if not cuda_runtime:
                print(
                    "CUDA runtime auto-detection failed. Pass `--cuda-version` "
                    "or enable `--cuda-runtime`."
                )
                raise SystemExit(1)
            cuda_version = DEFAULT_CUDA_VERSION
        else:
            print(f"Detected CUDA runtime version {cuda_version}.")

    if cuda_version not in CUDA_VERSIONS:
        print(f"Unsupported CUDA version `{cuda_version}`. Choice: {CUDA_VERSIONS}")
        raise SystemExit(1)

    archives = [f"cuopt-link-release-cu{cuda_version}-{architecture}.zip"]
    if cuda_runtime:
        archives.append(f"cu{cuda_version}-runtime-{architecture}.zip")

    urls = _get_asset_urls(archives, release_tag=release)
    installed_files = get_installed_files(gams_dir)
    _backup_config(gams_dir, installed_files)

    files = set(installed_files)
    with tempfile.TemporaryDirectory() as directory:
        for name, url in zip(archives, urls):
            archive_path = os.path.join(directory, name)
            _download(url, archive_path)
            extracted_files = _extract(archive_path, gams_dir)
            files.update(extracted_files)

    _merge_cuopt_config(gams_dir)

    files.discard(CUOPT_CONFIG_FILE)
    files.discard(CONFIG_FILE)

    _set_installed_files(gams_dir, sorted(files))
    print(f"Successfully installed `{SOLVER_NAME}` into `{gams_dir}`.")

    _test_installation(gams_dir)


def _execute_uninstall(gams_dir: str) -> None:
    _get_architecture()
    installed_files = get_installed_files(gams_dir)
    if not installed_files:
        print(f"`{SOLVER_NAME}` is not recorded as installed in `{gams_dir}`.")
        raise SystemExit(1)

    for name in installed_files:
        if name == CONFIG_FILE:
            continue
        try:
            os.unlink(os.path.join(gams_dir, name))
        except FileNotFoundError:
            pass

    _remove_cuopt_from_config(gams_dir)

    try:
        os.unlink(_get_manifest_path(gams_dir))
    except FileNotFoundError:
        pass

    print(f"Successfully uninstalled `{SOLVER_NAME}` from `{gams_dir}`.")


def _install_command(args: argparse.Namespace) -> None:
    gams_dir = args.gams_dir
    if gams_dir is None:
        gams_dir = _prompt_gams_dir()
    _execute_install(gams_dir, args.cuda_version, args.cuda_runtime, args.release)


def _uninstall_command(args: argparse.Namespace) -> None:
    gams_dir = args.gams_dir
    if gams_dir is None:
        gams_dir = _prompt_gams_dir()
    _execute_uninstall(gams_dir)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Manage cuOpt solver link for standalone GAMS system on Linux.",
    )
    subparsers = parser.add_subparsers(dest="command")

    install_parser = subparsers.add_parser(
        "install", help="Install the cuOpt link into a target GAMS system directory."
    )
    install_parser.add_argument(
        "--gams-dir", "-g", dest="gams_dir", default=None,
        help="Path to the GAMS system directory.",
    )
    install_parser.add_argument(
        "--cuda-version", "-c", dest="cuda_version", default=None,
        help="CUDA version (12 or 13). Auto-detected if omitted.",
    )
    install_parser.add_argument(
        "--cuda-runtime", dest="cuda_runtime", action="store_true", default=False,
        help="Download and unpack bundled CUDA runtime libraries.",
    )
    install_parser.add_argument(
        "--release", "-r", dest="release", default="latest",
        help="Tag name of cuoptlink-builder release (e.g. v1.0.0 or 'latest').",
    )
    install_parser.set_defaults(func=_install_command)

    uninstall_parser = subparsers.add_parser(
        "uninstall", help="Uninstall the cuOpt link from a target GAMS system directory."
    )
    uninstall_parser.add_argument(
        "--gams-dir", "-g", dest="gams_dir", default=None,
        help="Path to the GAMS system directory.",
    )
    uninstall_parser.set_defaults(func=_uninstall_command)

    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()
    if args.command is None:
        _run_interactive_install()
    else:
        args.func(args)


if __name__ == "__main__":
    main()