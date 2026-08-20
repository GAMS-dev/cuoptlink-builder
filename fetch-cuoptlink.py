from __future__ import annotations

import ctypes
import os
import platform
import shutil
import tempfile
import zipfile
from typing import Optional

import typer

app = typer.Typer(
    help="Manage cuOpt solver link for standalone GAMS system on Linux.",
    invoke_without_command=True,
)

SOLVER_NAME = "cuopt"
REPOSITORY = "GAMS-dev/cuoptlink-builder"
BASE_RELEASE_URL = f"https://api.github.com/repos/{REPOSITORY}/releases"
CUDA_VERSIONS = ("12", "13")
DEFAULT_CUDA_VERSION = "13"
CONFIG_FILE = "gamsconfig.yaml"
BACKUP_FILE = "gamsconfig.yaml.cuopt_backup"
MANIFEST_FILE = ".cuopt_installed_files.txt"


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
        typer.echo(f"`{SOLVER_NAME}` link installer only supports Linux systems.")
        raise typer.Exit(code=1)

    machine = platform.machine().lower()
    if machine in ("x86_64", "amd64"):
        return "x86_64"
    if machine in ("aarch64", "arm64"):
        return "arm64"

    typer.echo(f"`{SOLVER_NAME}` supports x86_64 and arm64, found {machine}.")
    raise typer.Exit(code=1)


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


def _get_asset_urls(names: list[str], release_tag: Optional[str] = None) -> list[str]:
    import requests

    if release_tag and release_tag.lower() != "latest":
        url = f"{BASE_RELEASE_URL}/tags/{release_tag}"
    else:
        url = f"{BASE_RELEASE_URL}/latest"

    try:
        response = requests.get(url, timeout=30)
    except requests.RequestException as e:
        typer.echo(f"Could not reach GitHub API ({url}): {e}")
        raise typer.Exit(code=1) from e

    if response.status_code != 200:
        typer.echo(
            f"Failed to fetch release info ({response.status_code}): {response.text}"
        )
        raise typer.Exit(code=1)

    release = response.json()
    assets = {
        asset["name"]: asset["browser_download_url"] for asset in release["assets"]
    }

    urls = []
    for name in names:
        if name not in assets:
            typer.echo(
                f"Release `{release['tag_name']}` does not contain `{name}`. "
                f"Available assets: {sorted(assets)}"
            )
            raise typer.Exit(code=1)
        urls.append(assets[name])

    typer.echo(f"Installing `{SOLVER_NAME}` from release `{release['tag_name']}`...")
    return urls


def _download(url: str, path: str) -> None:
    import requests
    from rich.progress import (
        BarColumn,
        DownloadColumn,
        Progress,
        TextColumn,
        TimeRemainingColumn,
    )

    name = os.path.basename(path)
    try:
        with requests.get(url, stream=True, timeout=60) as response:
            response.raise_for_status()
            total = int(response.headers.get("Content-Length", 0))

            with (
                Progress(
                    TextColumn("[progress.description]{task.description}"),
                    BarColumn(),
                    DownloadColumn(),
                    TimeRemainingColumn(),
                ) as progress,
                open(path, "wb") as file,
            ):
                task = progress.add_task(name, total=total or None)
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    _ = file.write(chunk)
                    progress.update(task, advance=len(chunk))
    except requests.RequestException as e:
        typer.echo(f"Could not download {url}: {e}")
        raise typer.Exit(code=1) from e


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
    typer.echo(f"Backed up original `{config_path}` to `{backup_path}`.")


def _run_interactive_install() -> None:
    default_gams_dir = _get_default_gams_dir()
    
    while True:
        prompt_kwargs = {}
        if default_gams_dir:
            prompt_kwargs["default"] = default_gams_dir

        gams_dir_input = typer.prompt("GAMS system directory path", **prompt_kwargs)
        gams_dir = os.path.abspath(os.path.expanduser(gams_dir_input))

        if os.path.isdir(gams_dir) and os.access(gams_dir, os.W_OK):
            break
        typer.echo(f"Error: Directory '{gams_dir}' does not exist or is not writable. Try again.\n")

    cuda_version = typer.prompt(
        "CUDA version",
        default=DEFAULT_CUDA_VERSION,
    )

    cuda_runtime = typer.confirm(
        "Download and install bundled CUDA runtime libraries?",
        default=False,
    )

    release = typer.prompt(
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
                typer.echo(
                    "CUDA runtime auto-detection failed. Pass `--cuda-version` "
                    "or enable `--cuda-runtime`."
                )
                raise typer.Exit(code=1)
            cuda_version = DEFAULT_CUDA_VERSION
        else:
            typer.echo(f"Detected CUDA runtime version {cuda_version}.")

    if cuda_version not in CUDA_VERSIONS:
        typer.echo(f"Unsupported CUDA version `{cuda_version}`. Choice: {CUDA_VERSIONS}")
        raise typer.Exit(code=1)

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
            files.update(_extract(archive_path, gams_dir))

    _set_installed_files(gams_dir, sorted(files))
    typer.echo(f"Successfully installed `{SOLVER_NAME}` into `{gams_dir}`.")


@app.callback(invoke_without_command=True)
def main(ctx: typer.Context) -> None:
    if ctx.invoked_subcommand is None:
        _run_interactive_install()


@app.command()
def install(
    gams_dir: str = typer.Option(
        ...,
        "--gams-dir",
        "-g",
        help="Path to the GAMS system directory.",
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=True,
        resolve_path=True,
    ),
    cuda_version: Optional[str] = typer.Option(
        None,
        "--cuda-version",
        "-c",
        help="CUDA version (12 or 13). Auto-detected if omitted.",
    ),
    cuda_runtime: bool = typer.Option(
        False,
        "--cuda-runtime",
        help="Download and unpack bundled CUDA runtime libraries.",
    ),
    release: Optional[str] = typer.Option(
        "latest",
        "--release",
        "-r",
        help="Tag name of cuoptlink-builder release (e.g. v1.0.0 or 'latest').",
    ),
) -> None:
    """Install the cuOpt link into a target GAMS system directory."""
    _execute_install(gams_dir, cuda_version, cuda_runtime, release)


@app.command()
def uninstall(
    gams_dir: str = typer.Option(
        ...,
        "--gams-dir",
        "-g",
        help="Path to the GAMS system directory.",
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=True,
        resolve_path=True,
    ),
) -> None:
    """Uninstall the cuOpt link from a target GAMS system directory."""
    _get_architecture()
    installed_files = get_installed_files(gams_dir)
    if not installed_files:
        typer.echo(f"`{SOLVER_NAME}` is not recorded as installed in `{gams_dir}`.")
        raise typer.Exit(code=1)

    for name in installed_files:
        try:
            os.unlink(os.path.join(gams_dir, name))
        except FileNotFoundError:
            pass

    backup_path = os.path.join(gams_dir, BACKUP_FILE)
    if os.path.isfile(backup_path):
        shutil.move(backup_path, os.path.join(gams_dir, CONFIG_FILE))

    try:
        os.unlink(_get_manifest_path(gams_dir))
    except FileNotFoundError:
        pass

    typer.echo(f"Successfully uninstalled `{SOLVER_NAME}` from `{gams_dir}`.")


if __name__ == "__main__":
    app()