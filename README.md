# GAMS solver link for NVIDIA cuOpt solver

[![Build cuOpt link for GAMS](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main-x86_64.yml/badge.svg)](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main-x86_64.yml) [![Build cuOpt link for GAMS](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main-arm64.yml/badge.svg)](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main-arm64.yml)

This project builds and packages the [GAMS](https://gams.com/) and [GAMSPy](https://gamspy.readthedocs.io/en/latest/index.html) solver link for the [NVIDIA cuOpt solver](https://github.com/NVIDIA/cuopt).

You can get more details and tips by reading the blog post ["GPU-Accelerated Optimization with GAMS and NVIDIA cuOpt"](https://www.gams.com/blog/2025/09/gpu-accelerated-optimization-with-gams-and-nvidia-cuopt/).

Supported model types are LP, MIP, RMIP, QCP, MIQCP, RMIQCP.

## Requirements

- **Operating System:** Linux, Windows 11 through WSL2
- **CPU architecture:** x86_64, arm64
- **GAMS:** Version 49 or newer
- **GAMSPy:** Version 1.12.1 or newer
- **NVIDIA GPU:** Volta architecture or better
- **CUDA Runtime Libraries:** 12 or 13

## Installation using `fetch-cuoptlink.py`

You can automatically download, install, test, and manage the cuOpt solver link using the provided `fetch-cuoptlink.py` script.

__Quickstart:__ Run the following commands to download and execute the script:

```bash
curl -O https://raw.githubusercontent.com/GAMS-dev/cuoptlink-builder/main/fetch-cuoptlink.py
pip install typer
python fetch-cuoptlink.py
```

which becomes a one-liner with [uv](https://docs.astral.sh/uv/):

```bash
uv run https://raw.githubusercontent.com/GAMS-dev/cuoptlink-builder/main/fetch-cuoptlink.py
```

### Interactive Mode

Running the script with no arguments launches an interactive prompt. It auto-detects your GAMS path (via `which gams`) and system CUDA version, prompting you for any missing options:

```bash
python fetch-cuoptlink.py
```

Calling `uninstall` without additional options will interactively prompt for the GAMS directory path:

```bash
python fetch-cuoptlink.py uninstall
```

### Non-interactive CLI Mode

You can also pass command-line arguments to automate installation and uninstallation:

```bash
# Basic installation using detected GAMS directory and CUDA runtime download
python fetch-cuoptlink.py install --gams-dir /opt/gams/gams55.0_linux_x64_64_sfx --cuda-runtime

# Specify a CUDA version and release tag explicitly
python fetch-cuoptlink.py install -g /opt/gams/gams55.0 -c 12 -r v0.0.8

# Uninstall the solver link from a GAMS system directory
python fetch-cuoptlink.py uninstall -g /opt/gams/gams55.0
```

> **Note:** Successful installations automatically verify the solver link by running the GAMS `trnsport` test model with `solver=cuopt`.

## Manual installation

- Make sure [CUDA runtime](https://developer.nvidia.com/cuda-downloads?target_os=Linux) is installed
- Download and unpack `cuopt-link-release-cu12-{x86_64,arm64}.zip` or `cuopt-link-release-cu13-{x86_64,arm64}.zip` (for CUDA 12 and 13 respectively) from the [releases page](https://github.com/GAMS-dev/cuoptlink-builder/releases):
    - Unpack the contents of `cuopt-link-release-cu*-*.zip` into your GAMS system directory. For GAMSPy, you can find out your system directory by running `gamspy show base`. So for example you can run `unzip -o cuopt-link-release-cu*-*.zip -d $(gamspy show base)`.
    - **Caution:** This will overwrite any existing `gamsconfig.yaml` file in that directory. The contained `gamsconfig.yaml` contains a `solverConfig` section to make cuOpt available to GAMS.

The neccessary files from the CUDA 12 or 13 runtime can also be downloaded as convenient archive `cu12-runtime-{x86_64,arm64}.zip` or `cu13-runtime-{x86_64,arm64}.zip` from the [releases page](https://github.com/GAMS-dev/cuoptlink-builder/releases).

## Test the setup

Get an example model and explicitly choose `cuopt` as `lp` or `mip` solver:
```
gamslib trnsport
gams trnsport lp cuopt
```

## Examples

### Notebooks

- [examples/trnsport_cuopt.ipynb](examples/trnsport_cuopt.ipynb) for CUDA 12 on x86_64
- [examples/trnsport_cuopt.ipynb](examples/trnsport_cuopt_cu13.ipynb) for CUDA 13 on x86_64

### GAMS models

Various GAMS models can be found in subfolder `examples/models` and are used to verify the solver link.