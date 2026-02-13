# GAMS solver link for NVIDIA cuOpt solver

[![Build cuOpt link for GAMS](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main-x86_64.yml/badge.svg)](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main-x86_64.yml) [![Build cuOpt link for GAMS](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main-arm64.yml/badge.svg)](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main-arm64.yml)

This project builds and packages the [GAMS](https://gams.com/) and [GAMSPy](https://gamspy.readthedocs.io/en/latest/index.html) solver link for the [NVIDIA cuOpt solver](https://github.com/NVIDIA/cuopt).

You can get more details and tips by reading the blog post ["GPU-Accelerated Optimization with GAMS and NVIDIA cuOpt"](https://www.gams.com/blog/2025/09/gpu-accelerated-optimization-with-gams-and-nvidia-cuopt/).

## Requirements

- **Operating System:** Linux, Windows 11 through WSL2
- **CPU architecture:** x86_64, arm64
- **GAMS:** Version 49 or newer
- **GAMSPy:** Version 1.12.1 or newer
- **NVIDIA GPU:** Volta architecture or better
- **CUDA Runtime Libraries:** 12 or 13

## Getting started / installation

- Make sure [CUDA runtime](https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64) is installed
- Download and unpack `cuopt-link-release-cu12.zip` or `cuopt-link-release-cu13-{x86_64,arm64}.zip` (for CUDA 12 and 13 respectively) from the [releases page](https://github.com/GAMS-dev/cuoptlink-builder/releases):
    - Unpack the contents of `cuopt-link-release-cu*-{x86_64,arm64}.zip` into your GAMS system directory. For GAMSPy, you can find out your system directory by running `gamspy show base`. So for example you can run `unzip -o cuopt-link-release-cu*.zip -d $(gamspy show base)`.
    - **Caution:** This will overwrite any existing `gamsconfig.yaml` file in that directory. The contained `gamsconfig.yaml` contains a `solverConfig` section to make cuOpt available to GAMS.

The neccessary files from the CUDA 12 or 13 runtime can also be downloaded as convenient archive `cu12-runtime-{x86_64,arm64}.zip` or `cu13-runtime-{x86_64,arm64}.zip` from the [releases page](https://github.com/GAMS-dev/cuoptlink-builder/releases).

## Test the setup

Get an example model and explicitly choose `cuopt` as `lp` or `mip` solver:
```
gamslib trnsport
gams trnsport lp cuopt
```

## Examples

- [examples/trnsport_cuopt.ipynb](examples/trnsport_cuopt.ipynb) for CUDA 12 on x86_64
- [examples/trnsport_cuopt.ipynb](examples/trnsport_cuopt_cu13.ipynb) for CUDA 13 on x86_64
