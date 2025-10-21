# GAMS solver link for NVIDIA cuOpt solver

[![Build cuOpt link for GAMS](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main.yml/badge.svg)](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main.yml)

This project builds and packages the [GAMS](https://gams.com/) and [GAMSPy](https://gamspy.readthedocs.io/en/latest/index.html) solver link for the [NVIDIA cuOpt solver](https://github.com/NVIDIA/cuopt).

You can get more details and tips by reading the blog post ["GPU-Accelerated Optimization with GAMS and NVIDIA cuOpt"](https://www.gams.com/blog/2025/09/gpu-accelerated-optimization-with-gams-and-nvidia-cuopt/).

## Requirements

- **Operating System:** Linux, Windows through WSL2
- **GAMS:** Version 49 or newer.
- **GAMSPy:** Version 1.12.1 or newer
- **NVIDIA GPU:** Volta architecture or better
- **CUDA Runtime Libraries:** 12.0+

## Getting started / installation

- Make sure [CUDA runtime](https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64) is installed
- Download and unpack `cuopt-link-release.zip` from the [releases page](https://github.com/GAMS-dev/cuoptlink-builder/releases):
    - Unpack the contents of `cuopt-link-release.zip` into your GAMS system directory. For GAMSPy, you can find out your system directory by running `gamspy show base`. So for example you can run `unzip -o cuopt-link-release.zip -d $(gamspy show base)`.
    - **Caution:** This will overwrite any existing `gamsconfig.yaml` file in that directory. The contained `gamsconfig.yaml` contains a `solverConfig` section to make cuOpt available to GAMS.

More specifically, the files from the CUDA runtime needed are
```
libcublas.so.12
libcublasLt.so.12
libcudss.so.0
libcurand.so.10
libcusolver.so.11
libnvJitLink.so.12
```
and can be installed e.g. via `pip install --extra-index-url=https://pypi.nvidia.com cuopt-cu12==25.5.* nvidia-cuda-runtime-cu12==12.8.* nvidia-nvjitlink-cu12` into a Python environment or downloaded as archive `cu12-runtime.zip` from the [releases page](https://github.com/GAMS-dev/cuoptlink-builder/releases).

## Test the setup

Get an example model and explicitly choose `cuopt` as `lp` or `mip` solver:
```
gamslib trnsport
gams trnsport lp cuopt
```

## Examples

- [examples/trnsport_cuopt.ipynb](examples/trnsport_cuopt.ipynb)
