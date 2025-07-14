# GAMS solver link for NVIDIA cuOpt solver

[![Build cuOpt link for GAMS](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main.yml/badge.svg)](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main.yml)

This project builds and packages the GAMS link for the NVIDIA cuOpt solver.

## Requirements

- **Operating System:** Linux
- **GAMS:** Version 49 or newer.
- **GAMSPy:** Version 1.12.1 or newer
- **NVIDIA GPU:** Volta architecture or better
- **CUDA Runtime Libraries:** 12.0+

## Getting started / installation

- Make sure CUDA runtime is installed
- Download and unpack `cuopt-link-release.zip` from the [releases page](https://github.com/GAMS-dev/cuoptlink-builder/releases):
    - Unpack the contents of `cuopt-link-release.zip` into your GAMS system directory.
    - **Caution:** This will overwrite any existing `gamsconfig.yaml` file in that directory. The contained `gamsconfig.yaml` contains a `solverConfig` section to make cuOpt available to GAMS.

## Test the setup

Get an example model and explicitly choose `cuopt` as `lp` or `mip` solver:
```
gamsdist/gamslib trnsport
gamsdist/gams trnsport lp cuopt
```
