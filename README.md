# GAMS solver link for NVIDIA cuOpt solver

[![Build cuOpt link for GAMS](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main-x86_64.yml/badge.svg)](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main-x86_64.yml) [![Build cuOpt link for GAMS](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main-arm64.yml/badge.svg)](https://github.com/GAMS-dev/cuoptlink-builder/actions/workflows/main-arm64.yml)

This project builds and packages the [GAMS](https://gams.com/) and [GAMSPy](https://gamspy.readthedocs.io/en/latest/index.html) solver link for the [NVIDIA cuOpt solver](https://github.com/NVIDIA/cuopt).

You can get more details and tips by reading the blog post ["GPU-Accelerated Optimization with GAMS and NVIDIA cuOpt"](https://www.gams.com/blog/2025/09/gpu-accelerated-optimization-with-gams-and-nvidia-cuopt/).

Supported model types are LP, MIP, RMIP, QCP, RMIQCP. QCP and RMIQCP models must be convex. Mixed-integer quadratic models (MIQCP) are not supported, since cuOpt's MIP solver only handles linear objectives and constraints.

## Quickstart

Install GAMSPy and the cuOpt link (including the CUDA runtime libraries) into GAMSPy's GAMS system directory:

```bash
pip install gamspy
curl -O https://raw.githubusercontent.com/GAMS-dev/cuoptlink-builder/main/fetch-cuoptlink.py
python fetch-cuoptlink.py install -g "$(gamspy show base)" --cuda-runtime
```

Then pick `cuopt` as solver:

```python
import gamspy as gp

gp.set_options({"SOLVER_VALIDATION": 0})  # cuOpt is added to GAMSPy by hand

m = gp.Container()
x = gp.Variable(m, "x", type="positive")
y = gp.Variable(m, "y", type="positive")
e1 = gp.Equation(m, "e1", definition=x + 2 * y >= 2)
e2 = gp.Equation(m, "e2", definition=3 * x + y >= 3)

model = gp.Model(m, "demo", equations=[e1, e2], problem="LP", sense="min", objective=x + y)
model.solve(solver="cuopt")
print(model.objective_value)  # 1.4
```

For a standalone GAMS system pass its directory to `-g` instead and see [Test the setup](#test-the-setup). To try it without a local setup, open one of the [example notebooks](#notebooks) in Google Colab.

## Performance

cuOpt's GPU-based PDLP method pays off on large LPs. In NVIDIA's [benchmark on Mittelmann's LP test set](https://developer.nvidia.com/blog/accelerate-large-linear-programming-problems-with-nvidia-cuopt/) (October 2024, H100 SXM GPU, no presolve), cuOpt was faster than a state-of-the-art CPU LP solver on 60% of the instances, more than 10x faster on 20%, and up to 5000x faster on a large multi-commodity flow instance, while 8 of the 49 public instances hit the one-hour time limit. On small models the GPU overhead usually dominates, so compare on your own instances. The default `method 0` (concurrent) runs PDLP, dual simplex and barrier in parallel.

## Requirements

- **Operating System:** Linux, Windows 11 through WSL2
- **CPU architecture:** x86_64, arm64
- **GAMS:** Version 54 or newer
- **GAMSPy:** Version 1.12.1 or newer
- **NVIDIA GPU:** Volta architecture or better
- **CUDA Runtime Libraries:** 12 or 13

## Installation using `fetch-cuoptlink.py`

You can automatically download, install, test, and manage the cuOpt solver link using the provided `fetch-cuoptlink.py` script.

__Quickstart:__ The script has no external dependencies, so you can download and run it directly:

```bash
curl -O https://raw.githubusercontent.com/GAMS-dev/cuoptlink-builder/main/fetch-cuoptlink.py
python fetch-cuoptlink.py
```

or as a one-liner with [uv](https://docs.astral.sh/uv/):

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

> **Note:** Successful installations automatically verify the solver link by solving a small embedded LP model with `solver=cuopt`.

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

- [examples/trnsport_cuopt.ipynb](examples/trnsport_cuopt.ipynb) for CUDA 12 [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GAMS-dev/cuoptlink-builder/blob/main/examples/trnsport_cuopt.ipynb)
- [examples/trnsport_cuopt_cu13.ipynb](examples/trnsport_cuopt_cu13.ipynb) for CUDA 13 [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GAMS-dev/cuoptlink-builder/blob/main/examples/trnsport_cuopt_cu13.ipynb)

In Google Colab, select a GPU runtime first (*Runtime → Change runtime type*). If unsure, start with the CUDA 12 notebook, which also works with older NVIDIA drivers.

### GAMS models

Various GAMS models can be found in subfolder `examples/models` and are used to verify the solver link.
### Regression tests

The self-checking models in `examples/models/regression_tests` cover dual signs and reduced costs, QP/QCQP marginals, RMIQCP, option handling, error reporting, LP limit points, GMO handling (e.g. `=N=` rows, `requestMarginals=2`), rejection of unsupported features (SOS, semi-integer, MIQCP) and solve/model status mapping. Each model aborts if a result deviates from the reference values (obtained with CPLEX). Run them all against the GAMS system found in your `PATH` (it needs a GPU and the installed solver link):

```
examples/models/regression_tests/run_tests.sh
```

The script prints `[PASS]` or `[FAIL]` per model and keeps the listing and log file of failed models for inspection. Its exit code is the number of failed models.
