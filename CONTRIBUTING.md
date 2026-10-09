# Contributing

Thanks for helping improve the GAMS/cuOpt solver link! Bug reports, fixes and new cuOpt options are all welcome.

## Questions and bug reports

- Usage questions about GAMS or GAMSPy are best asked in the [GAMS forum](https://forum.gams.com/) or via support@gams.com.
- Bugs of the cuOpt link itself go into a [GitHub issue](https://github.com/GAMS-dev/cuoptlink-builder/issues/new/choose). Please include the versions and the solver log the bug report form asks for, and a small model that reproduces the problem if possible.
- Please report security issues privately, see [SECURITY.md](SECURITY.md).

## Building the link locally

The link (`gmscuopt.c`) is built against the GAMS C API from a GAMS distribution and the cuOpt libraries from NVIDIA's Python wheels:

```bash
python -m venv .venv
.venv/bin/pip install --extra-index-url=https://pypi.nvidia.com 'cuopt-cu13==26.8.*'
GAMSDIST=/path/to/gams ./build-link.sh
```

`build-link.sh` compiles `gmscuopt.out` for CUDA 13, sets its RPATH (needs `gcc` and `patchelf`), and copies it together with the cuOpt and CUDA libraries and the files from `assets/` into `$GAMSDIST`. Since that includes `gamsconfig_cuopt.yaml`, merge it into `gamsconfig.yaml` as described in the [README](README.md#manual-installation). The CI workflows in `.github/workflows` build the release archives for CUDA 12 and 13 on x86_64 and arm64.

New or changed cuOpt options go into `assets/optcuopt.def` and have to be passed on to cuOpt in `gmscuopt.c`.

## Testing

| Command | What it checks | GPU needed |
|---|---|---|
| `python3 tests/test-fetch-cuoptlink.py -v` | Unit tests of the installer script `fetch-cuoptlink.py` | no |
| `examples/models/regression_tests/run_tests.sh` | Self-checking models for duals, statuses, options, error handling and capability checks | yes |
| `python3 tests/test-gamslib-cuopt.py` | Solves the models in `tests/baseline.txt` from the GAMS model library and compares the objectives with CPLEX | yes |

The GPU tests use the GAMS system found in your `PATH`, with the link installed.

Without a local GPU, open `tests/colab-gpu-test.ipynb` in Google Colab with a GPU runtime: it builds the link from a branch the same way as CI, installs it into a fresh GAMS system and runs a `trnsport` smoke test and `tests/test-gamslib-cuopt.py`.

## Pull requests

- Keep each pull request focused on one change and describe how you tested it.
- If you fix a bug in the link, add or extend a model in `examples/models/regression_tests` that would have caught it.
- By contributing you agree that your contribution is licensed under the [Apache License 2.0](LICENSE).
