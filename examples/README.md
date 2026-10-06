## GAMSPy cuOpt integration example notebook

### Run in Google Colab
- CUDA 12: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GAMS-dev/cuoptlink-builder/blob/main/examples/trnsport_cuopt.ipynb)
- CUDA 13: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GAMS-dev/cuoptlink-builder/blob/main/examples/trnsport_cuopt_cu13.ipynb)
- Select a GPU runtime first (*Runtime → Change runtime type*), then run all cells.

### How to run the notebook locally
- Setup and activate a Python environment e.g. with `python -m venv .venv` and `source .venv/bin/activate`
- Install requirements `pip install -r requirements.txt`
- Run this notebook from a GUI or update its cells in the command line via `jupyter nbconvert --to notebook --execute --inplace trnsport_cuopt.ipynb`
- The first cell takes care of installing cuOpt runtime, GAMSPy, and the solver link. Skips most parts when it's already there.
