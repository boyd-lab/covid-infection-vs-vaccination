# Serology in SARS-CoV-2 vaccination vs infection

## Installation

```bash
# If using a custom conda environment, create and activate it:
conda create -n serology-env python=3.7.10
conda activate serology-env

# Install requirements
pip install --upgrade pip wheel
pip install -r requirements.txt

# Install local package
pip install -e .

# Install pre-commit
pre-commit install

# Run tests
make test

# Run lint
make lint

# Register this conda environment's Jupyter notebook kernel in base environment Jupyter:
ipython kernel install --user --name py37-serology-env

# Install fonts as TTF files to this directory:
python -c "import matplotlib; from pathlib import Path; print( (Path(matplotlib.matplotlib_fname()) / '../fonts/ttf').resolve() )"

# Rebuild font manager cache:
python -c "import matplotlib.font_manager; matplotlib.font_manager._rebuild()"
```

Review `config.py` configuration.

## Runbook

Extract serology measurements from raw data sheets, then plot:

```bash
# reset outputs
rm -r data/generated
rm -r out
mkdir -p data/generated
mkdir -p out

# run notebooks
./run_notebooks.sh \
    notebooks/load_data.infection_cohort1.ipynb \
    notebooks/load_data.pfizer_original.coronavirus_plate.ipynb \
    notebooks/load_data.pfizer_full_data.coronavirus_plate.ipynb \
    notebooks/load_data.pfizer_full_data.variant_plate.ipynb \
    notebooks/load_data.mongolia_vaccines.ipynb \
    notebooks/load_data.variant_infections.ipynb \
    notebooks/load_data.infection_cohort2.coronavirus_plate.ipynb \
    notebooks/load_data.infection_cohort2.variant_plate.ipynb \
    notebooks/load_data.concatenate.ipynb \
    notebooks/load_data.concatenate_coronavirus_plate_only.ipynb \
    notebooks/plot.ipynb \
    notebooks/summary.ipynb;
```

Review `summary.ipynb`.

## Development

This is based on the [datascience-cookiecutter-starter](https://github.com/maximz/datascience-cookiecutter-starter) template. `covid_serology/` is a library imported by the notebooks in `notebooks/` (mirrored to `notebooks_src/`).

### Jupytext mirroring

Using jupytext, every `.ipynb` notebook in `notebooks/` has a paired `.py` script in `notebooks_src/`. This makes diffs much cleaner (nbdiff/nbdime is too slow and finnicky to be practical) and allows for bulk refactoring.

When you edit and save the `.ipynb` notebook in Jupyter Lab, Jupytext updates the `.py` paired script automatically. And when you edit the `.py` script in a text editor, reloading the paired `.ipynb` notebook in Jupyter Lab will sync the updates to the notebook.

Sometimes we fall back to the command line to sync these notebook/script pairs: after editing either element of the pair, you can `git add` it to the git staging area and then `make lint-staged` to run jupytext-sync (and other pre-commit hooks) on all staged files. We tend to do this after bulk refactors of `.py` scripts (`make lint-staged` will update the paired notebooks without having to open each one in Jupyter Lab) or after auto-generating notebooks from other template notebooks with `run_notebook_to_new_file.sh` (`make lint-staged` will then generate or update the paired script.)

### Common commands

```bash
# lint all files
make lint
# or lint staged files only - most frequently used
make lint-staged
# or lint files that have changed since upstream only
make lint-diff

# run all tests
make test

# recreate all jupytext generated scripts
make regen-jupytext
```
