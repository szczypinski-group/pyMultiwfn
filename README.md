# pyMultiwfn

[![DOI](https://zenodo.org/badge/1098975665.svg)](https://doi.org/10.5281/zenodo.23740372)

A Python wrapper for automating [Multiwfn 3.8 (07/01/2026)](http://sobereva.com/multiwfn/) batch calculations.
This project is currently in progress and is regularly maintained.
Documentation for `pymultiwfn` can be found [here](https://szczypinski-group.github.io/pyMultiwfn/).

## Installation

You can now install directly from [PyPI](https://pypi.org/project/pymultiwfn/) (via `pip` or your manager of choice, e.g., `uv`):

```
pip install pymultiwfn
```

## Running calculations

The most simple entry-point to analysis is through the `MultiwfnAnalysis` class.
Different **Multiwfn** menu entries can then by added as `analyses` and a `MultiwfnJob` will be created when the analysis is `run()`.
There is also a possibility to create `MultiwfnJob` objects directly, which could be useful when executing custom menu sequences for tailored applications.

The `Multiwfn` object currently holds a path to the **Multiwfn** executable.
If it is not specified, the bundled version is used.


```python
from pymultiwfn import MultiwfnAnalysis, Menu

analysis = MultiwfnAnalysis("benzene.wfn", analyses=[
    Menu.HIRSHFELD_CHARGE,
    Menu.MAYER_BOND_ORDER,
])

analysis.run()
print(analysis.results[0].to_dict())
```

The `MultiwfnResult` dataclasses hold different types of parsed results, which can easily be serialised into JSON.

## Tutorials

The `tutorials` folder holds three [marimo](https://marimo.io/) notebooks, one per case study.
Each runs a small batch of analyses on a real input file and walks through the results.
They are complete with a set of expected results to compare against.

| Tutorial | Input file | Analyses | Run time (approx.) |
|---|---|---|---|
| `chalcogen_dimer/chalcogen_bonding_tutorial.py` | `coord.molden` | `NCI_ANALYSIS` | 20 s |
| `kraken_ligand/kraken_ligand_tutorial.py` | `kraken_ligand.molden` | `VIEW_STRUCTURE`, `SURFACE_AREA_VOLUME`, `QMSA_ALIE`, `QMSA_ESP` | 2 min |
| `optoelectronic_materials/optoelectronic_screening_tutorial.py` | `optoelectronic.molden` | `VIEW_STRUCTURE` | a few seconds |

### Running a tutorial

The tutorials live in this repository, not in the PyPI package, and marimo is not installed with `pymultiwfn`.
From a clone of the repository (see [Development](#development)), open a tutorial with:

```
uv run --with marimo marimo edit tutorials/chalcogen_dimer/chalcogen_bonding_tutorial.py
```

This opens the notebook in your browser and runs it.
Use `marimo run` instead of `marimo edit` to view it as a read-only app.
If you installed `pymultiwfn` with `pip`, install marimo alongside it and open the notebook the same way:

```
pip install marimo
marimo edit tutorials/chalcogen_dimer/chalcogen_bonding_tutorial.py
```

You can launch a tutorial from any folder: each notebook finds its input file and writes its results relative to its own location.

### Tutorial results

Each tutorial writes its output to `tutorial_results/<input_file_name>.output/` inside its own folder: the captured Multiwfn output for every analysis, any files Multiwfn generated (such as cube files), and one JSON summary that the notebook reads its results from.
The copies committed to the repository are the expected results.

Analyses are cached: when a tutorial finds its results already in `tutorial_results/`, it reuses them instead of running Multiwfn again.
To recalculate and compare against the expected results, move that folder aside and run the tutorial again.

## Contributing

We actively welcome contributions from the community.
If you do decide to contribute please follow the guidelines below.
When using our code for other projects, please respect the [Mozilla Public License 2.0](https://www.mozilla.org/en-US/MPL/2.0/).
In particular:

- You must provide the original/modified files licensed under MPL, alongside a list of changes from the original files.
- If using for closed source development, you must probide original source code of the files ("file-based copyleft") under MPL alongside your software.


### Development

Make sure you have `uv` [installed](https://docs.astral.sh/uv/getting-started/installation/) on your system.
Then any `uv` command should create a fully functional local environment:

```
git clone git@github.com:szczypinski-group/pyMultiwfn.git
cd pyMultiwfn
uv sync
uvx pre-commit install
```

### Linting

Default linting settings and formatting settings (using [`ruff`](https://docs.astral.sh/ruff/)) have been created within `pyproject.toml` and will
be applied if the optional dependencies have been installed.

### Unit testing

We try to include positive and negative unit tests for each new function.
Ensure the tests pass before submitting a pull request:

```
uv run pytest
```

### Code standards

1. Use explicit imports wherever possible.
2. For class and function definition/calls, split arguments into multiple lines.
3. Always call functions with keyword arguments.

### Conventional commits

We try to follow the [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/).
In particular, please link to the issue that your commits are related to:

> docs: update use cases in README.md
>
> The use cases now follow the new code structure. Related to #24.
> Also included more detailed contribution guidelines.

## Contact

If you have any question or inquiries, please create a [new issue](https://github.com/szczypinski-group/pyMultiwfn/issues).

## Referencing

If you are using this package, please reference the original Multiwfn package:

1. Tian Lu, *J. Chem. Phys.*, 161, 082503 (2024) DOI: [10.1063/5.0216272](https://doi.org/10.1063/5.0216272).
