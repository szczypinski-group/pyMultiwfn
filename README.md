# pyMultiwfn

A Python wrapper for automating [Multiwfn](http://sobereva.com/multiwfn/) batch calculations.
This project is currently in progress and is regularly maintained.
Information in the README might be inaccurate.

All functions from the newest update are implemented with descriptive comments(multiwfn 3.8, 07/01/2026).

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
from pathlib import Path
from pymultiwfn import MultiwfnAnalysis, Menu

mwfn = Multiwfn(exe_path=Path("/opt/multiwfn/Multiwfn"))

analysis = MultiwfnAnalysis("benzene.wfn", analyses=[
    Menu.HIRSHFELD_CHARGE,
    Menu.MAYER_BOND_ORDER,
])
analysis.run(work_dir=Path("./output"))

print(analysis.results)
```

The `MultiwfnResult` dataclasses hold different types of parsed results, which can easily be serialised into JSON.

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
