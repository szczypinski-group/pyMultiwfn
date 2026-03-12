# Developper Guide

We actively welcome contributions from the community.
If you do decide to contribute please follow the guidelines below.

## Practical notes

### Development environment

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
4. At the moment, we follow numpy docstrings.

### Conventional commits

We try to follow the [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/).
In particular, please link to the issue that your commits are related to:

> docs: update use cases in README.md
>
> The use cases now follow the new code structure. Related to #24.
> Also included more detailed contribution guidelines.

## Adding a Function/Menu

If you are adding a new function/menu sequence, we ask that you ensure the following:

1. Add a new menu entry in `enums/menu.py`. This must be in the Enum format and must return an array of strings as input and be added in the correct category of calculation. A short description of what the sequence performs should also be added.
2. If relevant, include the new `pymultiwfn.Menu` Enum in the pre-defined analysis sets in `enums/analyses.py`.
3. Create a new `ParsedMultiwfnResult` class corresponding to the analysis results (if necessary) in `analysis/result.py`.
4. Include a new parsing regex and utility functions in `analysis/parsers.py`.
5. Ensure there exists a correct Menu <-> Parser mapping within `analysis.parsers.ParserRoute`.

## Use in other projects

When using our code for other projects, please respect the [Mozilla Public License 2.0](https://www.mozilla.org/en-US/MPL/2.0/).

In particular:

- You must provide the original/modified files licensed under MPL, alongside a list of changes from the original files.
- If using for closed source development, you must probide original source code of the files ("file-based copyleft") under MPL alongside your software.