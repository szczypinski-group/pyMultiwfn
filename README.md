# pyMultiwfn

This project is currently in progress.
Information in the README might be inaccurate.

## Development

Make sure you have `uv` [installed](https://docs.astral.sh/uv/getting-started/installation/) on your system.
Then any `uv` command should create a fully functional local environment:

```
git clone git@github.com:szczypinski-group/pyMultiwfn.git
cd pyMultiwfn
uv sync
uvx pre-commit install
```

A Python wrapper for automating [Multiwfn](http://sobereva.com/multiwfn/) batch calculations.

Default linting settings and formatting settings (using [`ruff`](https://docs.astral.sh/ruff/)) have been created within `pyproject.toml` and will
be applied if the optional dependencies have been installed.
