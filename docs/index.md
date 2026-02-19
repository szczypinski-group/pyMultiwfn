# pyMultiwfn
A Python wrapper for automating [Multiwfn](http://sobereva.com/multiwfn/) batch calculations.
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

## Loading Files
Appropriate file types can be loaded in the package using the load() function.

```
import pymultiwfn
mol = pymultiwfn.load("molecule.wfn")
```

## Running calculations
Calculations can be run in a number of different ways:  
    1. running the calculations individually by calling the desired function relating to a specific Multiwfn input  
    2. creating a custom menu, composed of the desired calculations to be run  
    3. Running calculation suites(all charges, all topologies, etc...)  
    4. Running all the avaialble calculations integrated in the package   
```
from pymultiwfn import Analysis, Menu
analysis = Analysis("molecule.wfn")
1. Single analysis
result = analysis.run(Menu.HIRSHFELD_CHARGE)
charges = result.parse_charges()

2. Multiple analyses in one session
menus = [Menu.HIRSHFELD_CHARGE, Menu.MAYER_BOND_ORDER]
result = analysis.run_multiple(menus)

3. All charges
results = analysis.run_charges()
hirshfeld = results["hirshfeld_charge"].parse_charges()

4. Everything
all_results = analysis.run_all()
```

## Linting
Default linting settings and formatting settings (using [`ruff`](https://docs.astral.sh/ruff/)) have been created within `pyproject.toml` and will
be applied if the optional dependencies have been installed.

## Referencing
If you are using this package, pelase reference:
1. Tian Lu, J. Chem. Phys., 161, 082503 (2024) DOI: 10.1063/5.0216272 
