# pyMultiwfn
A Python wrapper for automating [Multiwfn](http://sobereva.com/multiwfn/) batch calculations.
This project is currently in progress and is regularly maintained.
Information in the README might be inaccurate.   

All functions from the newest update are implemented with descriptive comments(multiwfn 3.8, 07/01/2026).   


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

## Contributing
We actively welcome contributions from the community. If you do decide to contribute please follow these rules:   
    1. Please ensure that any contributions pass unit testing before creating a pull request
        1.1 Please create unit tests in the appropriate testing file   
        1.2 For each function created, please write both a positive and negative unit test and any relevant edge cases   
    2. Please commit regularly   
    3. Ensure your contribution passes linting for easier code maintenance and consistency   
    4. Do not change the current architecture unless absolutely necessary. If you have any reccomendations for changes, please contact me at the email below.
## Contact
If you have any question or inquiries, please email me here: ypkh17@durham.ac.uk   
If you do send an email, please start the subject line with "pymultiwfn:" 

## Referencing
If you are using this package, please reference:
1. Tian Lu, J. Chem. Phys., 161, 082503 (2024) DOI: 10.1063/5.0216272 
