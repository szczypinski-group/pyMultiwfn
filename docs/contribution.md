# Developper Guide
## Adding a Function/Menu
If you are adding a new function/menu sequence, we ask that you ensure the following:
1. Add a new menu entry in enum/menu.py. This must be in the enum format and must return an array of strings as input and be added in the correct category of calculation. A short description of what the sequence performs should also be added.   
2. Include the new Menu enum to analysis.   
3. Write/include a new parsing function in parsers.py if the Multiwfn output is non-standard.   
4. Ensure Menu <-> Parser mapping   

## Contributing
We actively welcome contributions from the community. If you do decide to contribute please follow these rules:
1. Please ensure that any contributions pass unit testing before creating a pull request.   
2. Please create unit tests in the appropriate testing file.   
3. For each function created, please write both a positive and negative unit test and any relevant edge cases.   
4. Please commit regularly.   
5. Ensure your contribution passes linting for easier code maintenance and consistency.   
6. Do not change the current architecture unless absolutely necessary. If you have any reccomendations for changes, please contact me here: ypkh17@durham.ac.uk.   

## Code standards
1. Use explicit imports wherever possible.
2. For class and function definition/calls, split arguments into multiple lines.
3. Always call functions with keyword arguments.

## Further Points
1. Use Numpy docstrings and explicit typing
