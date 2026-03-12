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
from pymultiwfn import MultiwfnAnalysis, Menu

analysis = MultiwfnAnalysis("benzene.wfn", analyses=[
    Menu.HIRSHFELD_CHARGE,
    Menu.MAYER_BOND_ORDER,
])

analysis.run()
print(analysis.results[0].to_dict())
```

The `MultiwfnResult` dataclasses hold different types of parsed results, which can easily be serialised into JSON.

## Contact

If you have any question or inquiries, please create a [new issue](https://github.com/szczypinski-group/pyMultiwfn/issues).

## Referencing

If you are using this package, please reference the original Multiwfn package:

1. Tian Lu, *J. Chem. Phys.*, 161, 082503 (2024) DOI: [10.1063/5.0216272](https://doi.org/10.1063/5.0216272).
