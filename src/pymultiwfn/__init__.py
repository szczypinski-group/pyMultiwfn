"""Python wrapper for Multiwfn batch automation.

Usage
-----
>>> import pymultiwfn
>>> mol = pymultiwfn.load("molecule.wfn")
>>> charges = mol.run_charges()
>>> result = mol.run(pymultiwfn.Menu.HIRSHFELD_CHARGE)
"""

from importlib.metadata import version

from pymultiwfn.analysis.parsers import (
    BondOrderParser,
    ChargeParser,
    CriticalPointParser,
    OutputParser,
    SpectrumParser,
)
from pymultiwfn.analysis.result import MultiwfnResult, ResultStore
from pymultiwfn.api.exceptions import InvalidInputFileError, MultiwfnError
from pymultiwfn.api.job import MultiwfnJob
from pymultiwfn.api.multiwfn import Multiwfn
from pymultiwfn.enums.menu import Menu

__version__ = version("pymultiwfn")
__all__ = [
    # Version
    "__version__",
    # Main entry point
    "load",
    # Analysis class
    "MultiwfnAnalysis",
    # Core classes
    "Menu",
    "MultiwfnJob",
    "MultiwfnResult",
    "Multiwfn",
    "MultiwfnError",
    "InvalidInputFileError",
    # Parsers
    "OutputParser",
    "ChargeParser",
    "BondOrderParser",
    "CriticalPointParser",
    "SpectrumParser",
    # Storage
    "ResultStore",
]


def __getattr__(name: str) -> object:
    """Lazily resolve heavy imports to keep package import robust."""
    if name == "MultiwfnAnalysis":
        from pymultiwfn.analysis.analysis import MultiwfnAnalysis

        return MultiwfnAnalysis
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
