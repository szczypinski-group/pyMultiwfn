"""Python wrapper for Multiwfn batch automation.

Usage
-----
>>> import pymultiwfn
>>> mol = pymultiwfn.load("molecule.wfn")
>>> charges = mol.run_charges()
>>> result = mol.run(pymultiwfn.Menu.HIRSHFELD_CHARGE)
"""

__version__ = "0.5.2"

from pymultiwfn.analysis.analysis import MultiwfnAnalysis
from pymultiwfn.analysis.parsers import (
    BondOrderParser,
    ChargeParser,
    CriticalPointParser,
    OutputParser,
    SpectrumParser,
)
from pymultiwfn.analysis.result import MultiwfnResult, ResultStore
from pymultiwfn.api.exceptions import MultiwfnError
from pymultiwfn.api.job import MultiwfnJob
from pymultiwfn.api.multiwfn import Multiwfn
from pymultiwfn.enums.menu import Menu

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
    # Parsers
    "OutputParser",
    "ChargeParser",
    "BondOrderParser",
    "CriticalPointParser",
    "SpectrumParser",
    # Storage
    "ResultStore",
]
