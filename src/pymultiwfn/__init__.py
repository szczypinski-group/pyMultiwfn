"""Python wrapper for Multiwfn batch automation.

Basic Usage
-----------
>>> from pymultiwfn import MultiwfnJob, menu
>>> job = MultiwfnJob("molecule.wfn")
>>> job.add_menu_sequence(menu.hirshfeld_charge)
>>> result = job.run()
>>> charges = result.parse_charges()
"""

__version__ = "0.1.0"

from pymultiwfn.config import MultiwfnConfig, MultiwfnError
from pymultiwfn.job import MultiwfnJob, MultiwfnJobBuilder
from pymultiwfn.menu import Menu
from pymultiwfn.parsers import (
    BondOrderParser,
    ChargeParser,
    CriticalPointParser,
    OutputParser,
    SpectrumParser,
)
from pymultiwfn.result import MultiwfnResult

__all__ = [
    "__version__",
    "Menu",
    # Classes
    "MultiwfnJob",
    "MultiwfnJobBuilder",
    "MultiwfnResult",
    "MultiwfnConfig",
    # Exception
    "MultiwfnError",
    # Parsers
    "OutputParser",
    "ChargeParser",
    "BondOrderParser",
    "CriticalPointParser",
    "SpectrumParser",
    # Submodules
    "menu",
]
