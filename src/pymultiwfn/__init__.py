"""Python wrapper for Multiwfn batch automation.

Basic Usage
-----------
>>> from pymultiwfn import MultiwfnJob, menu
>>> job = MultiwfnJob("molecule.wfn")
>>> job.add_menu_sequence(menu.hirshfeld_charge)
>>> result = job.run()
>>> charges = result.parse_charges()
"""

__version__ = "0.2.0"

# Config (includes MultiwfnError)
# Menu submodule
from . import menu
from .config import MultiwfnConfig, MultiwfnError

# Job
from .job import MultiwfnJob, MultiwfnJobBuilder

# Parsers
from .parsers import (
    BondOrderParser,
    ChargeParser,
    CriticalPointParser,
    OutputParser,
    SpectrumParser,
)

# Result
from .result import MultiwfnResult

__all__ = [
    "__version__",
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
