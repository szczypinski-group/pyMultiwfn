"""
pyMultiwfn - A Python wrapper for Multiwfn batch automation.

Basic Usage
-----------
>>> from pyMultiwfn import MultiwfnJob, menu
>>> job = MultiwfnJob("molecule.wfn")
>>> job.add_menu_sequence(menu.hirshfeld_charge)
>>> result = job.run()
>>> charges = result.parse_charges()
"""

__version__ = "0.2.0"

# Config (includes MultiwfnError)
from .config import MultiwfnConfig, MultiwfnError

# Parsers
from .parsers import (
    BondOrderParser,
    ChargeParser,
    CriticalPointParser,
    OutputParser,
    ParserRegistry,
    SpectrumParser,
)

# Result
from .result import MultiwfnResult

# Job
from .job import MultiwfnJob, MultiwfnJobBuilder

# Menu submodule
from . import menu

__all__ = [
    '__version__',
    # Classes
    'MultiwfnJob',
    'MultiwfnJobBuilder',
    'MultiwfnResult',
    'MultiwfnConfig',
    # Exception
    'MultiwfnError',
    # Parsers
    'ParserRegistry',
    'OutputParser',
    'ChargeParser',
    'BondOrderParser',
    'CriticalPointParser',
    'SpectrumParser',
    # Submodules
    'menu',
]