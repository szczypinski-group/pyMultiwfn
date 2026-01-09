"""
pyMultiwfn - A Python wrapper for Multiwfn batch automation.

Basic Usage
-----------
>>> from pyMultiwfn import MultiwfnJob, menu
>>> job = MultiwfnJob("molecule.wfn")
>>> job.add_menu_sequence(menu.hirshfeld_charge)
>>> result = job.run()
>>> charges = result.parse_charges()

>>> # Run ANY menu function by name
>>> from pyMultiwfn import run_analysis, list_functions
>>> list_functions()  # See all 190+ available functions
>>> result = run_analysis("mol.wfx", "igm_analysis")
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

# Analysis functions
from .analysis import (
    list_functions,
    get_function,
    run_analysis,
    run_menu_sequence,
    run_custom_sequence,
    quick_analysis,
    charge_analysis,
    bond_order_analysis,
    topology_analysis,
    nci_analysis,
)

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
    # Analysis
    'list_functions',
    'get_function',
    'run_analysis',
    'run_menu_sequence',
    'run_custom_sequence',
    'quick_analysis',
    'charge_analysis',
    'bond_order_analysis',
    'topology_analysis',
    'nci_analysis',
    # Submodules
    'menu',
]