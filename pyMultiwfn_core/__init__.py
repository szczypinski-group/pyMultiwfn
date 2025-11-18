"""
pyMultiwfn - Python wrapper for Multiwfn wavefunction analysis
================================================================

A comprehensive Python library for automating Multiwfn batch calculations.

Basic Usage
-----------
>>> from pyMultiwfn import MultiwfnJob
>>> from pyMultiwfn import menus
>>>
>>> # Create a job
>>> job = MultiwfnJob("molecule.molden")
>>>
>>> # Add analyses
>>> job.add_menu_sequence(menus.hirshfeld_charge)
>>> job.add_menu_sequence(menus.mayer_bond_order)
>>>
>>> # Run and get results
>>> job.run()
>>> charges = job.parse_charges()
>>> bonds = job.parse_bond_orders()

Quick Analysis
--------------
>>> from pyMultiwfn import quick_analysis
>>>
>>> # Run a quick charge analysis
>>> results = quick_analysis("molecule.molden", "charges", method="hirshfeld")
>>> print(results['charges'])
"""

from .core import (
    MultiwfnJob,
    MultiwfnError,
    MultiwfnExecutionError,
    MultiwfnOutputParser,
    generate_batch_file,
    run_multiwfn,
    quick_analysis,
    __version__,
)

# Import menu functions for convenience
from . import menus

__all__ = [
    # Main classes
    'MultiwfnJob',
    'MultiwfnError',
    'MultiwfnExecutionError',
    'MultiwfnOutputParser',
    
    # Legacy functions
    'generate_batch_file',
    'run_multiwfn',
    
    # Convenience
    'quick_analysis',
    'menus',
    
    # Version
    '__version__',
]

# Package metadata
__author__ = "Your Name"
__email__ = "your.email@example.com"
__license__ = "MIT"