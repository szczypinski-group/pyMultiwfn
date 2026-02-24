"""Python wrapper for Multiwfn batch automation.

Usage
-----
>>> import pymultiwfn
>>> mol = pymultiwfn.load("molecule.wfn")
>>> charges = mol.run_charges()
>>> result = mol.run(pymultiwfn.Menu.HIRSHFELD_CHARGE)
"""

__version__ = "0.1.1"

from pymultiwfn.analysis import MultiwfnAnalysis
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

# Alias for convenience
Analysis = MultiwfnAnalysis

#TODO(fs): maybe improve the entry point(multiwfnjob instead of multiwfnanalysis?)
def load(
    filepath: str,
    config: MultiwfnConfig | None = None,
) -> MultiwfnAnalysis:
    """Load a wavefunction file for analysis.

    Parameters
    ----------
    filepath : str or Path
        Path to wavefunction file (.wfn, .fchk, .molden, etc.)
    config : MultiwfnConfig, optional
        Configuration for Multiwfn execution

    Returns
    -------
    MultiwfnAnalysis
        Analysis object with methods to run calculations

    Examples
    --------
    >>> import pymultiwfn
    >>> mol = pymultiwfn.load("molecule.wfn")
    >>> charges = mol.run_charges()
    >>> result = mol.run(pymultiwfn.Menu.HIRSHFELD_CHARGE)
    """
    return MultiwfnAnalysis(filepath, config)


__all__ = [
    # Version
    "__version__",
    # Main entry point
    "load",
    # Analysis class
    "MultiwfnAnalysis",
    "Analysis",
    # Core classes
    "Menu",
    "MultiwfnJob",
    "MultiwfnJobBuilder",
    "MultiwfnResult",
    "MultiwfnConfig",
    "MultiwfnError",
    # Parsers
    "OutputParser",
    "ChargeParser",
    "BondOrderParser",
    "CriticalPointParser",
    "SpectrumParser",
]
