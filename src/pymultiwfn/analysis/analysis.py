"""Module to plan and execute Multiwfn analyses."""

from pathlib import Path

from pymultiwfn.enums.menu import Menu


class MultiwfnAnalysis:
    """Base class for Multiwfn analyses.

    This class serves as a template for specific analyses. Each analysis type
    should inherit from this base class and implement the get_menu_sequence()
    method to provide the appropriate menu commands for that analysis.

    Attributes
    ----------
    input_file : str or Path
        Path to the wavefunction file to be analyzed.

    Examples
    --------
    >>> # Example of a specific analysis inheriting from MultiwfnAnalysis
    >>> class ElectronDensityAnalysis(MultiwfnAnalysis):
    ...     def get_menu_sequence(self) -> list[str]:
    ...         return ["1", "2", "3"]

    """

    def __init__(
        self,
        input_file: str | Path,
        analyses: list[Menu] | None = None,
    ) -> None:
        self.input_file = input_file
        self.analyses = analyses if analyses is not None else []
