"""Dealing with wavefunction files."""


class WfnFile:
    """Class for wavefunction files."""

    def __init__(self, file_path: str) -> None:
        """Initialize the WfnFile object."""
        self.file_path = file_path
