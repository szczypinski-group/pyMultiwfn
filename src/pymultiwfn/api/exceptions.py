"""Exceptions for pymultiwfn interface."""


class MultiwfnError(Exception):
    """Base exception for pyMultiwfn errors."""

    pass


class InvalidInputFileError(MultiwfnError):
    """Raised when an input file extension is unsupported by Multiwfn."""

    def __init__(
        self,
        input_path: str,
        extension: str,
        supported_extensions: tuple[str, ...],
    ) -> None:
        supported = ", ".join(supported_extensions)
        message = (
            f"Unsupported Multiwfn input file: '{input_path}' (extension "
            f"'{extension or '<none>'}'). Supported extensions include: "
            f"{supported}. See Multiwfn manual input-format section."
        )
        super().__init__(message)
