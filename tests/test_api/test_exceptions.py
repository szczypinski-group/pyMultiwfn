"""Tests for pymultiwfn.api.exceptions."""

import pytest

from pymultiwfn.api.exceptions import InvalidInputFileError, MultiwfnError


class TestMultiwfnError:
    """Tests for the MultiwfnError exception class."""

    def test_is_exception(self) -> None:
        assert issubclass(MultiwfnError, Exception)

    def test_raise_with_message(self) -> None:
        with pytest.raises(MultiwfnError, match="test error"):
            raise MultiwfnError("test error")

    def test_catch_and_inspect(self) -> None:
        try:
            raise MultiwfnError("caught")
        except MultiwfnError as e:
            assert str(e) == "caught"


class TestInvalidInputFileError:
    """Tests for the InvalidInputFileError exception class."""

    def test_is_multiwfn_error(self) -> None:
        assert issubclass(InvalidInputFileError, MultiwfnError)

    def test_message_contains_context(self) -> None:
        with pytest.raises(InvalidInputFileError, match="Unsupported"):
            raise InvalidInputFileError(
                input_path="bad.txt",
                extension=".txt",
                supported_extensions=(".wfn", ".fchk"),
            )
