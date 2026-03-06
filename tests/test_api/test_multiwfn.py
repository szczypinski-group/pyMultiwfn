"""Tests for pymultiwfn.api.multiwfn — Multiwfn executable config."""

from pathlib import Path

import pytest

from pymultiwfn.api.exceptions import MultiwfnError
from pymultiwfn.api.multiwfn import Multiwfn


class TestMultiwfnInit:
    """Tests for Multiwfn.__init__() and _parse_exe()."""

    def test_with_valid_exe(self, mock_executable: Path) -> None:
        mwfn = Multiwfn(exe_path=mock_executable)
        assert mwfn.exe_path == mock_executable.resolve()

    def test_with_nonexistent_exe(self, temp_dir: Path) -> None:
        with pytest.raises(MultiwfnError, match="not found"):
            Multiwfn(exe_path=temp_dir / "nonexistent")

    def test_default_raises_without_bundled(self) -> None:
        """Default init will raise if no bundled exe exists."""
        with pytest.raises(MultiwfnError):
            Multiwfn()

    def test_str(self, mock_executable: Path) -> None:
        mwfn = Multiwfn(exe_path=mock_executable)
        assert str(mock_executable.resolve()) in str(mwfn)

    def test_repr(self, mock_executable: Path) -> None:
        mwfn = Multiwfn(exe_path=mock_executable)
        assert "pymultiwfn.Multiwfn" in repr(mwfn)
        assert "exe_path=" in repr(mwfn)
