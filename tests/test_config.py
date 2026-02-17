"""Tests for config.py - MultiwfnConfig and MultiwfnError."""

from pathlib import Path

import pytest

from pymultiwfn import MultiwfnConfig, MultiwfnError


class TestMultiwfnError:
    """Tests for MultiwfnError exception."""

    def test_error_is_exception(self) -> None:
        """MultiwfnError should be an Exception."""
        assert issubclass(MultiwfnError, Exception)

    def test_error_can_be_raised(self) -> None:
        """MultiwfnError can be raised with message."""
        with pytest.raises(MultiwfnError, match="test error"):
            raise MultiwfnError("test error")

    def test_error_can_be_caught(self) -> None:
        """MultiwfnError can be caught."""
        try:
            raise MultiwfnError("caught")
        except MultiwfnError as e:
            assert str(e) == "caught"


class TestMultiwfnConfig:
    """Tests for MultiwfnConfig class."""

    def test_default_config(self) -> None:
        """Default config should have sensible defaults."""
        config = MultiwfnConfig()

        assert config.exe_path is None
        assert config.working_dir == Path.cwd()
        assert config.timeout is None
        assert config.verbose is False

    def test_custom_config(self, temp_dir: Path) -> None:
        """Config accepts custom values."""
        config = MultiwfnConfig(
            exe_path=Path("/custom/path"),
            working_dir=temp_dir,
            timeout=300,
            verbose=True,
        )

        assert config.exe_path == Path("/custom/path")
        assert config.working_dir == temp_dir
        assert config.timeout == 300
        assert config.verbose is True

    def test_explicit_exe_path_exists(self, mock_executable: Path) -> None:
        """Config uses explicit exe_path when it exists."""
        config = MultiwfnConfig(exe_path=mock_executable)
        assert config.executable == mock_executable

    def test_explicit_exe_path_not_found(self, temp_dir: Path) -> None:
        """Config raises error when explicit exe_path doesn't exist."""
        config = MultiwfnConfig(exe_path=temp_dir / "nonexistent")

        with pytest.raises(MultiwfnError, match="not found"):
            _ = config.executable

    def test_executable_cached(self, mock_executable: Path) -> None:
        """Executable resolution is cached."""
        config = MultiwfnConfig(exe_path=mock_executable)

        exe1 = config.executable
        exe2 = config.executable

        assert exe1 is exe2

    def test_executable_not_found_anywhere(
        self, temp_dir: Path, monkeypatch: pytest
    ) -> None:
        """Config raises error when executable not found anywhere."""
        # Ensure shutil.which returns None
        monkeypatch.setattr("shutil.which", lambda x: None)
        # Ensure Path.exists returns False for all candidates
        original_exists = Path.exists

        def fake_exists(self: Path) -> bool:
            if "Multiwfn" in str(self):
                return False
            return original_exists(self)

        monkeypatch.setattr(Path, "exists", fake_exists)

        config = MultiwfnConfig(working_dir=temp_dir)

        with pytest.raises(MultiwfnError, match="not found"):
            _ = config.executable


class TestMultiwfnConfigWithRealExecutable:
    """Tests using real Multiwfn executable."""

    @pytest.mark.requires_multiwfn
    def test_finds_real_executable(self, real_executable: Path) -> None:
        """Config can find real Multiwfn executable."""
        config = MultiwfnConfig(exe_path=real_executable)
        assert config.executable.exists()
        assert "Multiwfn" in config.executable.name

    @pytest.mark.requires_multiwfn
    def test_auto_discovery(self, project_root: Path) -> None:
        """Config auto-discovers executable in bin/ folder."""
        # This test relies on the walk-up logic finding bin/
        config = MultiwfnConfig()
        try:
            exe = config.executable
            assert exe.exists()
        except MultiwfnError:
            pytest.skip("Multiwfn not in expected location")
