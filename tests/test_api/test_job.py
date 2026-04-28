"""Tests for pymultiwfn.api.job — MultiwfnJob."""

from pathlib import Path

import pytest

from pymultiwfn.api.exceptions import InvalidInputFileError, MultiwfnError
from pymultiwfn.api.job import MultiwfnJob
from pymultiwfn.api.multiwfn import Multiwfn
from pymultiwfn.enums.menu import Menu


@pytest.fixture
def multiwfn(mock_executable: Path) -> Multiwfn:
    return Multiwfn(exe_path=mock_executable)


class TestJobInit:
    """Tests for the __init__ method of MultiwfnJob, including edge cases."""

    def test_create_with_valid_file(
        self, mock_wfn_file: Path, multiwfn: Multiwfn
    ) -> None:
        job = MultiwfnJob(
            input_file=mock_wfn_file,
            analysis=Menu.HIRSHFELD_CHARGE,
            multiwfn=multiwfn,
        )
        assert job.input_file == mock_wfn_file.resolve()
        assert not job.executed

    def test_file_not_found(self, multiwfn: Multiwfn) -> None:
        with pytest.raises(FileNotFoundError):
            MultiwfnJob(
                input_file="/nonexistent/file.wfn",
                analysis=Menu.HIRSHFELD_CHARGE,
                multiwfn=multiwfn,
            )

    def test_unsupported_input_extension_raises(
        self, temp_dir: Path, multiwfn: Multiwfn
    ) -> None:
        unsupported = temp_dir / "bad_input.txt"
        unsupported.write_text("dummy")
        with pytest.raises(InvalidInputFileError, match="Unsupported"):
            MultiwfnJob(
                input_file=unsupported,
                analysis=Menu.HIRSHFELD_CHARGE,
                multiwfn=multiwfn,
            )

    def test_none_analysis_creates_empty_commands(
        self, mock_wfn_file: Path, multiwfn: Multiwfn
    ) -> None:
        job = MultiwfnJob(
            input_file=mock_wfn_file,
            analysis=None,
            multiwfn=multiwfn,
        )
        assert job.commands == []

    def test_default_multiwfn_raises(self, mock_wfn_file: Path) -> None:
        """Default Multiwfn() either resolves bundled exe or raises."""
        try:
            job = MultiwfnJob(input_file=mock_wfn_file, analysis=None)
            assert isinstance(job, MultiwfnJob)
        except MultiwfnError:
            pass

    def test_invalid_timeout(
        self, mock_wfn_file: Path, multiwfn: Multiwfn
    ) -> None:
        with pytest.raises(ValueError, match="positive"):
            MultiwfnJob(
                input_file=mock_wfn_file,
                analysis=None,
                multiwfn=multiwfn,
                timeout=-1,
            )

    def test_zero_timeout(
        self, mock_wfn_file: Path, multiwfn: Multiwfn
    ) -> None:
        with pytest.raises(ValueError, match="positive"):
            MultiwfnJob(
                input_file=mock_wfn_file,
                analysis=None,
                multiwfn=multiwfn,
                timeout=0,
            )


class TestJobProperties:
    """Tests for properties of MultiwfnJob."""

    def test_commands_from_menu(
        self, mock_wfn_file: Path, multiwfn: Multiwfn
    ) -> None:
        job = MultiwfnJob(
            input_file=mock_wfn_file,
            analysis=Menu.HIRSHFELD_CHARGE,
            multiwfn=multiwfn,
        )
        cmds = job.commands
        assert "7" in cmds  # charges menu
        assert isinstance(cmds, list)

    def test_commands_is_copy(
        self, mock_wfn_file: Path, multiwfn: Multiwfn
    ) -> None:
        job = MultiwfnJob(
            input_file=mock_wfn_file,
            analysis=None,
            multiwfn=multiwfn,
        )
        cmds = job.commands
        cmds.append("junk")
        assert "junk" not in job.commands

    def test_stdout_stderr_before_run(
        self, mock_wfn_file: Path, multiwfn: Multiwfn
    ) -> None:
        job = MultiwfnJob(
            input_file=mock_wfn_file,
            analysis=None,
            multiwfn=multiwfn,
        )
        assert job.stdout == ""
        assert job.stderr is None
        assert job.return_code is None
        assert job.execution_time is None
        assert job.success is None

    def test_timeout_property(
        self, mock_wfn_file: Path, multiwfn: Multiwfn
    ) -> None:
        job = MultiwfnJob(
            input_file=mock_wfn_file,
            analysis=None,
            multiwfn=multiwfn,
            timeout=120,
        )
        assert job.timeout == 120
        job.timeout = 300
        assert job.timeout == 300

    def test_verbose_property(
        self, mock_wfn_file: Path, multiwfn: Multiwfn
    ) -> None:
        job = MultiwfnJob(
            input_file=mock_wfn_file,
            analysis=None,
            multiwfn=multiwfn,
            verbose=True,
        )
        assert job.verbose is True


class TestJobAddCommands:
    """Tests for the add_commands() method of MultiwfnJob."""

    def test_add_custom_commands(
        self, mock_wfn_file: Path, multiwfn: Multiwfn
    ) -> None:
        job = MultiwfnJob(
            input_file=mock_wfn_file,
            analysis=None,
            multiwfn=multiwfn,
        )
        job.add_commands(["7", "1", "1"])
        assert "7" in job.commands

    def test_add_commands_extends(
        self, mock_wfn_file: Path, multiwfn: Multiwfn
    ) -> None:
        job = MultiwfnJob(
            input_file=mock_wfn_file,
            analysis=Menu.HIRSHFELD_CHARGE,
            multiwfn=multiwfn,
        )
        initial_len = len(job.commands)
        job.add_commands(["extra"])
        assert len(job.commands) == initial_len + 1


class TestJobFromFile:
    """Tests for the MultiwfnJob.from_file() factory method."""

    def test_from_file(self, mock_wfn_file: Path, multiwfn: Multiwfn) -> None:
        job = MultiwfnJob.from_file(
            input_file=mock_wfn_file,
            analysis=Menu.MAYER_BOND_ORDER,
            multiwfn=multiwfn,
        )
        assert isinstance(job, MultiwfnJob)
        assert "9" in job.commands

    def test_from_file_no_analysis(
        self, mock_wfn_file: Path, multiwfn: Multiwfn
    ) -> None:
        job = MultiwfnJob.from_file(
            input_file=mock_wfn_file,
            multiwfn=multiwfn,
        )
        assert job.commands == []


class TestJobStrRepr:
    """Tests for the __str__ and __repr__ methods of MultiwfnJob."""

    def test_str(self, mock_wfn_file: Path, multiwfn: Multiwfn) -> None:
        job = MultiwfnJob(
            input_file=mock_wfn_file,
            analysis=Menu.HIRSHFELD_CHARGE,
            multiwfn=multiwfn,
        )
        assert "mock.wfn" in str(job)

    def test_repr(self, mock_wfn_file: Path, multiwfn: Multiwfn) -> None:
        job = MultiwfnJob(
            input_file=mock_wfn_file,
            analysis=Menu.HIRSHFELD_CHARGE,
            multiwfn=multiwfn,
        )
        r = repr(job)
        assert "MultiwfnJob(" in r
        assert "input_file=" in r