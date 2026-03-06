"""Tests for pymultiwfn.api.outcome — MultiwfnJobOutcome."""

from pathlib import Path

import pytest

from pymultiwfn.api.outcome import MultiwfnJobOutcome


@pytest.fixture
def success_outcome() -> MultiwfnJobOutcome:
    return MultiwfnJobOutcome(
        stdout="Hirshfeld charges computed",
        stderr="",
        return_code=0,
        execution_time=1.5,
    )


@pytest.fixture
def failed_outcome() -> MultiwfnJobOutcome:
    return MultiwfnJobOutcome(
        stdout="partial output",
        stderr="error: cannot open file",
        return_code=1,
        execution_time=0.3,
    )


@pytest.fixture
def signal_killed_outcome() -> MultiwfnJobOutcome:
    return MultiwfnJobOutcome(
        stdout="", stderr="", return_code=-9, execution_time=0.0
    )


class TestIsSuccessful:
    def test_success(self, success_outcome: MultiwfnJobOutcome) -> None:
        assert success_outcome.is_successful() is True

    def test_failure_stderr_error(
        self, failed_outcome: MultiwfnJobOutcome
    ) -> None:
        assert failed_outcome.is_successful() is False

    def test_failure_negative_rc(
        self, signal_killed_outcome: MultiwfnJobOutcome
    ) -> None:
        assert signal_killed_outcome.is_successful() is False

    def test_failure_none_rc(self) -> None:
        outcome = MultiwfnJobOutcome(
            stdout="", stderr="", return_code=None, execution_time=0.0
        )
        assert outcome.is_successful() is False

    def test_nonzero_rc_but_clean_stderr(self) -> None:
        """Multiwfn often returns non-zero on success in batch mode."""
        outcome = MultiwfnJobOutcome(
            stdout="result", stderr="", return_code=1, execution_time=1.0
        )
        assert outcome.is_successful() is True

    def test_fatal_in_stderr(self) -> None:
        outcome = MultiwfnJobOutcome(
            stdout="", stderr="fatal: segfault", return_code=0,
            execution_time=0.0,
        )
        assert outcome.is_successful() is False

    @pytest.mark.parametrize(
        "indicator",
        ["error", "fatal", "cannot open", "not found", "failed to"],
    )
    def test_all_error_indicators(self, indicator: str) -> None:
        outcome = MultiwfnJobOutcome(
            stdout="", stderr=f"something {indicator} something",
            return_code=0, execution_time=0.0,
        )
        assert outcome.is_successful() is False


class TestSaveOutput:
    def test_save_and_read(
        self, success_outcome: MultiwfnJobOutcome, temp_dir: Path
    ) -> None:
        path = temp_dir / "output.log"
        success_outcome.save_output(path)
        assert path.exists()
        assert path.read_text() == success_outcome.stdout