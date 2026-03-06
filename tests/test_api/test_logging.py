"""Tests for pymultiwfn.api.logging — BatchLogger."""

from pathlib import Path

import pytest

from pymultiwfn.api.logging import (
    BatchLogger,
    JobLogEntry,
    JobStatus,
)
from pymultiwfn.api.outcome import MultiwfnJobOutcome
from pymultiwfn.enums.menu import Menu


@pytest.fixture
def logger(temp_dir: Path) -> BatchLogger:
    return BatchLogger(log_dir=temp_dir)


@pytest.fixture
def success_outcome() -> MultiwfnJobOutcome:
    return MultiwfnJobOutcome(
        stdout="Charges computed\nAtom 1(C): 0.05",
        stderr="",
        return_code=0,
        execution_time=2.5,
    )


@pytest.fixture
def failed_outcome() -> MultiwfnJobOutcome:
    return MultiwfnJobOutcome(
        stdout="partial",
        stderr="error: cannot open file.wfn",
        return_code=1,
        execution_time=0.1,
    )


class TestBatchLoggerLifecycle:
    """Tests for the lifecycle of BatchLogger."""

    def test_log_file_created(self, logger: BatchLogger) -> None:
        assert logger.log_path.exists()

    def test_log_filename_auto(self, logger: BatchLogger) -> None:
        assert logger.log_path.name.startswith("batch_")
        assert logger.log_path.suffix == ".log"

    def test_custom_filename(self, temp_dir: Path) -> None:
        lg = BatchLogger(log_dir=temp_dir, log_filename="custom.log")
        assert lg.log_path.name == "custom.log"

    def test_start_and_end_batch(self, logger: BatchLogger) -> None:
        logger.start_batch(
            files=["mol1.wfn"], analyses=[Menu.HIRSHFELD_CHARGE]
        )
        logger.end_batch()
        content = logger.log_path.read_text()
        assert "PYMULTIWFN BATCH LOG" in content
        assert "BATCH SUMMARY" in content
        assert "mol1.wfn" in content
        assert "HIRSHFELD_CHARGE" in content

    def test_context_manager(self, temp_dir: Path) -> None:
        with BatchLogger(log_dir=temp_dir) as lg:
            lg.start_batch(files=["a.wfn"])
        content = lg.log_path.read_text()
        assert "BATCH SUMMARY" in content


class TestJobLogging:
    """Tests for logging individual jobs."""

    def test_log_successful_job(
        self, logger: BatchLogger, success_outcome: MultiwfnJobOutcome
    ) -> None:
        logger.start_batch(files=["mol.wfn"])
        entry = logger.log_job_start(
            input_file="mol.wfn",
            analysis_name="HIRSHFELD_CHARGE",
            commands=["7", "1", "1", "n", "0"],
        )
        logger.log_job_end(entry=entry, outcome=success_outcome)
        logger.end_batch()

        assert entry.status == JobStatus.SUCCESS
        assert entry.execution_time == pytest.approx(2.5)
        assert entry.return_code == 0

        content = logger.log_path.read_text()
        assert "SUCCESS" in content
        assert "7 -> 1 -> 1 -> n -> 0" in content

    def test_log_failed_job(
        self, logger: BatchLogger, failed_outcome: MultiwfnJobOutcome
    ) -> None:
        logger.start_batch(files=["mol.wfn"])
        entry = logger.log_job_start(
            input_file="mol.wfn",
            analysis_name="HIRSHFELD_CHARGE",
            commands=["7", "1"],
        )
        logger.log_job_end(entry=entry, outcome=failed_outcome)
        logger.end_batch()

        assert entry.status == JobStatus.FAILED
        assert entry.error_message is not None
        assert "cannot open" in entry.error_message

        content = logger.log_path.read_text()
        assert "FAILED" in content
        assert "FAILED / TIMED-OUT JOBS" in content

    def test_log_job_with_exception(self, logger: BatchLogger) -> None:
        from pymultiwfn.api.exceptions import MultiwfnError

        logger.start_batch(files=["mol.wfn"])
        entry = logger.log_job_start(
            input_file="mol.wfn",
            analysis_name="CUBE_DENSITY",
            commands=["5", "1"],
        )
        error = MultiwfnError("timed out after 60s")
        logger.log_job_end(entry=entry, error=error)
        logger.end_batch()

        assert entry.status == JobStatus.TIMEOUT
        assert entry.error_message is not None
        assert "MultiwfnError" in entry.error_message

    def test_log_skipped_job(self, logger: BatchLogger) -> None:
        logger.start_batch(files=["mol.wfn"])
        logger.log_job_skipped(
            input_file="mol.wfn",
            analysis_name="HIRSHFELD_CHARGE",
            reason="cached",
        )
        logger.end_batch()

        assert len(logger.entries) == 1
        assert logger.entries[0].status == JobStatus.SKIPPED

        content = logger.log_path.read_text()
        assert "SKIPPED" in content

    def test_no_outcome_no_error(self, logger: BatchLogger) -> None:
        logger.start_batch(files=["mol.wfn"])
        entry = logger.log_job_start(
            input_file="mol.wfn",
            analysis_name="TEST",
            commands=[],
        )
        logger.log_job_end(entry=entry, outcome=None, error=None)
        logger.end_batch()
        assert entry.status == JobStatus.FAILED
        assert entry.error_message is not None
        assert "no outcome" in entry.error_message.lower()


class TestBatchSummary:
    """Tests that the batch summary correctly counts successes."""

    def test_summary_counts(
        self,
        logger: BatchLogger,
        success_outcome: MultiwfnJobOutcome,
        failed_outcome: MultiwfnJobOutcome,
    ) -> None:
        logger.start_batch(files=["a.wfn", "b.wfn"])

        e1 = logger.log_job_start("a.wfn", "HIRSHFELD_CHARGE", ["7"])
        logger.log_job_end(e1, outcome=success_outcome)

        e2 = logger.log_job_start("b.wfn", "HIRSHFELD_CHARGE", ["7"])
        logger.log_job_end(e2, outcome=failed_outcome)

        logger.log_job_skipped("a.wfn", "MAYER_BOND_ORDER")

        logger.end_batch()

        content = logger.log_path.read_text()
        assert "Succeeded         : 1" in content
        assert "Failed            : 1" in content
        assert "Skipped (cached)  : 1" in content
        assert "Total jobs        : 3" in content


class TestJobStatusEnum:
    """Tests for the JobStatus enum."""

    def test_all_values(self) -> None:
        assert len(JobStatus) == 6
        assert JobStatus.SUCCESS.value == "SUCCESS"
        assert JobStatus.TIMEOUT.value == "TIMEOUT"


class TestJobLogEntry:
    """Tests for the JobLogEntry dataclass."""

    def test_defaults(self) -> None:
        entry = JobLogEntry(input_file="a.wfn", analysis_name="TEST")
        assert entry.status == JobStatus.PENDING
        assert entry.commands == []
        assert entry.start_time is None
        assert entry.error_message is None
