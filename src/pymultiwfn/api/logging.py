"""Batch calculation logging for pymultiwfn.

This module is an internal implementation detail. Users do not create or
interact with :class:`BatchLogger` directly — it is instantiated
automatically when :meth:`MultiwfnAnalysis.run` executes and a log file
is written alongside the calculation outputs.

Uses Python's standard :mod:`logging` package with a dedicated
``FileHandler`` per batch run.

"""

import logging
import textwrap
import time
from collections.abc import Sequence
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any
from uuid import uuid4

from pymultiwfn.api.outcome import MultiwfnJobOutcome

# ─────────────────────────────────────────────────────────────────────────────
# Data containers
# ─────────────────────────────────────────────────────────────────────────────


class JobStatus(Enum):
    """Status of a single calculation within a batch."""

    PENDING = "PENDING"
    RUNNING = "RUNNING"
    SUCCESS = "SUCCESS"
    FAILED = "FAILED"
    TIMEOUT = "TIMEOUT"
    SKIPPED = "SKIPPED"


@dataclass
class JobLogEntry:
    """Record for a single Multiwfn job execution within a batch."""

    input_file: str
    analysis_name: str
    commands: list[str] = field(default_factory=list)
    status: JobStatus = JobStatus.PENDING
    start_time: datetime | None = None
    end_time: datetime | None = None
    execution_time: float | None = None
    return_code: int | None = None
    error_message: str | None = None
    stderr_snippet: str | None = None
    stdout_snippet: str | None = None


@dataclass
class BatchLogSummary:
    """Aggregate statistics for the entire batch."""

    total_files: int = 0
    total_jobs: int = 0
    succeeded: int = 0
    failed: int = 0
    timed_out: int = 0
    skipped: int = 0
    total_wall_time: float = 0.0
    total_execution_time: float = 0.0


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

_SNIPPET_LENGTH = 2000
_SEP_WIDTH = 78
_LOG_FMT = "%(asctime)s | %(levelname)-8s | %(message)s"
_LOG_DATEFMT = "%Y-%m-%d %H:%M:%S"


# ─────────────────────────────────────────────────────────────────────────────
# BatchLogger
# ─────────────────────────────────────────────────────────────────────────────


class BatchLogger:
    """Internal logger that writes a plain-text log file for a batch run.

    This class is **not** part of the public API. It is created and managed
    automatically by :class:`MultiwfnAnalysis` during execution.

    Wraps a :class:`logging.Logger` with a dedicated
    :class:`logging.FileHandler` that is removed and closed when the
    batch ends.

    """

    def __init__(
        self,
        log_dir: Path,
        log_filename: str | None = None,
        snippet_length: int = _SNIPPET_LENGTH,
    ) -> None:
        self._log_dir = Path(log_dir)
        self._log_dir.mkdir(parents=True, exist_ok=True)

        if log_filename is None:
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            short_id = uuid4().hex[:8]
            log_filename = f"batch_{ts}_{short_id}.log"

        self._log_path = self._log_dir / log_filename
        self._snippet_length = snippet_length

        self._entries: list[JobLogEntry] = []
        self._batch_start: datetime | None = None
        self._batch_end: datetime | None = None
        self._batch_wall_start: float | None = None

        self._input_files: list[str] = []
        self._analysis_names: list[str] = []

        # Namespaced logger — unique per batch so concurrent runs don't
        # collide.
        logger_name = f"pymultiwfn.batch.{uuid4().hex[:8]}"
        self._logger = logging.getLogger(logger_name)
        self._logger.setLevel(logging.DEBUG)
        self._logger.propagate = False

        formatter = logging.Formatter(fmt=_LOG_FMT, datefmt=_LOG_DATEFMT)
        self._file_handler = logging.FileHandler(
            self._log_path,
            mode="w",
            encoding="utf-8",
        )
        self._file_handler.setLevel(logging.DEBUG)
        self._file_handler.setFormatter(formatter)
        self._logger.addHandler(self._file_handler)

    # ── public properties ────────────────────────────────────────────────

    @property
    def log_path(self) -> Path:
        return self._log_path

    @property
    def entries(self) -> list[JobLogEntry]:
        return list(self._entries)

    # ── batch lifecycle ──────────────────────────────────────────────────

    def start_batch(
        self,
        files: Sequence[str | Path],
        analyses: Sequence[Any] | None = None,
    ) -> None:
        self._batch_start = datetime.now()
        self._batch_wall_start = time.monotonic()
        self._input_files = [str(f) for f in files]

        if analyses is not None:
            self._analysis_names = [
                a.name if hasattr(a, "name") else str(a) for a in analyses
            ]

        self._write_header()

    def end_batch(self) -> None:
        self._batch_end = datetime.now()
        wall_time = (
            time.monotonic() - self._batch_wall_start
            if self._batch_wall_start is not None
            else 0.0
        )

        summary = self._compute_summary(wall_time)
        self._write_summary(summary)
        self._write_failed_job_details()
        self._write_footer()
        self._close()

    # ── per-job logging ──────────────────────────────────────────────────

    def log_job_start(
        self,
        input_file: str | Path,
        analysis_name: str,
        commands: list[str],
    ) -> JobLogEntry:
        entry = JobLogEntry(
            input_file=str(input_file),
            analysis_name=analysis_name,
            commands=list(commands),
            status=JobStatus.RUNNING,
            start_time=datetime.now(),
        )
        self._entries.append(entry)

        seq_str = " -> ".join(entry.commands) if entry.commands else "(none)"
        self._logger.info(
            "JOB START | %s | file: %s",
            entry.analysis_name,
            entry.input_file,
        )
        self._logger.info("  Sequence: %s", seq_str)

        return entry

    def log_job_end(
        self,
        entry: JobLogEntry,
        outcome: MultiwfnJobOutcome | None = None,
        error: Exception | None = None,
    ) -> None:
        entry.end_time = datetime.now()

        if error is not None:
            entry.status = self._classify_error(error)
            entry.error_message = self._describe_error(error)
            if outcome is not None:
                self._populate_from_outcome(entry, outcome)
        elif outcome is not None:
            self._populate_from_outcome(entry, outcome)
            if outcome.is_successful():
                entry.status = JobStatus.SUCCESS
            else:
                entry.status = JobStatus.FAILED
                entry.error_message = self._describe_outcome_failure(outcome)
        else:
            entry.status = JobStatus.FAILED
            entry.error_message = (
                "Job produced no outcome and no exception was raised."
            )

        self._write_job_end(entry)

    def log_job_skipped(
        self,
        input_file: str | Path,
        analysis_name: str,
        reason: str = "cached result exists",
    ) -> None:
        entry = JobLogEntry(
            input_file=str(input_file),
            analysis_name=analysis_name,
            status=JobStatus.SKIPPED,
            start_time=datetime.now(),
            end_time=datetime.now(),
            error_message=reason,
        )
        self._entries.append(entry)
        self._logger.info(
            "SKIPPED | %s | file: %s | reason: %s",
            analysis_name,
            str(input_file),
            reason,
        )

    # ── error helpers ────────────────────────────────────────────────────

    @staticmethod
    def _classify_error(error: Exception) -> JobStatus:
        msg = str(error).lower()
        if "timeout" in msg or "timed out" in msg:
            return JobStatus.TIMEOUT
        return JobStatus.FAILED

    @staticmethod
    def _describe_error(error: Exception) -> str:
        cls_name = type(error).__qualname__
        return f"[{cls_name}] {error}"

    @staticmethod
    def _describe_outcome_failure(outcome: MultiwfnJobOutcome) -> str:
        parts: list[str] = []
        if outcome.return_code is not None and outcome.return_code < 0:
            parts.append(
                f"Process terminated with signal {-outcome.return_code}."
            )
        error_indicators = [
            "error",
            "fatal",
            "cannot open",
            "not found",
            "failed to",
        ]
        found = [
            ind for ind in error_indicators if ind in outcome.stderr.lower()
        ]
        if found:
            parts.append(
                f"stderr contains error indicator(s): {', '.join(found)}."
            )
            for line in outcome.stderr.splitlines():
                if any(ind in line.lower() for ind in found):
                    parts.append(f"  >> {line.strip()}")
                    break
        if not parts:
            parts.append(
                "Unknown failure (no specific error indicator found)."
            )
        return " ".join(parts)

    def _populate_from_outcome(
        self,
        entry: JobLogEntry,
        outcome: MultiwfnJobOutcome,
    ) -> None:
        entry.return_code = outcome.return_code
        entry.execution_time = outcome.execution_time
        if outcome.stderr:
            entry.stderr_snippet = outcome.stderr[: self._snippet_length]
        if outcome.stdout:
            entry.stdout_snippet = outcome.stdout[-self._snippet_length :]

    # ── summary ──────────────────────────────────────────────────────────

    def _compute_summary(self, wall_time: float) -> BatchLogSummary:
        summary = BatchLogSummary(
            total_files=len(set(e.input_file for e in self._entries)),
            total_jobs=len(self._entries),
            total_wall_time=wall_time,
        )
        for entry in self._entries:
            if entry.status == JobStatus.SUCCESS:
                summary.succeeded += 1
            elif entry.status == JobStatus.FAILED:
                summary.failed += 1
            elif entry.status == JobStatus.TIMEOUT:
                summary.timed_out += 1
            elif entry.status == JobStatus.SKIPPED:
                summary.skipped += 1
            if entry.execution_time is not None:
                summary.total_execution_time += entry.execution_time
        return summary

    # ── writing via logging ──────────────────────────────────────────────

    def _write_header(self) -> None:
        self._logger.info("=" * _SEP_WIDTH)
        self._logger.info("  PYMULTIWFN BATCH LOG")
        self._logger.info("=" * _SEP_WIDTH)
        self._logger.info("  Log file: %s", self._log_path)
        self._logger.info("")
        self._logger.info("  Input files (%d):", len(self._input_files))
        for f in self._input_files:
            self._logger.info("    - %s", f)
        if self._analysis_names:
            self._logger.info("")
            self._logger.info(
                "  Analyses (%d):",
                len(self._analysis_names),
            )
            for a in self._analysis_names:
                self._logger.info("    - %s", a)
        self._logger.info("=" * _SEP_WIDTH)
        self._logger.info("")

    def _write_job_end(self, entry: JobLogEntry) -> None:
        elapsed = (
            f"{entry.execution_time:.2f}s"
            if entry.execution_time is not None
            else "n/a"
        )

        if entry.status == JobStatus.SUCCESS:
            self._logger.info(
                "JOB END   | %s | file: %s | %s | time: %s | rc: %s",
                entry.analysis_name,
                entry.input_file,
                entry.status.value,
                elapsed,
                entry.return_code,
            )
        else:
            self._logger.error(
                "JOB END   | %s | file: %s | %s | time: %s | rc: %s",
                entry.analysis_name,
                entry.input_file,
                entry.status.value,
                elapsed,
                entry.return_code,
            )
            self._write_error_block(entry)

    def _write_error_block(self, entry: JobLogEntry) -> None:
        if entry.error_message:
            for line in textwrap.wrap(
                entry.error_message,
                width=_SEP_WIDTH - 4,
            ):
                self._logger.error("  %s", line)
        if entry.stderr_snippet:
            self._logger.debug(
                "  STDERR (first %d chars):",
                self._snippet_length,
            )
            for line in entry.stderr_snippet.splitlines()[:30]:
                self._logger.debug("    | %s", line)
        if entry.stdout_snippet:
            self._logger.debug(
                "  STDOUT (last %d chars):",
                self._snippet_length,
            )
            for line in entry.stdout_snippet.splitlines()[-30:]:
                self._logger.debug("    | %s", line)

    def _write_summary(self, summary: BatchLogSummary) -> None:
        self._logger.info("=" * _SEP_WIDTH)
        self._logger.info("  BATCH SUMMARY")
        self._logger.info("=" * _SEP_WIDTH)
        self._logger.info("  Total files       : %d", summary.total_files)
        self._logger.info("  Total jobs        : %d", summary.total_jobs)
        self._logger.info("  Succeeded         : %d", summary.succeeded)
        self._logger.info("  Failed            : %d", summary.failed)
        self._logger.info("  Timed out         : %d", summary.timed_out)
        self._logger.info("  Skipped (cached)  : %d", summary.skipped)
        self._logger.info(
            "  Wall time         : %.2fs",
            summary.total_wall_time,
        )
        self._logger.info(
            "  Multiwfn CPU time : %.2fs",
            summary.total_execution_time,
        )
        self._logger.info("")

        files_seen: dict[str, dict[str, int]] = {}
        for entry in self._entries:
            f = entry.input_file
            files_seen.setdefault(f, {s.value: 0 for s in JobStatus})
            files_seen[f][entry.status.value] += 1

        self._logger.info("  Per-file breakdown:")
        for fname, counts in files_seen.items():
            parts = [f"{k}={v}" for k, v in counts.items() if v > 0]
            self._logger.info("    %s: %s", fname, ", ".join(parts))
        self._logger.info("")

    def _write_failed_job_details(self) -> None:
        failures = [
            e
            for e in self._entries
            if e.status in (JobStatus.FAILED, JobStatus.TIMEOUT)
        ]
        if not failures:
            return

        self._logger.warning("-" * _SEP_WIDTH)
        self._logger.warning(
            "  FAILED / TIMED-OUT JOBS (%d)",
            len(failures),
        )
        self._logger.warning("-" * _SEP_WIDTH)
        for i, entry in enumerate(failures, 1):
            self._logger.warning(
                "  %d. [%s] %s | %s",
                i,
                entry.status.value,
                entry.analysis_name,
                entry.input_file,
            )
            if entry.error_message:
                for line in textwrap.wrap(
                    entry.error_message,
                    width=_SEP_WIDTH - 8,
                ):
                    self._logger.warning("       %s", line)
        self._logger.warning("")

    def _write_footer(self) -> None:
        self._logger.info("=" * _SEP_WIDTH)
        end_str = (
            self._batch_end.strftime("%Y-%m-%d %H:%M:%S")
            if self._batch_end is not None
            else "unknown"
        )
        self._logger.info("  Finished: %s", end_str)
        self._logger.info("=" * _SEP_WIDTH)

    # ── teardown ─────────────────────────────────────────────────────────

    def _close(self) -> None:
        self._file_handler.flush()
        self._file_handler.close()
        self._logger.removeHandler(self._file_handler)

    def __enter__(self) -> "BatchLogger":
        return self

    def __exit__(
        self,
        exc_type: type | None,
        exc_val: BaseException | None,
        exc_tb: Any,
    ) -> None:
        if self._batch_start is not None and self._batch_end is None:
            self.end_batch()
        else:
            self._close()
