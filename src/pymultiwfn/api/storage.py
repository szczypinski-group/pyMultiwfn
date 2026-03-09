"""Per-molecule JSON result store for pymultiwfn.

Each wavefunction file gets a companion ``.json`` file (e.g.
``benzene.wfn.json``) that accumulates the **parsed** output from every
analysis run on that molecule.

Only parsed representations (from ``parsers.py``) are stored — never raw
stdout and never large binary artefacts (cube files, grids, etc.).
Spectrum peak data is stored in the JSON; spectrum *figures* are written
as separate image files.

This module is an internal implementation detail. It is called
automatically by :class:`MultiwfnAnalysis` during :meth:`run`.

"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any

from pymultiwfn.analysis.routing import ParserRoute
from pymultiwfn.enums.menu import Menu

# ─────────────────────────────────────────────────────────────────────────────
# Parser routing table
# ─────────────────────────────────────────────────────────────────────────────
# Maps each Menu member to a callable ``(stdout: str) -> Any`` that returns
# the parsed result to be stored in the JSON.  Menu items that produce only
# non-standard output (cubes, GUI-only, interactive) are mapped to ``None``
# and will be skipped.
#
# When a Menu member is absent from this dict it is also skipped.
# ─────────────────────────────────────────────────────────────────────────────


class ResultStore:
    """Per-molecule JSON result store.

    Manages a ``<molecule>.json`` file that accumulates parsed results
    keyed by analysis name. The file is read on construction (if it
    exists) and written back after every new result is added.

    This class is an internal implementation detail.

    """

    def __init__(self, input_file: Path, work_dir: Path) -> None:
        self._input_file = Path(input_file)
        self._work_dir = Path(work_dir)
        self._work_dir.mkdir(parents=True, exist_ok=True)

        # e.g. benzene.wfn -> benzene.wfn.json
        self._json_path = self._work_dir / f"{self._input_file.name}.json"
        self._data: dict[str, Any] = self._load()

    @property
    def json_path(self) -> Path:
        return self._json_path

    @property
    def data(self) -> dict[str, Any]:
        """Return a copy of the stored data."""
        return dict(self._data)

    # ── persistence ──────────────────────────────────────────────────────

    def _load(self) -> dict[str, Any]:
        """Load existing JSON or return an empty structure."""
        if self._json_path.exists():
            with Path.open(self._json_path, encoding="utf-8") as f:
                return json.load(f)
        return {
            "input_file": str(self._input_file),
            "analyses": {},
        }

    def _save(self) -> None:
        """Write the current data to disk."""
        with Path.open(self._json_path, "w", encoding="utf-8") as f:
            json.dump(self._data, f, indent=2, default=str)

    # ── read / write analyses ────────────────────────────────────────────

    def has_result(self, analysis: Menu) -> bool:
        """Check whether a parsed result already exists for *analysis*."""
        return analysis.name in self._data.get("analyses", {})

    def get_result(self, analysis: Menu) -> dict[str, Any] | None:
        """Retrieve a previously stored parsed result, or ``None``."""
        return self._data.get("analyses", {}).get(analysis.name)

    def store_result(
        self,
        analysis: Menu,
        stdout: str,
    ) -> dict[str, Any] | None:
        """Parse *stdout* for *analysis* and persist the result.

        Returns the parsed dict, or ``None`` if no parser is available
        for this analysis type.
        """
        parser_fn = ParserRoute.ROUTE_TABLE.get(analysis)
        if parser_fn is None:
            return None

        parsed = parser_fn(stdout)
        if not parsed:
            return None

        entry: dict[str, Any] = {
            "parsed": parsed,
            "timestamp": datetime.now().isoformat(),
        }

        self._data.setdefault("analyses", {})[analysis.name] = entry
        self._save()

        return parsed
