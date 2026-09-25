"""Discover and parse the output of a Multiwfn analysis directory.

Every analysis for a molecule writes into one flat output directory
(``<input_file_name>.output/``): its captured stdout as
``<ANALYSIS_NAME>.txt`` (named after the ``Menu`` sequence that
produced it), plus whatever other files Multiwfn generated for it
(cubes, images, exported structures, etc.). Rather than parsing stdout
while a job is running, :func:`scan_output_directory` does a
*subsequent* pass over that directory:

* Each stdout file whose stem matches a ``Menu`` member name is mapped
  back to that member and re-parsed from disk using the existing
  regex-based parsers.
* Every other file is only ever recorded by path. Cube files, images,
  and any other non-stdout output are not opened or interpreted —
  only their location is stored.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pymultiwfn.analysis.result import MultiwfnResult
from pymultiwfn.enums.menu import Menu

RAW_STDOUT_SUFFIX = ".txt"


def scan_output_directory(
    directory: Path,
    exclude: set[Path] | None = None,
) -> dict[str, Any]:
    """Discover every file Multiwfn wrote into *directory*.

    Parameters
    ----------
    directory
        The molecule's flat output directory to scan.
    exclude
        Paths to skip (e.g. the per-molecule JSON file, which lives
        inside the same directory it describes).

    Returns
    -------
    A JSON-safe dict with two keys:

    ``generated_files``
        Paths (relative to *directory*) of every file found.
    ``analyses``
        Mapping of ``Menu`` name to that analysis's stdout-derived
        parsed result, re-parsed from its saved ``<ANALYSIS_NAME>.txt``
        file. All other files are only ever listed in
        ``generated_files``, never parsed.
    """
    exclude = exclude or set()
    generated_files: list[str] = []
    analyses: dict[str, Any] = {}

    if not directory.exists():
        return {"generated_files": generated_files, "analyses": analyses}

    for path in sorted(directory.rglob("*")):
        if not path.is_file() or path in exclude:
            continue
        rel = str(path.relative_to(directory))
        generated_files.append(rel)

        if not path.name.endswith(RAW_STDOUT_SUFFIX):
            continue

        menu_name = path.name[: -len(RAW_STDOUT_SUFFIX)]
        try:
            menu = Menu[menu_name]
        except KeyError:
            continue
        stdout = path.read_text(encoding="utf-8", errors="replace")
        mwfn_result = MultiwfnResult(analysis=menu)
        mwfn_result.parse(stdout)
        if mwfn_result.result:
            analyses[menu_name] = {
                "parsed": mwfn_result.to_dict(),
                "raw_stdout_file": rel,
            }

    return {"generated_files": generated_files, "analyses": analyses}
