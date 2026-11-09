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
* The one exception is ``NCI_ANALYSIS``: Multiwfn prints no NCI numbers
  to stdout, so its scatter file (``output.txt``) is read to add the
  main interaction regions ("spikes") to its result, and the
  sign(lambda2)rho cube's header is read to name the nearest atoms.
  The scatter file is assumed to belong to ``NCI_ANALYSIS``.
"""

from __future__ import annotations

import logging
import math
from collections import deque
from dataclasses import asdict
from pathlib import Path
from typing import Any

from pymultiwfn.analysis.result import MultiwfnResult, NCISpike
from pymultiwfn.enums.menu import Menu

logger = logging.getLogger(__name__)

RAW_STDOUT_SUFFIX = ".txt"

# NCI spike detection (NCI_ANALYSIS only): keep scatter points with
# RDG < NCI_RDG_CUTOFF and |sign(lambda2)rho| < NCI_DENSITY_CUTOFF (a.u.),
# group touching grid points into regions, drop regions smaller than
# NCI_MIN_POINTS. |sign(lambda2)rho| < NCI_VDW_THRESHOLD is van der Waals;
# beyond it, negative is attractive and positive repulsive.
NCI_RDG_CUTOFF = 0.5
NCI_DENSITY_CUTOFF = 0.05
NCI_VDW_THRESHOLD = 0.01
NCI_MIN_POINTS = 5
BOHR_PER_ANGSTROM = 1.0 / 0.529177210903
_ELEMENTS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe "
    "Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn "
    "Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W "
    "Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf "
    "Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()


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
            parsed = mwfn_result.to_dict()
            if menu is Menu.NCI_ANALYSIS:
                entry = parsed["results"][0]
                entry["spikes"] = _nci_spikes(entry, directory)
            analyses[menu_name] = {
                "parsed": parsed,
                "raw_stdout_file": rel,
            }

    return {"generated_files": generated_files, "analyses": analyses}


def _nci_spikes(entry: dict[str, Any], directory: Path) -> list[dict] | None:
    """Find the interaction regions in an NCI_ANALYSIS scatter file.

    Returns the spikes sorted by sign(lambda2)rho at their lowest-RDG
    point (most attractive first), or ``None`` if the scatter file or
    grid is missing or the file cannot be read.
    """
    scatter, origin, spacing = (
        entry.get("scatter_file"),
        entry.get("grid_origin_bohr"),
        entry.get("grid_spacing_bohr"),
    )
    if not scatter or origin is None or spacing is None:
        return None
    if not (directory / scatter).is_file():
        logger.warning("NCI_ANALYSIS: scatter file %s not found.", scatter)
        return None

    # Columns: X, Y, Z (Angstrom), sign(lambda2)rho, RDG (a.u.).
    points: list[tuple[float, float, float, float, float]] = []
    try:
        with open(directory / scatter, encoding="utf-8") as fh:
            for line in fh:
                fields = line.split()
                if len(fields) < 5:
                    continue
                value, rdg = float(fields[3]), float(fields[4])
                if rdg < NCI_RDG_CUTOFF and abs(value) < NCI_DENSITY_CUTOFF:
                    x, y, z = (float(v) for v in fields[:3])
                    points.append((x, y, z, value, rdg))
    except (OSError, ValueError) as exc:
        logger.warning("NCI_ANALYSIS: could not read %s: %s", scatter, exc)
        return None

    # Group touching grid points (26 neighbours). Neighbours join only if
    # they share a sign or both lie in the van der Waals band, so the
    # attractive and repulsive parts of one contact stay separate.
    cells = [
        tuple(round((p[i] * BOHR_PER_ANGSTROM - origin[i]) / spacing[i]) for i in range(3))
        for p in points
    ]
    index = {cell: n for n, cell in enumerate(cells)}
    offsets = [
        (i, j, k)
        for i in (-1, 0, 1)
        for j in (-1, 0, 1)
        for k in (-1, 0, 1)
        if (i, j, k) != (0, 0, 0)
    ]

    def joinable(a: int, b: int) -> bool:
        va, vb = points[a][3], points[b][3]
        same_sign = (va >= 0) == (vb >= 0)
        both_vdw = abs(va) < NCI_VDW_THRESHOLD and abs(vb) < NCI_VDW_THRESHOLD
        return same_sign or both_vdw

    seen = [False] * len(points)
    regions: list[list[int]] = []
    for start in range(len(points)):
        if seen[start]:
            continue
        seen[start] = True
        members, queue = [start], deque([start])
        while queue:
            cur = queue.popleft()
            ci, cj, ck = cells[cur]
            for di, dj, dk in offsets:
                nb = index.get((ci + di, cj + dj, ck + dk))
                if nb is not None and not seen[nb] and joinable(cur, nb):
                    seen[nb] = True
                    members.append(nb)
                    queue.append(nb)
        if len(members) >= NCI_MIN_POINTS:
            regions.append(members)

    cube = entry.get("sign_lambda2_rho_cube")
    atoms = _cube_atoms(directory / cube) if cube else []
    cell_volume = spacing[0] * spacing[1] * spacing[2]

    spikes = []
    for members in regions:
        region = [points[n] for n in members]
        low = min(region, key=lambda p: p[4])
        value = low[3]
        if value < -NCI_VDW_THRESHOLD:
            kind = "attractive"
        elif value > NCI_VDW_THRESHOLD:
            kind = "repulsive"
        else:
            kind = "van_der_waals"
        nearest = sorted(
            (math.dist([c * BOHR_PER_ANGSTROM for c in low[:3]], xyz) / BOHR_PER_ANGSTROM, name)
            for name, xyz in atoms
        )[:2]
        spikes.append(
            NCISpike(
                interaction_type=kind,
                sign_lambda2_rho_at_min=value,
                rdg_min=low[4],
                sign_lambda2_rho_range=[
                    min(p[3] for p in region),
                    max(p[3] for p in region),
                ],
                n_points=len(region),
                volume_bohr3=round(len(region) * cell_volume, 6),
                min_position_angstrom=list(low[:3]),
                centroid_angstrom=[
                    round(sum(p[i] for p in region) / len(region), 6) for i in range(3)
                ],
                nearest_atoms=[
                    {"atom": name, "distance_angstrom": round(d, 4)}
                    for d, name in nearest
                ]
                or None,
            )
        )
    spikes.sort(key=lambda s: s.sign_lambda2_rho_at_min)
    return [asdict(s) for s in spikes]


def _cube_atoms(cube_path: Path) -> list[tuple[str, list[float]]]:
    """Read atom labels (e.g. ``"Se1"``) and Bohr positions from a cube header."""
    try:
        with open(cube_path, encoding="utf-8") as fh:
            lines = [fh.readline() for _ in range(6)]
            n_atoms = abs(int(lines[2].split()[0]))
            atoms = []
            for i in range(n_atoms):
                fields = fh.readline().split()
                z = int(fields[0])
                symbol = _ELEMENTS[z] if 0 < z < len(_ELEMENTS) else f"Z{z}"
                atoms.append((f"{symbol}{i + 1}", [float(v) for v in fields[2:5]]))
    except (OSError, ValueError, IndexError) as exc:
        logger.warning("Could not read atoms from %s: %s", cube_path, exc)
        return []
    return atoms
