import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import marimo as mo

    mo.md(
        r"""
    # High-throughput optoelectronic screening with pyMultiwfn

    Optoelectronic materials — and dyes for dye-sensitised solar cells in
    particular — are a major target for computational materials design,
    because a good dye needs several properties at once: a small optical
    gap, a favourable HOMO&ndash;LUMO gap, a small exciton binding energy, and
    a matched ionization potential / electron affinity. Screening thousands
    of candidate dyes for this combination by hand is exactly the kind of
    repetitive, multi-step analysis that high-throughput virtual screening
    workflows are built for — see Heath-Apostolopoulos *et al.*, who screened
    ~45,000 candidate dyes built from a set of scaffolds and moieties.

    Computing each of those properties in plain Multiwfn means driving a
    *different* interactive menu by hand for every molecule. `pyMultiwfn`
    turns that into one call: give it a file and the list of `Menu` entries
    you want, and it drives Multiwfn through all of them, parses every
    stdout into structured data, and writes one JSON summary per molecule.

    This notebook shows the first step of that pipeline on a real dye from
    an ORCA TD-DFT calculation: reading its structure and its
    HOMO&ndash;LUMO gap from the wavefunction file
    (`optoelectronic.molden`), the starting point before any of the other
    properties above.
    """
    )
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Before you start

    A few things that matter for every analysis, not just this one:

    - **Two files from one ORCA job.** This folder holds both outputs of the
      same TD-DFT/TDA calculation, and each carries different information:
        - `tddft.out` is ORCA's text log. It has the geometry and the
          excited states (excitation energies and oscillator strengths), but
          Multiwfn reads *only the geometry* from it — not the molecular
          orbitals. Run `Menu.VIEW_STRUCTURE` on it and you get the atoms but
          no HOMO&ndash;LUMO gap.
        - `optoelectronic.molden` is the wavefunction, exported from the
          job's `.gbw` file with `orca_2mkl optoelectric -molden`. It has the
          molecular orbitals and their energies, which is what
          `Menu.VIEW_STRUCTURE` needs to report the HOMO&ndash;LUMO gap, and
          what any density-based analysis (surface ESP, hole&ndash;electron
          analysis) needs too.
    - **File type matters.** Multiwfn doesn't accept every input for every
      calculation, so pick the file by what the analysis needs: orbitals and
      densities from a `.molden`/`.wfn`/`.wfx`/`.fchk` file, excited-state
      data (UV&ndash;Vis spectra) from a Gaussian or ORCA `.out`/`.log` file.
      `pyMultiwfn` checks the extension up front and raises
      `InvalidInputFileError` immediately if it won't work, rather than
      letting Multiwfn fail halfway through a calculation.
    - **Every analysis writes to one place.** All output for a given input
      file — raw Multiwfn stdout, any cube/graphics files, and a combined
      JSON summary — lands in `<input_file_name>.output/` next to wherever
      you tell it to run. You never have to hunt across several working
      directories for one molecule's results.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Reading the structure and HOMO&ndash;LUMO gap

    Point this at `optoelectronic.molden` (or any other wavefunction file)
    and run the cell below.
    """)
    return


@app.cell
def _(mo):
    import json
    from collections import Counter

    from pymultiwfn import Menu, MultiwfnAnalysis

    # Paths are relative to this notebook's folder, wherever marimo is
    # launched from.
    notebook_dir = mo.notebook_dir()

    dye_analysis = MultiwfnAnalysis(
        notebook_dir / "optoelectronic.molden",
        analyses=[Menu.VIEW_STRUCTURE],
    )
    dye_analysis.run(work_dir=notebook_dir / "tutorial_results")
    return Counter, dye_analysis, json, notebook_dir


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Reading the results

    Every run writes one JSON summary per molecule, inside
    `tutorial_results/<input_file_name>.output/`. It is rebuilt from that
    folder on every run, so it holds the same data whether the analysis was
    just computed or loaded from the cache — read results from there.
    `Menu.VIEW_STRUCTURE` stores one entry per atom plus one entry with the
    frontier orbitals.
    """)
    return


@app.cell
def _(Counter, dye_analysis, json, mo, notebook_dir):
    dye_summary_path = (
        notebook_dir
        / "tutorial_results"
        / f"{dye_analysis.input_file.name}.output"
        / f"{dye_analysis.input_file.stem}.json"
    )
    structure = json.loads(dye_summary_path.read_text())["analyses"][
        "VIEW_STRUCTURE"
    ]["parsed"]["results"]

    # One entry per atom, keyed by atom number.
    dye_atoms = {
        entry["atom_id"]: entry for entry in structure if "atom_id" in entry
    }
    composition = Counter(atom["element"] for atom in dye_atoms.values())
    formula = " ".join(
        f"{element}{count}" for element, count in sorted(composition.items())
    )
    frontier = next(entry for entry in structure if "gap_eV" in entry)

    mo.md(
        f"""
    - **Structure**: `Menu.VIEW_STRUCTURE` read {len(dye_atoms)} atoms
      ({formula}) from the wavefunction file. Check this against the
      structure you meant to run before trusting anything computed on it.
    - **HOMO** (orbital {frontier["homo_index"]}):
      {frontier["homo_energy_eV"]:.3f} eV
    - **LUMO** (orbital {frontier["lumo_index"]}):
      {frontier["lumo_energy_eV"]:.3f} eV
    - **HOMO&ndash;LUMO gap**: {frontier["gap_eV"]:.3f} eV
      ({frontier["gap_kJ_mol"]:.1f} kJ/mol)

    For a dye this gap is a first filter: it tracks how easily the molecule
    can be excited and how its frontier levels sit against the electrodes it
    has to inject charge into. It is a ground-state orbital-energy
    difference, though, not the optical gap: the energy of the first bright
    excitation comes from the TD-DFT excited states in `tddft.out`, and the
    two are not the same number (see below).
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Where to go from here

    The HOMO&ndash;LUMO gap covers one of the properties a dye screen needs.
    Each step below adds another, using the two files from the same ORCA
    job:

    - **Optical gap and absorption spectrum** — run
      `Menu.PLOT_UV_VIS_SPECTRUM` on `tddft.out`. It builds the absorption
      spectrum from the excitation energies and oscillator strengths, which
      shows whether the dye absorbs where the solar spectrum is strongest.
      The optical gap is the lowest-energy bright excitation, so read it
      from the excitation list rather than from the broadened band maximum,
      which neighbouring excitations can shift. Comparing it with the
      HOMO&ndash;LUMO gap above is a first look at how much the excited
      electron and hole are bound to each other, the exciton binding
      energy a good dye keeps small.
    - **Exciton character** — hole&ndash;electron analysis measures how far
      the electron moves from the hole on excitation. A charge-transfer
      excitation (large separation) is what separates charge efficiently
      in a solar cell. It needs both files: the orbitals in
      `optoelectronic.molden` and the excitations in `tddft.out`, which
      Multiwfn asks for partway through. The current
      `Menu.HOLE_ELECTRON_ANALYSIS` preset stops at that prompt, so it has
      to be extended to supply the `tddft.out` path before this step runs
      unattended.
    - **Surface electrostatic potential** — `Menu.QMSA_ESP` on
      `optoelectronic.molden` maps where the dye is electron-rich and
      electron-poor, which governs how it anchors to the semiconductor
      surface and packs with neighbouring molecules.
    - **Ionization potential and electron affinity** — the
      `AnalysisClasses.CDFT` category covers these and other conceptual-DFT
      descriptors. They decide whether the dye's levels line up with the
      semiconductor and electrolyte for charge injection and regeneration.
    - **Run a whole category at once.** Instead of listing individual `Menu`
      members, pass a category from `AnalysisClasses` — e.g.
      `AnalysisClasses.CHARGES` runs *every* supported charge method in one
      call. Useful once you know you want "all of the X", not a specific
      subset.
    - **Caching.** `MultiwfnAnalysis` is cached by default
      (`cached=True`) — re-running the same analysis for a molecule you've
      already processed reuses the stored result instead of re-invoking
      Multiwfn, which matters a lot once you're iterating over thousands of
      candidate dyes rather than one.
    - For the full list of available `Menu` entries and what each one maps
      to, see the enum docstrings in `pymultiwfn.enums.menu`, or the
      [Multiwfn manual](http://sobereva.com/multiwfn) itself for what each
      underlying calculation does chemically.
    """)
    return


if __name__ == "__main__":
    app.run()
