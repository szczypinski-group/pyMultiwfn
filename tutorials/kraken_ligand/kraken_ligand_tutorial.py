import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import marimo as mo

    mo.md(
        r"""
    # Batch descriptor calculation for phosphorus ligand design

    The [Kraken platform](https://doi.org/10.1126/sciadv.abj3528) (Gensch
    *et al.*) tackles ligand design for metal catalysis by treating it as a
    data problem: it catalogues phosphorus ligands, generates conformers for
    each, runs DFT on those conformers to extract physical and electronic
    descriptors, and trains a machine-learning model on the result &mdash;
    over 1,500 ligands' worth. That gives chemists a curated dataset and
    model for *selecting* and *generating* useful phosphorus ligands, rather
    than guessing.

    The part that doesn't scale by hand is step three: computing a
    consistent set of descriptors for every single ligand. `pyMultiwfn` is
    built exactly for that — Multiwfn already computes both the electronic
    and steric descriptors this kind of screening needs, `pyMultiwfn` just
    makes running the same batch across hundreds or thousands of ligands
    trivial instead of tedious.

    This notebook reproduces a small, non-exhaustive proof-of-concept from
    that case study on a real phosphine ligand (`kraken_ligand.molden`, a
    50-atom structure with a phosphorus centre): a steric descriptor pair
    (van der Waals surface area and enclosed volume, plus the grid data
    needed for a Pauli-repulsion-style potential) and two electronic ones
    (V<sub>S,min</sub> &mdash; the most negative electrostatic potential on
    the surface &mdash; and the frontier-orbital-related ALIE surface).
    """
    )
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Before you start

    - **Why only four descriptors.** Several other electronic descriptors
      were considered for this case study and dropped — not because
      Multiwfn can't compute them, but because they either need more than
      one input file, or need P-specific adaptation (e.g. an NMR chemical
      shielding reference for phosphorus needs either a linear-regression
      scaling slope or a second optimized molecule in the same solvent just
      to calibrate against). The four descriptors here only need the one
      ligand geometry/wavefunction file.
    - **The Pauli potential is not computed here.** It stops at an
      interactive prompt in plain Multiwfn rather than running
      non-interactively. The surface analyses below export their grid data
      to a `.cube` file (`surf.cub`), which is what a *second*, separate
      analysis would read for a Pauli-repulsion-style potential.
      V<sub>S,min</sub>, by contrast, comes straight out of
      `Menu.QMSA_ESP`'s surface analysis with no extra step.
    - **Steric and electronic descriptors in the same call.** Unlike the
      optoelectronic case study, every descriptor here reads from the same
      plain `.molden`/`.wfn`/`.wfx` geometry — no excited-state output
      needed, so this one runs on literally any ground-state wavefunction
      file you already have.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## The ligand descriptor batch

    Point this at `kraken_ligand.molden` (or any other `.molden`/`.wfn`/
    `.wfx` ligand file) and run the cell below.
    """)
    return


@app.cell
def _(mo):
    import json

    from pymultiwfn import Menu, MultiwfnAnalysis

    # Paths are relative to this notebook's folder, wherever marimo is
    # launched from.
    notebook_dir = mo.notebook_dir()

    ligand_analysis = MultiwfnAnalysis(
        notebook_dir / "kraken_ligand.molden",
        analyses=[
            Menu.VIEW_STRUCTURE,
            Menu.SURFACE_AREA_VOLUME,
            Menu.QMSA_ALIE,
            Menu.QMSA_ESP,
        ],
    )
    ligand_analysis.run(work_dir=notebook_dir / "tutorial_results")
    return json, ligand_analysis, notebook_dir


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Reading the results

    Every run writes one JSON summary per molecule, inside
    `tutorial_results/<input_file_name>.output/`. It is rebuilt from that
    folder on every run, so it holds the same data whether an analysis was
    just computed or loaded from the cache — read results from there.
    """)
    return


@app.cell
def _(json, ligand_analysis, mo, notebook_dir):
    ligand_summary_path = (
        notebook_dir
        / "tutorial_results"
        / f"{ligand_analysis.input_file.name}.output"
        / f"{ligand_analysis.input_file.stem}.json"
    )
    ligand_results = {
        name: entry["parsed"]["results"]
        for name, entry in json.loads(ligand_summary_path.read_text())[
            "analyses"
        ].items()
    }
    geometry = ligand_results["SURFACE_AREA_VOLUME"][0]["geometry"]
    esp_result = ligand_results["QMSA_ESP"][0]
    alie_result = ligand_results["QMSA_ALIE"][0]

    # The P lone pair: the ESP surface minimum closest to phosphorus.
    bohr_to_angstrom = 0.529177210903
    phosphorus = next(
        atom
        for atom in ligand_results["VIEW_STRUCTURE"]
        if atom.get("element") == "P"
    )
    p_xyz = [
        phosphorus[axis] * bohr_to_angstrom
        for axis in ("x_bohr", "y_bohr", "z_bohr")
    ]

    def distance_to_p(point: dict) -> float:
        return (
            sum(
                (point[f"{axis}_angstrom"] - p) ** 2
                for axis, p in zip("xyz", p_xyz, strict=True)
            )
            ** 0.5
        )

    p_lone_pair = min(esp_result["minima"], key=distance_to_p)

    mo.md(
        f"""
    - **Van der Waals volume and surface area** (one call, both numbers):
      {geometry["volume_angstrom3"]:.2f} Å<sup>3</sup> enclosed by
      {geometry["area_angstrom2"]:.2f} Å<sup>2</sup> of surface — the
      steric descriptor pair, straight from `Menu.SURFACE_AREA_VOLUME`.
    - **V<sub>S,min</sub> at the phosphorus lone pair**:
      {p_lone_pair["value_kcal_mol"]:.2f} kcal/mol, at the ESP surface
      minimum {distance_to_p(p_lone_pair):.2f} Å from P
      (P{phosphorus["atom_id"]}). For ligand design this is the
      V<sub>S,min</sub> that matters: it describes the lone pair that
      donates to the metal.
    - **Global V<sub>S,min</sub>** (most negative ESP anywhere on the
      surface): {esp_result["statistics"]["global_min_kcal_mol"]:.2f}
      kcal/mol. In this ligand it is *not* at phosphorus, so always check
      where a minimum sits before using it as a descriptor.
    - **ALIE minimum**: {alie_result["statistics"]["global_min_au"]:.4f}
      a.u. — the most easily ionized electrons on the surface.
    - The full ALIE and ESP surface analyses also carry every surface
      minimum/maximum (with positions) and the area/polarity statistics —
      not just the headline numbers above.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Scaling this to a real ligand library

    Kraken's own dataset is 1,500+ ligands; a single call per ligand, same
    four `Menu` entries, is what makes that tractable:

    - **Point it at a directory of ligands.** Loop over your conformer
      files and construct one `MultiwfnAnalysis` per file — each writes to
      its own `<input_file_name>.output/`, so results never collide between
      ligands even when run from the same working directory.
    - **Caching matters at this scale.** `MultiwfnAnalysis` is cached by
      default (`cached=True`): re-running the same descriptor for a ligand
      you've already processed reads the stored JSON instead of
      re-invoking Multiwfn, which is the difference between reprocessing a
      library of 1,500 ligands taking minutes versus hours if you need to
      add one more descriptor later.
    - **The `.cube` file produced here is reusable.** The `surf.cub` file
      written alongside the ESP/ALIE surface analyses is the input the
      Pauli-potential-style follow-up calculation needs — feed it into a
      second `MultiwfnAnalysis` as its own input file rather than
      re-deriving it from scratch.
    - For the full set of electronic and steric descriptors Multiwfn can
      compute (beyond this proof-of-concept's four), see the `Menu`
      docstrings in `pymultiwfn.enums.menu`, or the
      [Multiwfn manual](http://sobereva.com/multiwfn) for what each one
      means chemically.
    """)
    return


if __name__ == "__main__":
    app.run()
