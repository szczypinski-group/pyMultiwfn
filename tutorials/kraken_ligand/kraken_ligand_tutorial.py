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
    - **Two of these can't run to completion on their own.** The Pauli
      potential and V<sub>S,min</sub> calculations both stop at an
      interactive prompt in plain Multiwfn, rather than running
      non-interactively. The workaround, and what `pyMultiwfn` does under
      the hood, is to export the surface grid data to a `.cube` file
      instead — that `.cube` file is what a *second*, separate analysis
      would then read to get the final number, rather than trying to force
      one interactive calculation to behave non-interactively.
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
def _():
    from pathlib import Path

    from pymultiwfn import Menu, MultiwfnAnalysis

    ligand_analysis = MultiwfnAnalysis(
        "kraken_ligand.molden",
        analyses=[
            Menu.VIEW_STRUCTURE,
            Menu.SURFACE_AREA_VOLUME,
            Menu.QMSA_ALIE,
            Menu.QMSA_ESP,
        ],
    )
    ligand_analysis.run(work_dir=Path("tutorial_results"))
    return (ligand_analysis,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Reading the results

    `ligand_analysis.results` holds one parsed result per analysis, in the
    order requested.
    """)
    return


@app.cell
def _(ligand_analysis, mo):
    area_volume_result = ligand_analysis.results[1].to_dict()["results"][0]
    esp_result = ligand_analysis.results[3].to_dict()["results"][0]
    geometry = area_volume_result["geometry"]

    mo.md(
        f"""
    - **Van der Waals volume and surface area** (one call, both numbers):
      {geometry["volume_angstrom3"]:.2f} Å<sup>3</sup> enclosed by
      {geometry["area_angstrom2"]:.2f} Å<sup>2</sup> of surface — the
      steric descriptor pair, straight from `Menu.SURFACE_AREA_VOLUME`.
    - **V<sub>S,min</sub>** (most negative ESP on the surface):
      {esp_result["statistics"]["global_min_kcal_mol"]:.2f} kcal/mol — read
      off `Menu.QMSA_ESP`'s surface statistics, no separate `.cube` step
      needed for this particular quantity.
    - The full ALIE and ESP surface analyses (`results[2]` and `results[3]`)
      also carry every surface minimum/maximum and the area/polarity
      statistics — not just the one headline number above.
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
