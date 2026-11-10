import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import marimo as mo

    mo.md(
        r"""
    # NCI analysis for chalcogen-bonded systems

    Benchmarking DFT methods — functionals, dispersion corrections, basis
    sets — against high-level reference data is how computational chemists
    work out which combination actually suits their system. Mehta *et al.*'s
    [CHAL336](https://doi.org/10.1021/acs.jctc.1c00711) benchmark set does
    exactly this for chalcogen bonds (CBs): 336 dimers, split into four
    categories by Lewis-acid character and by whether the interaction runs
    through a σ-hole or a π-hole, with over 1,000 single-point calculations
    per dimer.

    In a chalcogen bond, a group 16 atom (O, S, Se, Te) acts as the
    Lewis acid, through an electron-deficient σ- or π-hole, towards a
    Lewis base. These contacts show up throughout supramolecular and
    biological chemistry: they steer crystal packing and help stabilize
    protein structure. Understanding *why* a given system binds the way it
    does, beyond just the interaction energy, starts with seeing where the
    noncovalent interaction actually is and how strong it is — which is
    what NCI (non-covalent interaction, or reduced density gradient)
    analysis shows, and what Multiwfn can compute and `pyMultiwfn` can
    batch.

    This notebook runs NCI analysis on a real selenium/sulfur-containing
    molecule (`coord.molden`) and reads off the interaction regions it
    finds.
    """
    )
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Before you start

    - **Interaction energy is not computed here.** For a DFT benchmark set
      like CHAL336 the headline quantity is the interaction energy itself —
      and Multiwfn *can* get at it, through the sobEDA method. But sobEDA
      needs an external interfacing script that `pyMultiwfn` does not
      currently wrap, so it's out of scope for this notebook. What follows
      maps the interaction, not its energy.
    - **What NCI analysis produces.** Multiwfn evaluates the reduced
      density gradient (RDG) and sign(λ₂)ρ on a grid around the molecule
      and writes them to a scatter file (`output.txt`) and two cube files
      (`func1.cub` for sign(λ₂)ρ, `func2.cub` for RDG). `pyMultiwfn` reads
      the scatter file and reports each interaction region ("spike") it
      finds: low RDG at low density. The sign of sign(λ₂)ρ says whether a
      region is attractive (negative), van der Waals (near zero) or
      repulsive (positive), and its magnitude says how strong it is.
    - **One molecule, not a non-covalent dimer.** `coord.molden` is a single
      8-atom molecule — a Se(O)₃ group bonded through sulfur to an
      N–C fragment — not two separate monomers held together across a
      gap. The Se···S contact found below is *within* this one molecule,
      not an intermolecular contact. For an actual CHAL336-style dimer, the
      same `Menu` entry applies unchanged — just point the path below at a
      real two-monomer geometry.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Running the NCI analysis

    Point this at `coord.molden` (or your own dimer's `.molden`/`.wfn`/
    `.wfx` file) and run the cell below.
    """)
    return


@app.cell
def _(mo):
    import json

    from pymultiwfn import Menu, MultiwfnAnalysis

    # Paths are relative to this notebook's folder, wherever marimo is
    # launched from.
    notebook_dir = mo.notebook_dir()

    dimer_analysis = MultiwfnAnalysis(
        notebook_dir / "coord.molden",
        analyses=[Menu.NCI_ANALYSIS],
    )
    dimer_analysis.run(work_dir=notebook_dir / "tutorial_results")
    return dimer_analysis, json, notebook_dir


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Reading the results

    Every run writes one JSON summary per molecule, inside
    `tutorial_results/<input_file_name>.output/`. It is rebuilt from that
    folder on every run, so it holds the same data whether the analysis was
    just computed or loaded from the cache — read results from there.
    """)
    return


@app.cell
def _(dimer_analysis, json, mo, notebook_dir):
    dimer_summary_path = (
        notebook_dir
        / "tutorial_results"
        / f"{dimer_analysis.input_file.name}.output"
        / f"{dimer_analysis.input_file.stem}.json"
    )
    dimer_summary = json.loads(dimer_summary_path.read_text())
    nci = dimer_summary["analyses"]["NCI_ANALYSIS"]["parsed"]["results"][0]

    spike_rows = []
    for spike in nci["spikes"]:
        contacts = ", ".join(
            f"{atom['atom']} ({atom['distance_angstrom']:.2f} Å)"
            for atom in spike["nearest_atoms"] or []
        )
        spike_rows.append(
            f"| {spike['interaction_type']} "
            f"| {spike['sign_lambda2_rho_at_min']:+.4f} "
            f"| {spike['rdg_min']:.3f} "
            f"| {spike['n_points']} "
            f"| {contacts} |"
        )
    spike_table = "\n".join(
        [
            "| Type | sign(λ₂)ρ at lowest RDG (a.u.) | Lowest RDG "
            "| Grid points | Nearest atoms |",
            "|---|---|---|---|---|",
            *spike_rows,
        ]
    )
    grid = " × ".join(str(n) for n in nci["grid_points_xyz"])

    mo.vstack(
        [
            mo.md(
                f"""
    `Menu.NCI_ANALYSIS` evaluated RDG and sign(λ₂)ρ on a {grid} grid
    ({nci["grid_points_total"]:,} points, {nci["grid_spacing_bohr"][0]:.3f}
    Bohr spacing) for {nci["total_electrons"]:.0f} electrons, and found
    these interaction regions:
    """
            ),
            mo.md(spike_table),
            mo.md(
                f"""
    The attractive region sits between Se1 and S5: that is the Se···S
    chalcogen contact. Its sign(λ₂)ρ of about −0.04 a.u. is the electron
    density at the contact, well above typical van der Waals values (below
    about 0.01 a.u.), so this is a strong, partly covalent interaction. The
    weaker repulsive regions alongside it are the steric part of the same
    contact.

    To see the interaction region in 3D, open `{nci["rdg_cube"]}` (RDG) as
    an isosurface coloured by `{nci["sign_lambda2_rho_cube"]}`
    (sign(λ₂)ρ) in VMD or any other cube viewer. Both are in
    `{dimer_summary_path.parent.relative_to(notebook_dir)}`.
    """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Scaling this to a benchmark set like CHAL336

    CHAL336 is 336 dimers with 1,000+ single-point calculations each — the
    NCI analysis above is one call per dimer, not per single-point:

    - **Loop over dimer geometries, one `MultiwfnAnalysis` each.** Just like
      the ligand-screening case study, every dimer writes to its own
      `<input_file_name>.output/`, so results never collide.
    - **Caching matters here too.** `MultiwfnAnalysis` is cached by default
      (`cached=True`) — re-running the batch after adding one more `Menu`
      entry reuses everything already computed instead of recomputing all
      336 dimers from scratch. Read the results from each molecule's JSON
      summary, as above, so cached and freshly computed analyses look the
      same.
    - **Interaction energies still need sobEDA run separately.** Whatever
      DFT-method comparison this feeds into, remember the energetics
      themselves come from outside this pipeline — NCI analysis shows
      *where* the interaction is and how it is distributed, not how strong
      it is energetically.
    - For the other weak-interaction menus Multiwfn exposes (promolecular
      NCI, IRI, DORI, IGM and its variants among them), and for QTAIM and
      ESP descriptors, see the `Menu` docstrings in
      `pymultiwfn.enums.menu`, or the
      [Multiwfn manual](http://sobereva.com/multiwfn) for what each one
      means chemically.
    """)
    return


if __name__ == "__main__":
    app.run()
