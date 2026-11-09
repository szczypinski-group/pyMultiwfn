import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import marimo as mo

    mo.md(
        r"""
    # QTAIM, ESP, and NCI descriptors for chalcogen-bonded systems

    Benchmarking DFT methods — functionals, dispersion corrections, basis
    sets — against high-level reference data is how computational chemists
    work out which combination actually suits their system. Mehta *et al.*'s
    [CHAL336](https://doi.org/10.1021/acs.jctc.1c00711) benchmark set does
    exactly this for chalcogen bonds (CBs): 336 dimers, split into four
    categories by Lewis-acid character and by whether the interaction runs
    through a σ-hole or a π-hole, with over 1,000 single-point calculations
    per dimer.

    Chalcogen bonds — a group 14–17 atom bound to an electron-withdrawing
    group acting as a Lewis-acid center — show up throughout supramolecular
    and biological chemistry: they steer crystal packing and help stabilize
    protein structure. Understanding *why* a given system binds the way it
    does, beyond just the interaction energy, means characterizing the bond
    critical point at the chalcogen contact (QTAIM), locating the σ- or
    π-hole that drives the interaction (ESP), and mapping the noncovalent
    interaction region itself (NCI) — all of which Multiwfn can compute, and
    `pyMultiwfn` can batch.

    This notebook runs that descriptor trio on a real selenium/sulfur-
    containing molecule (`coord.molden`) to see what each one actually
    reports.
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
      covers the structural and electronic descriptors around the bond —
      QTAIM, ESP, NCI — not the energetics.
    - **QTAIM via bond critical points, not a full topology menu.** The
      source material for this notebook listed `Menu.PLOT_UV_VIS_SPECTRUM`
      here, which computes a UV-Vis spectrum, not anything to do with bond
      critical points — a leftover from a different case study. The actual
      QTAIM descriptor this section needs is `Menu.TOPOLOGY_SEARCH_CPS`,
      which searches for every critical point (nuclear, bond, ring, cage)
      and is what's used below instead.
    - **Frontier orbitals aren't re-derived here.** The paragraph this
      notebook is based on also calls out frontier orbitals as relevant to
      CB character. That's the same ALIE/orbital-composition territory the
      Kraken ligand tutorial already covers (`Menu.QMSA_ALIE`, and
      `Menu.HOLE_ELECTRON_ANALYSIS` from the optoelectronic tutorial) rather
      than something new — add either `Menu` entry to the list below if you
      want it alongside QTAIM/ESP/NCI for a given system.
    - **One molecule, not a non-covalent dimer.** `coord.molden` is a single
      8-atom molecule — a Se(O)₃ group bonded through sulfur to an
      N–C fragment — not two separate monomers held together across a
      gap. The chalcogen-relevant bond critical points below (Se···S and
      S···N) are real QTAIM bond paths *within* this one molecule, not an
      intermolecular contact. For an actual CHAL336-style dimer, the same
      four `Menu` entries apply unchanged — just point the path below at a
      real two-monomer geometry.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## The descriptor batch

    Point this at `coord.molden` (or your own dimer's `.molden`/`.wfn`/
    `.wfx` file) and run the cell below.
    """)
    return


@app.cell
def _():
    from pathlib import Path

    from pymultiwfn import Menu, MultiwfnAnalysis

    dimer_analysis = MultiwfnAnalysis(
        "coord.molden",
        analyses=[
            Menu.VIEW_STRUCTURE,
            Menu.NCI_ANALYSIS,
            Menu.TOPOLOGY_SEARCH_CPS,
            Menu.QMSA_ESP,
        ],
    )
    dimer_analysis.run(work_dir=Path("tutorial_results"))
    return (dimer_analysis,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Reading the results

    `dimer_analysis.results` holds one parsed result per analysis, in the
    order requested.
    """)
    return


@app.cell
def _(dimer_analysis, mo):
    atoms = dimer_analysis.results[0].to_dict()["results"]
    elements = {a["atom_id"]: a["element"] for a in atoms if "atom_id" in a}

    topo_results = dimer_analysis.results[2].to_dict()["results"]
    bond_cps = [cp for cp in topo_results if cp.get("type") == "bond"]
    poincare_hopf = next(r for r in topo_results if "satisfied" in r)
    chalcogen_bonds = [
        cp
        for cp in bond_cps
        if {
            elements.get(cp["bonded_atom1_id"]),
            elements.get(cp["bonded_atom2_id"]),
        }
        & {"Se", "S"}
    ]

    esp_statistics = dimer_analysis.results[3].to_dict()["results"][0]["statistics"]
    nci_result = dimer_analysis.results[1].to_dict()["results"]

    chalcogen_lines = "\n".join(
        f"  - {elements.get(cp['bonded_atom1_id'], '?')}"
        f"{cp['bonded_atom1_id']}–"
        f"{elements.get(cp['bonded_atom2_id'], '?')}{cp['bonded_atom2_id']}"
        for cp in chalcogen_bonds
    )

    mo.md(
        f"""
    - **QTAIM bond critical points**: `Menu.TOPOLOGY_SEARCH_CPS` found
      {len(bond_cps)} (3,-1) bond critical points in total. Filtering for
      ones touching Se or S — the chalcogen atoms this whole analysis is
      about — narrows that down to the bonds that actually matter here:
    {chalcogen_lines}
      The Poincaré–Hopf check is
      {"satisfied" if poincare_hopf["satisfied"] else "**not** satisfied"}
      ({poincare_hopf["nuclear"]} nuclear + {poincare_hopf["bond"]} bond +
      {poincare_hopf["ring"]} ring + {poincare_hopf["cage"]} cage critical
      points), confirming the topology search didn't miss any.
    - **σ-hole proxy** (most positive ESP on the surface):
      {esp_statistics["global_max_kcal_mol"]:.2f} kcal/mol — the surface
      site most likely to act as an electrophilic chalcogen-bond donor,
      from `Menu.QMSA_ESP`'s surface statistics. The most negative value,
      {esp_statistics["global_min_kcal_mol"]:.2f} kcal/mol, marks the most
      nucleophilic site instead — a candidate Lewis-base acceptor region.
    - **NCI analysis**: `Menu.NCI_ANALYSIS` produced
      {nci_result[0]["cube_names"] if nci_result else "no cube files"} —
      the sign(λ₂)ρ and reduced-density-gradient grids that define the NCI
      isosurface. Open those two `.cub` files in VMD (or any other
      isosurface viewer) to see the actual noncovalent interaction region;
      `pyMultiwfn` parses the file names but doesn't attempt to reconstruct
      the 3D isosurface itself.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Scaling this to a benchmark set like CHAL336

    CHAL336 is 336 dimers with 1,000+ single-point calculations each — the
    descriptors above are one call per dimer, not per single-point:

    - **Loop over dimer geometries, one `MultiwfnAnalysis` each.** Just like
      the ligand-screening case study, every dimer writes to its own
      `<input_file_name>.output/`, so results never collide.
    - **Caching matters here too.** `MultiwfnAnalysis` is cached by default
      (`cached=True`) — re-running the batch after adding one more `Menu`
      entry reads the stored JSON for everything already computed instead
      of recomputing all 336 dimers from scratch.
    - **Interaction energies still need sobEDA run separately.** Whatever
      DFT-method comparison this descriptor set is feeding into, remember
      the energetics themselves come from outside this pipeline — the
      descriptors here explain *why* a bond behaves as it does, not how
      strong it is.
    - For the full set of QTAIM, ESP, and NCI-family menus Multiwfn
      exposes (promolecular NCI, IRI, DORI, IGM and its variants among
      them), see the `Menu` docstrings in `pymultiwfn.enums.menu`, or the
      [Multiwfn manual](http://sobereva.com/multiwfn) for what each one
      means chemically.
    """)
    return


if __name__ == "__main__":
    app.run()
