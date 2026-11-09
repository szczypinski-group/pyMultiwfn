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
    *different* interactive menu by hand for every molecule:

    - **HOMO&ndash;LUMO gap** → main function `0` (view structure / orbitals)
    - **Optical gap** → a TD-DFT output, main function `11` → `3` (UV&ndash;Vis
      absorption spectrum)
    - **Nature of the transition / exciton analysis** → main function `18`
      (hole&ndash;electron analysis)

    `pyMultiwfn` turns that into one call: give it a wavefunction file and
    the list of `Menu` entries you want, and it drives Multiwfn through all
    of them, parses every stdout into structured data, and writes one JSON
    summary per molecule. This notebook walks through that, first with a
    tiny bundled example so you can see the mechanics work, then with a real
    ORCA TD-DFT output (`tddft.out`) from the optoelectronic-screening
    pattern in the case study.
    """
    )
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Before you start

    A few things that matter for every analysis below, not just this one:

    - **File type matters.** Multiwfn doesn't accept every input for every
      calculation — ground-state properties (charges, bond orders, basic
      structure) work fine from a `.molden`/`.wfn`/`.wfx`/`.fchk` file, but
      anything that needs **excited-state data** (UV&ndash;Vis spectra,
      hole&ndash;electron analysis) needs a Gaussian or ORCA `.out`/`.log`
      file instead, because only those contain the excitation-state block.
      `pyMultiwfn` checks the extension up front and raises
      `InvalidInputFileError` immediately if it won't work, rather than
      letting Multiwfn fail halfway through a calculation.
    - **A TD-DFT `.out` file gets you the spectrum, not the density.**
      `tddft.out` below is a genuine ORCA 6.1.0 TD-DFT/TDA run with real
      excitation energies and oscillator strengths, so
      `Menu.PLOT_UV_VIS_SPECTRUM` works on it directly. But Multiwfn only
      ever reads the *geometry* and the excitation list out of a plain
      `.out` log — not the full MO coefficient/basis-function data a
      density-based analysis needs. Point `Menu.QMSA_ESP` or
      `Menu.HOLE_ELECTRON_ANALYSIS` at this same file and Multiwfn either
      refuses outright ("the input file does not contain basis function
      information") or silently computes nonsense (a bare-nuclear-charge
      ESP with no electron screening, off by orders of magnitude). For
      those, export the actual wavefunction from the same ORCA job with
      `orca_2mkl jobname -molden` and feed *that* `.molden` file in
      instead — exactly how `kraken_ligand.molden` and `coord.molden` in
      the other two tutorials were made in the first place.
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
    ## 1. Quickstart: one call, structured output
    """)
    return


@app.cell
def _():
    from pathlib import Path

    from pymultiwfn import Menu, MultiwfnAnalysis

    return Menu, MultiwfnAnalysis, Path


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The repository ships a small example wavefunction
    (`tests/test_data/M062X_TZVPP_D30.molden`) purely so this cell runs for
    you out of the box. Swap in any of your own `.molden`/`.wfn`/`.wfx` files
    to try it on real data — the pattern doesn't change.
    """)
    return


@app.cell
def _(Menu, MultiwfnAnalysis, Path):
    quickstart_input = Path("tests/test_data/M062X_TZVPP_D30.molden")

    quickstart = MultiwfnAnalysis(
        quickstart_input,
        analyses=[Menu.VIEW_STRUCTURE, Menu.HIRSHFELD_CHARGE],
    )
    quickstart.run(work_dir=Path("tutorial_results"))

    # One MultiwfnResult per analysis, in the order you requested them.
    quickstart.results[1].to_dict()
    return (quickstart,)


@app.cell(hide_code=True)
def _(mo, quickstart):
    mo.md(rf"""
    `quickstart.run()` drove Multiwfn through *both* requested menus in one
    process, and `quickstart.results` holds one parsed
    [`MultiwfnResult`][result] per analysis, in the order you asked for
    them. Everything it produced — including a combined
    `{quickstart.input_file.stem}.json` — is sitting in
    `tutorial_results/{quickstart.input_file.name}.output/` for you to
    inspect or hand off to another tool.

    [result]: https://szczypinski-group.github.io/pyMultiwfn/
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 2. The optoelectronic screening case study

    This is the pattern from the dye-screening case study itself, run
    against a real TD-DFT output (`tddft.out`). Unlike the other two
    tutorials, this batch only asks for the two `Menu` entries that
    actually work from a plain TD-DFT `.out` file — see "Before you start"
    above for why `Menu.QMSA_ESP` and `Menu.HOLE_ELECTRON_ANALYSIS` are
    left out here rather than quietly returning wrong numbers.
    """)
    return


@app.cell
def _(Menu, MultiwfnAnalysis, Path):
    dye_analysis = MultiwfnAnalysis(
        "tddft.out",
        analyses=[Menu.VIEW_STRUCTURE, Menu.PLOT_UV_VIS_SPECTRUM, Menu.QMSA_ESP],
    )
    dye_analysis.run(work_dir=Path("tutorial_results"))
    return (dye_analysis,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Reading the results

    `dye_analysis.results` holds one parsed result per analysis, in the
    order requested.
    """)
    return


@app.cell
def _(dye_analysis, mo):
    atoms = {
        a["atom_id"]
        for a in dye_analysis.results[0].to_dict()["results"]
        if "atom_id" in a
    }
    spectrum_extrema = dye_analysis.results[1].to_dict()["results"][0]["extrema"]
    peak = max(spectrum_extrema, key=lambda e: e["value"])

    mo.md(
        f"""
    - **Structure check**: `Menu.VIEW_STRUCTURE` parsed {len(atoms)} atoms
      from the TD-DFT job's own geometry — the same sanity check you'd want
      before trusting any descriptor computed on it.
    - **Optical gap**: `Menu.PLOT_UV_VIS_SPECTRUM` built the simulated
      absorption curve from every excitation energy and oscillator strength
      in the file, and its strongest peak lands at
      {peak["x"]:.1f} nm — this *is* the optical gap, read directly off the
      broadened TD-DFT spectrum rather than picked by hand from the raw
      state list.
    - **What you'd add with the wavefunction export**: V<sub>S,min</sub> /
      V<sub>S,max</sub> (`Menu.QMSA_ESP`) and the hole/electron centroid
      distance for exciton character (`Menu.HOLE_ELECTRON_ANALYSIS`) both
      need the actual MO coefficients — export `tddft.molden` from the same
      ORCA job with `orca_2mkl` and add those two `Menu` entries pointed at
      *that* file instead.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 3. Where to go from here

    - **Run a whole category at once.** Instead of listing individual `Menu`
      members, pass a category from `AnalysisClasses` — e.g.
      `AnalysisClasses.CHARGES` runs *every* supported charge method in one
      call. Useful once you know you want "all of the X", not a specific
      subset.
    - **Ionization potential / electron affinity and other CDFT
      descriptors** live under the `AnalysisClasses.CDFT` category if your
      screening needs those too.
    - **Caching.** `MultiwfnAnalysis` is cached by default
      (`cached=True`) — re-running the same analysis for a molecule you've
      already processed loads the stored result instead of re-invoking
      Multiwfn, which matters a lot once you're iterating over thousands of
      structures rather than one.
    - For the full list of available `Menu` entries and what each one maps
      to, see the enum docstrings in `pymultiwfn.enums.menu`, or the
      [Multiwfn manual](http://sobereva.com/multiwfn) itself for what each
      underlying calculation does chemically.
    """)
    return


if __name__ == "__main__":
    app.run()
