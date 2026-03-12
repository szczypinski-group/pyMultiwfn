# pymultiwfn — Use Case Guide

## What pymultiwfn does

pymultiwfn automates the Multiwfn wavefunction analysis program. Multiwfn is an
interactive command-line tool that expects keystroke sequences typed into menus
and submenus. pymultiwfn translates analysis requests into those keystroke
sequences, feeds them to the Multiwfn executable via stdin, captures stdout, and
parses the text output into structured Python objects. Parsed results are
persisted to a per-molecule JSON file so that identical analyses are not re-run.


## Class map

```
┌─────────────────────────────────────────────────────────────────┐
│                      MultiwfnAnalysis                           │
│  The orchestrator. Holds a queue of Menu items, a list of       │
│  completed MultiwfnResult objects, and a ResultStore for        │
│  JSON persistence.                                              │
│                                                                 │
│  .analyses : list[Menu]          ← what to run                  │
│  .results  : list[MultiwfnResult]← parsed output per analysis   │
│  .jobs     : list[MultiwfnJob]   ← raw execution records        │
│  ._store   : ResultStore         ← JSON cache (lazy-init)       │
├─────────────────────────────────────────────────────────────────┤
│  For each Menu in .analyses, .run() either:                     │
│    • loads a cached result from the JSON store, or              │
│    • creates a MultiwfnJob, executes it, parses stdout into     │
│      a MultiwfnResult, and writes that result to the store.     │
└────────────┬──────────────────────────────────┬─────────────────┘
             │                                  │
             ▼                                  ▼
┌────────────────────────┐       ┌──────────────────────────────┐
│      MultiwfnJob       │       │       MultiwfnResult         │
│                        │       │                              │
│  Runs Multiwfn once.   │       │  Container for one analysis. │
│  Accepts an input file │       │                              │
│  and a single Menu.    │       │  .analysis : Menu            │
│  Translates the Menu   │       │  .result   : list[Parsed...] │
│  into a keystroke      │       │                              │
│  sequence, writes it   │       │  .parse(stdout) looks up the │
│  to a temp .inp file,  │       │  right OutputParser subclass │
│  and pipes it to the   │       │  via ParserRoute, calls      │
│  Multiwfn subprocess.  │       │  parse_for_result(), and     │
│                        │       │  stores the returned objects │
│  After execution:      │       │  in .result.                 │
│  .stdout  : str        │       │                              │
│  .stderr  : str        │       │  .to_dict() serialises       │
│  .success : bool       │       │  everything to plain dicts.  │
└────────────────────────┘       └──────────────┬───────────────┘
                                                │
             ┌──────────────────────────────────┘
             ▼
┌──────────────────────────────────────────────────────────────────┐
│                         OutputParser                             │
│  Base class. Each subclass (ChargeParser, BondOrderParser, etc.) │
│  owns a parse_for_result(analysis, stdout) classmethod.          │
│                                                                  │
│  parse_for_result consults an internal case_map dict keyed by    │
│  Menu to decide which static parsing methods to call, then       │
│  flattens the results through _collect().                        │
│                                                                  │
│  ParserRoute.ROUTE_TABLE maps every supported Menu member to     │
│  the correct OutputParser subclass. MultiwfnResult.parse() uses  │
│  this table as its sole dispatch mechanism.                      │
└──────────────────────────────────────────────────────────────────┘

┌────────────────────────┐       ┌──────────────────────────────┐
│        Multiwfn        │       │       ResultStore            │
│                        │       │                              │
│  Locates the Multiwfn  │       │  Manages <file>.json on      │
│  executable on disk.   │       │  disk. Keyed by Menu.name.   │
│  Passed to             │       │                              │
│  MultiwfnJob so the    │       │  .has_result(menu) → bool    │
│  subprocess knows      │       │  .store(result)   → write    │
│  which binary to call. │       │  .get_result(menu)→ dict     │
└────────────────────────┘       └──────────────────────────────┘

┌────────────────────────┐       ┌──────────────────────────────┐
│         Menu           │       │      AnalysisClasses         │
│                        │       │                              │
│  Enum. Each member's   │       │  Enum of pre-defined lists   │
│  value is a tuple of   │       │  of Menu items.              │
│  strings — the exact   │       │                              │
│  keystrokes to type    │       │  CHARGES, BOND_ORDERS,       │
│  into Multiwfn.        │       │  TOPOLOGY, SPECTRA, etc.     │
│                        │       │                              │
│  .get_sequence()       │       │  Pass one to add_menu() to   │
│  returns that tuple.   │       │  queue an entire category.   │
└────────────────────────┘       └──────────────────────────────┘
```


## How a job runs end-to-end

This section traces what happens from the moment you call `analysis.run()` to
the moment you read results back.

**1. Queue analyses.** You create a `MultiwfnAnalysis` with an input file and
one or more `Menu` members. Each `Menu` member encodes the literal keystrokes
Multiwfn expects — for example `Menu.HIRSHFELD_CHARGE` stores
`("7", "1", "1", "n", "0")`, meaning "enter main menu 7, pick option 1,
choose built-in atomic densities, decline file output, return to main menu."

**2. Check the cache.** When `run()` starts iterating over the queued `Menu`
items, it first checks the `ResultStore`. If the JSON file already contains
a parsed entry for that `Menu` name, the analysis is skipped and no subprocess
is launched.

**3. Build and execute a job.** For uncached analyses, `_create_and_run`
creates a `MultiwfnJob`. The job translates the `Menu` value into a command
list, writes that list to a temporary `.inp` file, and launches the Multiwfn
binary as a subprocess with the `.inp` file piped to stdin and the wavefunction
file passed as a positional argument. Stdout and stderr are captured in a
`MultiwfnJobOutcome` dataclass.

**4. Parse stdout.** After the subprocess finishes, `_create_and_run` creates
a `MultiwfnResult` bound to the same `Menu` and calls its `parse(stdout)`
method. Inside `parse`:

   - `ParserRoute.ROUTE_TABLE` is consulted to find the correct `OutputParser`
     subclass for this `Menu` member. For example, every charge-related `Menu`
     maps to `ChargeParser`, every bond-order `Menu` maps to
     `BondOrderParser`, and so on.
   - The subclass's `parse_for_result(analysis, stdout)` classmethod is called.
     This method consults an internal `case_map` dictionary to decide which
     static regex-based helper methods apply to the specific `Menu` variant,
     then calls `_collect()` to run them all and flatten the results into a
     single list.
   - The returned `list[ParsedMultiwfnResult]` is stored in
     `MultiwfnResult.result`.

**5. Persist to JSON.** The `MultiwfnResult` is handed to
`ResultStore.store()`, which calls `to_dict()` on the result (serialising each
`ParsedMultiwfnResult` dataclass via `dataclasses.asdict`) and writes it into
the JSON file under the `Menu` member's name. The file is saved to disk
immediately.

**6. Access results.** After `run()` completes, `analysis.results` contains
one `MultiwfnResult` per queued `Menu`. Each `MultiwfnResult.result` is a flat
list of typed dataclasses — `Charge`, `BondOrder`, `CriticalPoint`,
`Spectrum`, etc. — that you can iterate, filter by type, or serialise further.


## Parsed result types

Every parser produces instances of `ParsedMultiwfnResult` subclasses. These
are plain `@dataclass` objects with a `to_dict()` method. The subclass you get
depends on which `Menu` was analysed:

| Parser class                | Menu category          | Result dataclasses produced                                |
|-----------------------------|------------------------|------------------------------------------------------------|
| `ChargeParser`              | Menu 7 — charges       | `Charge`, `Dipole`                                         |
| `OrbitalCompositionParser`  | Menu 8 — orbitals      | `OrbitalComponent`, `OxidationState`                       |
| `BondOrderParser`           | Menu 9 — bond orders   | `BondOrder`, `Valence`, `MultiCenterBondOrder`, `BondOrderDecomposition` |
| `CriticalPointParser`       | Menu 2 — topology      | `CriticalPoint`, `BondPath`                                |
| `DOSParser`                 | Menu 10 — DOS          | `DensityOfStates`, `OrbitalEnergy`                         |
| `SpectrumParser`            | Menu 11 — spectra      | `Spectrum`, `Transition`, `Color`                          |
| `SurfaceParser`             | Menu 12 — surfaces     | `SurfaceAnalysis`, `SurfaceExtremum`                       |
| `FuzzySpaceParser`          | Menu 15 — fuzzy space  | `FuzzyAtomicProperty`, `DelocalizationIndex`, `AromaticityIndex` |
| `BasinParser`               | Menu 17 — basins       | `Basin`, `Charge`                                          |
| `ExcitationParser`          | Menu 18 — excitations  | `HoleElectron`, `ChargeTransfer`, `DeltaR`, `LambdaIndex`  |
| `WeakInteractionParser`     | Menu 20 — NCI / IGM    | `WeakInteraction`                                          |
| `EDAParser`                 | Menu 21 — EDA          | `EnergyDecompositionAnalysis`, `DispersionContribution`    |
| `CDFTParser`                | Menu 22 — CDFT         | `Reactivity`, `CondensedFukui`, `DualDescriptor`           |
| `PolarizabilityParser`      | Menu 24 — polarizability | `Polarizability`                                         |
| `AromaticityParser`         | Menu 25 — aromaticity  | `Aromaticity`, `NICSScan`                                  |
| `WavefunctionParser`        | Menu 6 — wavefunction  | `Orbital`                                                  |
| `CubeParser`                | Menu 5 — cubes         | `Cube`                                                     |
| `UtilityParser`             | Menu 100/200/300       | `BondLength`, `BondAngle`, `DihedralAngle`, `DipoleMoment`, `QuadrupoleMoment`, `CoordinationNumber`, `BLA_BOA` |


## How parser dispatch works

Not every `Menu` member within a parser class calls the same regex helpers. For
example, `BondOrderParser` handles both pairwise bond orders and multicenter
bond orders, which appear in completely different stdout formats. The dispatch
works in two layers:

**Layer 1 — `ParserRoute.ROUTE_TABLE`.** This is a flat dictionary mapping
every supported `Menu` member to the `OutputParser` subclass responsible for
it. `MultiwfnResult.parse()` does a single lookup here.

**Layer 2 — `case_map` inside `parse_for_result`.** Each subclass defines a
dictionary mapping specific `Menu` members to lists of parsing callables. If
the `Menu` is found in the map, only those callables run. Otherwise a default
list runs. The shared `_collect(stdout, parsers)` method iterates the list,
calls each function with `stdout`, and accumulates results — `None` return
values are skipped, single objects are appended, and lists are extended.

Example from `BondOrderParser`:

```python
case_map = {
    Menu.MULTICENTER_BOND_ORDER:        [cls.parse_multicenter],
    Menu.MULTICENTER_BOND_ORDER_NAO:    [cls.parse_multicenter],
    Menu.MULLIKEN_BOND_ORDER_DECOMPOSE: [cls.parse_decomposition],
    Menu.WIBERG_DECOMPOSITION:          [cls.parse_decomposition],
}
default = [cls.parse, cls.parse_valence]
```

When `Menu.MAYER_BOND_ORDER` arrives, it is not in `case_map`, so the default
runs: `parse()` extracts pairwise bond orders, then `parse_valence()` extracts
total and free valences. When `Menu.MULTICENTER_BOND_ORDER` arrives, only
`parse_multicenter()` runs.


## The JSON store

`ResultStore` creates a file named `<input_file>.json` in the working
directory. Its structure:

```json
{
  "input_file": "benzene.wfn",
  "analyses": {
    "HIRSHFELD_CHARGE": {
      "parsed": {
        "analysis": "HIRSHFELD_CHARGE",
        "results": [
          {"atom_id": 1, "charge": 0.032},
          {"atom_id": 2, "charge": -0.032},
          {"x": 0.0, "y": 0.0, "z": 0.12, "total": 0.12}
        ]
      },
      "timestamp": "2026-03-10T14:22:01.123456"
    },
    "MAYER_BOND_ORDER": {
      "parsed": { ... },
      "timestamp": "2026-03-10T14:22:03.789012"
    }
  }
}
```

Each key under `"analyses"` is the `Menu` member name. The `"parsed"` value is
the output of `MultiwfnResult.to_dict()`. The file is rewritten on every new
analysis. On the next run, `has_result()` checks whether a key exists and
skips re-execution if `cached=True` (the default).


## Usage examples

### Run a single analysis

```python
from pathlib import Path
from pymultiwfn import MultiwfnAnalysis, Multiwfn, Menu

mwfn = Multiwfn(exe_path=Path("/opt/multiwfn/Multiwfn"))

analysis = MultiwfnAnalysis("benzene.wfn", analyses=Menu.HIRSHFELD_CHARGE)
analysis.run(multiwfn=mwfn, work_dir=Path("./output"))

result = analysis.results[0]          # MultiwfnResult for HIRSHFELD_CHARGE
for item in result.result:
    print(item)                        # Charge(...) and Dipole(...) objects
```

### Queue several analyses

```python
analysis = MultiwfnAnalysis("benzene.wfn")
analysis.add_menu(Menu.HIRSHFELD_CHARGE)
analysis.add_menu(Menu.MAYER_BOND_ORDER)
analysis.add_menu(Menu.TOPOLOGY_SEARCH_CPS)
analysis.run(multiwfn=mwfn, timeout=120, work_dir=Path("./output"))

# results[0] → charges,  results[1] → bond orders,  results[2] → critical points
```

### Queue an entire category

```python
from pymultiwfn.enums.analyses import AnalysisClasses

analysis = MultiwfnAnalysis("benzene.wfn")
analysis.add_menu(AnalysisClasses.CHARGES)   # queues all 19 charge methods
analysis.run(multiwfn=mwfn, work_dir=Path("./output"))
```

### Filter results by type

```python
from pymultiwfn.analysis.result import Charge, Dipole

for r in analysis.results:
    charges = [x for x in r.result if isinstance(x, Charge)]
    dipoles = [x for x in r.result if isinstance(x, Dipole)]
```

### Re-run with caching

```python
# Second run loads from benzene.wfn.json instead of re-executing Multiwfn.
analysis2 = MultiwfnAnalysis("benzene.wfn", analyses=Menu.HIRSHFELD_CHARGE)
analysis2.run(multiwfn=mwfn, work_dir=Path("./output"))
# No subprocess launched — result loaded from cache.
```

### Disable caching

```python
analysis = MultiwfnAnalysis("benzene.wfn", analyses=Menu.HIRSHFELD_CHARGE, cached=False)
analysis.run(multiwfn=mwfn, work_dir=Path("./output"))
# Always re-runs Multiwfn, even if the JSON file has a stored result.
```

### Use ResultStore directly

```python
from pymultiwfn.analysis.result import ResultStore, MultiwfnResult

store = ResultStore(input_file=Path("benzene.wfn"), work_dir=Path("./output"))

# Parse raw stdout from some other source and persist it:
result = store.store_from_stdout(Menu.MAYER_BOND_ORDER, raw_stdout_string)

# Check what is stored:
print(store.has_result(Menu.MAYER_BOND_ORDER))  # True
print(store.get_result(Menu.MAYER_BOND_ORDER))  # dict with "parsed" and "timestamp"
```

### Use MultiwfnJob directly

`MultiwfnJob` is the low-level class that actually calls the Multiwfn binary.
`MultiwfnAnalysis` creates jobs internally, but you can use `MultiwfnJob` on
its own if you need manual control over the command sequence or want to skip
the analysis/caching layer entirely.

```python
from pymultiwfn import MultiwfnJob, Multiwfn, Menu
from pymultiwfn.analysis.result import MultiwfnResult

job = MultiwfnJob(
    input_file="benzene.wfn",
    analysis=Menu.HIRSHFELD_CHARGE,
    multiwfn=Multiwfn(exe_path=Path("/opt/multiwfn/Multiwfn")),
    timeout=60,
)
job = job.run()

print(job.stdout[:200])     # raw Multiwfn output
print(job.success)          # True if no error indicators in stderr

# Parse manually:
result = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
result.parse(job.stdout)
print(result.result)        # [Charge(...), Charge(...), Dipole(...)]
```

### Search for a Menu member

```python
from pymultiwfn import Menu

matches = Menu.search("HIRSHFELD")
# [Menu.HIRSHFELD_CHARGE, Menu.HIRSHFELD_I_CHARGE, Menu.HIRSHFELD_SURFACE, ...]
```