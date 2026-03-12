# Getting started

Here we have some simple example use cases.

## Run a single analysis

```python
from pathlib import Path
from pymultiwfn import MultiwfnAnalysis, Multiwfn, Menu

analysis = MultiwfnAnalysis("benzene.wfn", analyses=Menu.HIRSHFELD_CHARGE)
analysis.run(multiwfn=mwfn, work_dir=Path("./output"))

result = analysis.results[0]          # MultiwfnResult for HIRSHFELD_CHARGE
for item in result.result:
    print(item)                        # Charge(...) and Dipole(...) objects
```

## Queue several analyses

```python
analysis = MultiwfnAnalysis("benzene.wfn")
analysis.add_menu(Menu.HIRSHFELD_CHARGE)
analysis.add_menu(Menu.MAYER_BOND_ORDER)
analysis.add_menu(Menu.TOPOLOGY_SEARCH_CPS)
analysis.run(multiwfn=mwfn, timeout=120, work_dir=Path("./output"))

# results[0] → charges,  results[1] → bond orders,  results[2] → critical points
```

## Queue an entire category

```python
from pymultiwfn.enums.analyses import AnalysisClasses

analysis = MultiwfnAnalysis("benzene.wfn")
analysis.add_menu(AnalysisClasses.CHARGES)   # queues all 19 charge methods
analysis.run(multiwfn=mwfn, work_dir=Path("./output"))
```

## Filter results by type

```python
from pymultiwfn.analysis.result import Charge, Dipole

for r in analysis.results:
    charges = [x for x in r.result if isinstance(x, Charge)]
    dipoles = [x for x in r.result if isinstance(x, Dipole)]
```

## Re-run with caching

```python
# Second run loads from benzene.wfn.json instead of re-executing Multiwfn.
analysis2 = MultiwfnAnalysis("benzene.wfn", analyses=Menu.HIRSHFELD_CHARGE)
analysis2.run(multiwfn=mwfn, work_dir=Path("./output"))
# No subprocess launched — result loaded from cache.
```

## Disable caching

```python
analysis = MultiwfnAnalysis("benzene.wfn", analyses=Menu.HIRSHFELD_CHARGE, cached=False)
analysis.run(multiwfn=mwfn, work_dir=Path("./output"))
# Always re-runs Multiwfn, even if the JSON file has a stored result.
```

## Use ResultStore directly

```python
from pymultiwfn.analysis.result import ResultStore, MultiwfnResult

store = ResultStore(input_file=Path("benzene.wfn"), work_dir=Path("./output"))

# Parse raw stdout from some other source and persist it:
result = store.store_from_stdout(Menu.MAYER_BOND_ORDER, raw_stdout_string)

# Check what is stored:
print(store.has_result(Menu.MAYER_BOND_ORDER))  # True
print(store.get_result(Menu.MAYER_BOND_ORDER))  # dict with "parsed" and "timestamp"
```

## Use MultiwfnJob directly

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

## Search for a Menu member

```python
from pymultiwfn import Menu

matches = Menu.search("HIRSHFELD")
# [Menu.HIRSHFELD_CHARGE, Menu.HIRSHFELD_I_CHARGE, Menu.HIRSHFELD_SURFACE, ...]
```