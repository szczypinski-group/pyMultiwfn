# Currently WIP
# pyMultiwfn

A Python wrapper for automating [Multiwfn](http://sobereva.com/multiwfn/) batch calculations.

## Installation
Currently not available on pypi yet
```bash
pip install pyMultiwfn
```

**Requirements:**
- Python 3.8+
- Multiwfn executable (place in `bin/` directory or specify path)

## Quick Start

```python
from pyMultiwfn import MultiwfnJob, Menu

job = MultiwfnJob("molecule.wfn")
job.add_menu(Menu.HIRSHFELD_CHARGE)
job.add_menu(Menu.MAYER_BOND_ORDER)

result = job.run()
charges = result.parse_charges()
bonds = result.parse_bond_orders()
```

## Usage

### Basic Job

```python
from pyMultiwfn import MultiwfnJob, Menu

job = MultiwfnJob("molecule.wfn")
job.add_menu(Menu.HIRSHFELD_CHARGE)
result = job.run()


### Builder Pattern

```python
from pyMultiwfn import MultiwfnJob, MultiwfnConfig, TimeoutConfig, Menu
from pathlib import Path

config = MultiwfnConfig(
    working_dir=Path("./output"),
    timeout=TimeoutConfig(default=120, topology=600),
    verbose=True
)

job = (
    MultiwfnJob.builder("molecule.wfn")
    .with_config(config)
    .with_menu(Menu.HIRSHFELD_CHARGE)
    .with_menu(Menu.TOPOLOGY_SEARCH_CPS)
    .with_menu(Menu.NCI_ANALYSIS)
    .build()
)

result = job.run()
```

### Available Analyses

```python
from pyMultiwfn import Menu

# List all available menu items
print(Menu.list_all())

# Search for specific analyses
charge_menus = Menu.search("charge")
bond_menus = Menu.search("bond")
```

### Common Menu Items

**Charges:**
- `Menu.HIRSHFELD_CHARGE`
- `Menu.MULLIKEN_POPULATION`
- `Menu.ADCH_CHARGE`
- `Menu.RESP_CHARGE`
- `Menu.CHELPG_CHARGE`

**Bond Orders:**
- `Menu.MAYER_BOND_ORDER`
- `Menu.WIBERG_BOND_ORDER`
- `Menu.FUZZY_BOND_ORDER`

**Topology:**
- `Menu.TOPOLOGY_SEARCH_CPS`
- `Menu.TOPOLOGY_ANALYSIS_COMPLETE`

**Weak Interactions:**
- `Menu.NCI_ANALYSIS`
- `Menu.IGM_ANALYSIS`
- `Menu.IRI_ANALYSIS`

**Cube Files:**
- `Menu.CUBE_DENSITY`
- `Menu.CUBE_ELF`
- `Menu.CUBE_ESP`

### Parsing Results

```python
result = job.run()

# Atomic charges
charges = result.parse_charges()

# Bond orders
bonds = result.parse_bond_orders()

# Critical points (QTAIM)
cps = result.parse_critical_points()

# Spectrum data
spectrum = result.parse_spectrum()

# Raw output
print(result.stdout)

# Save output
result.save_output("analysis.log")
```

### Custom Commands

```python
job = MultiwfnJob("molecule.wfn")
job.add_custom_commands(["7", "1", "1"])  # Manual menu navigation
result = job.run()
```

### Timeout Configuration for Lengthy Calculations

```python
from pyMultiwfn import TimeoutConfig

# Analysis-specific timeouts
timeout = TimeoutConfig(
    default=120,      # Simple analyses
    topology=600,     # QTAIM/topology
    batch=1200,       # Multiple analyses
    cube=300          # Cube generation
)
```

## Supported File Formats

- `.wfn` / `.wfx` - Wavefunction files
- `.fchk` - Gaussian formatted checkpoint
- `.molden` - Molden format
- `.cub` - cube format (Future)

## License

MIT License