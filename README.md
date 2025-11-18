# pyMultiwfn

[![PyPI version](https://badge.fury.io/py/pyMultiwfn.svg)](https://badge.fury.io/py/pyMultiwfn)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyMultiwfn.svg)](https://pypi.org/project/pyMultiwfn/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive Python wrapper for [Multiwfn](http://sobereva.com/multiwfn/) wavefunction analysis automation. This library provides a clean, Pythonic interface to Multiwfn's extensive functionality, enabling high-throughput batch processing and seamless integration into computational chemistry workflows.

## Features

- 🚀 **200+ Menu Functions**: Complete coverage of Multiwfn's functionality
- 🔧 **Method Chaining**: Fluent API for building complex workflows
- 📊 **Output Parsing**: Automatic extraction of charges, bond orders, critical points, etc.
- 🎯 **Type Hints**: Full typing support for better IDE integration
- 📦 **PyPI Ready**: Easy installation and distribution
- 🔄 **Batch Processing**: Efficient handling of multiple molecules
- 💾 **Flexible I/O**: Support for all Multiwfn file formats

## Installation

```bash
pip install pyMultiwfn
```

### Requirements

- Python ≥ 3.8
- Multiwfn (must be installed separately from [official website](http://sobereva.com/multiwfn/))

## Quick Start

### Basic Usage

```python
from pyMultiwfn import MultiwfnJob
from pyMultiwfn import menus

# Create a job with your wavefunction file
job = MultiwfnJob("molecule.molden")

# Add analyses using method chaining
job.add_menu_sequence(menus.hirshfeld_charge)
job.add_menu_sequence(menus.mayer_bond_order)

# Execute and parse results
job.run()
charges = job.parse_charges("Hirshfeld")
bonds = job.parse_bond_orders()

print(f"Atom 1 charge: {charges[1]:.4f}")
print(f"Bond order (1-2): {bonds[(1, 2)]:.4f}")
```

### Quick Analysis

For common analyses, use the convenience function:

```python
from pyMultiwfn import quick_analysis

# ADCH charges
results = quick_analysis("molecule.molden", "charges", method="adch")
print(results['charges'])

# Wiberg bond orders
results = quick_analysis("molecule.molden", "bonds", method="wiberg")
print(results['bond_orders'])

# Topology analysis
results = quick_analysis("molecule.molden", "topology")
print(results['critical_points'])
```

## Comprehensive Examples

### 1. Population Analysis Suite

```python
from pyMultiwfn import MultiwfnJob, menus

job = MultiwfnJob("molecule.wfn")

# Compare multiple charge methods
charge_methods = [
    (menus.hirshfeld_charge, "Hirshfeld"),
    (menus.adch_charge, "ADCH"),
    (menus.resp_charge, "RESP"),
    (menus.cm5_charge, "CM5"),
]

results = {}
for menu_func, method_name in charge_methods:
    job_temp = MultiwfnJob("molecule.wfn")
    job_temp.add_menu_sequence(menu_func).run()
    results[method_name] = job_temp.parse_charges(method_name)

# Compare results
for method, charges in results.items():
    print(f"\n{method} Charges:")
    for atom_idx, charge in charges.items():
        print(f"  Atom {atom_idx}: {charge:+.4f}")
```

### 2. Complete Topology Analysis

```python
from pyMultiwfn import MultiwfnJob, menus

job = MultiwfnJob("molecule.fch")

# Full QTAIM analysis
job.add_menu_sequence(menus.topology_analysis_complete)
job.run()

# Extract critical points
cps = job.parse_critical_points()

# Classify critical points
for cp in cps:
    cp_type = cp['type']
    position = cp['position']
    
    if cp_type == (3, -3):
        print(f"Nucleus at {position}")
    elif cp_type == (3, -1):
        print(f"Bond CP at {position}")
    elif cp_type == (3, +1):
        print(f"Ring CP at {position}")
    elif cp_type == (3, +3):
        print(f"Cage CP at {position}")
```

### 3. Weak Interaction Analysis

```python
from pyMultiwfn import MultiwfnJob, menus

job = MultiwfnJob("complex.molden")

# Multiple weak interaction methods
job.add_menu_sequence(menus.nci_analysis)
job.add_menu_sequence(menus.igm_analysis)
job.add_menu_sequence(menus.iri_analysis)

job.run(verbose=True)

# Output will generate visualization files
# Look for .cub, .pdb, and .png files in working directory
```

### 4. Spectroscopy Suite

```python
from pyMultiwfn import MultiwfnJob, menus

# IR spectrum
job_ir = MultiwfnJob("molecule.out")
job_ir.add_menu_sequence(menus.plot_ir_spectrum).run()

# UV-Vis spectrum
job_uv = MultiwfnJob("molecule.out")
job_uv.add_menu_sequence(menus.plot_uv_vis_spectrum).run()

# NMR spectrum
job_nmr = MultiwfnJob("molecule.out")
job_nmr.add_menu_sequence(menus.plot_nmr_spectrum).run()

# Color prediction from UV-Vis
job_color = MultiwfnJob("molecule.out")
job_color.add_menu_sequence(menus.predict_color).run()
print(job_color.get_output())
```

### 5. Aromaticity Analysis

```python
from pyMultiwfn import MultiwfnJob, menus

job = MultiwfnJob("aromatic.molden")

# Multiple aromaticity indices
job.add_menu_sequence(menus.homa_index)
job.add_menu_sequence(menus.nics_1d_scan)
job.add_menu_sequence(menus.flu_aromaticity)
job.add_menu_sequence(menus.pdi_aromaticity)

job.run()

# Parse output for aromaticity values
output = job.get_output()
# ... extract specific values based on your needs
```

### 6. Orbital Analysis

```python
from pyMultiwfn import MultiwfnJob, menus

job = MultiwfnJob("molecule.fch")

# Analyze HOMO composition
job.add_menu_sequence(menus.orbital_composition_nao, orbital_index="HOMO")
job.run()

# Localize orbitals
job2 = MultiwfnJob("molecule.fch")
job2.add_menu_sequence(menus.pipek_mezey_localization)
job2.run()

# Generate cube files for visualization
job3 = MultiwfnJob("molecule.fch")
job3.add_menu_sequence(menus.cube_orbital, orbital_index=25)
job3.run()
```

### 7. Batch Processing Multiple Molecules

```python
from pathlib import Path
from pyMultiwfn import MultiwfnJob, menus

# Process all .molden files in directory
molden_files = Path("./molecules").glob("*.molden")

results = {}
for molfile in molden_files:
    try:
        job = MultiwfnJob(molfile)
        job.add_menu_sequence(menus.hirshfeld_charge)
        job.add_menu_sequence(menus.mayer_bond_order)
        job.run(timeout=300)  # 5 minute timeout
        
        results[molfile.stem] = {
            'charges': job.parse_charges(),
            'bonds': job.parse_bond_orders()
        }
        
    except Exception as e:
        print(f"Failed to process {molfile}: {e}")
        continue

# Save results
import json
with open("batch_results.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
```

### 8. Custom Workflow with Error Handling

```python
from pyMultiwfn import MultiwfnJob, MultiwfnError, menus

def analyze_molecule(input_file, output_dir):
    """Complete analysis workflow with error handling."""
    try:
        job = MultiwfnJob(input_file, working_dir=output_dir)
        
        # Add multiple analyses
        analyses = [
            menus.adch_charge,
            menus.wiberg_bond_order,
            menus.fukui_function,
            menus.nci_analysis,
        ]
        
        for analysis in analyses:
            job.add_menu_sequence(analysis)
        
        # Run with timeout
        job.run(timeout=600, verbose=True)
        
        # Parse all results
        results = {
            'charges': job.parse_charges("ADCH"),
            'bonds': job.parse_bond_orders(),
            'stdout': job.get_output(),
        }
        
        # Save output
        job.save_output(output_dir / "analysis.log")
        
        return results
        
    except MultiwfnError as e:
        print(f"Multiwfn error: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        raise

# Usage
results = analyze_molecule("molecule.molden", Path("./output"))
```

## Available Menu Functions

### Population Analysis
- `hirshfeld_charge()` - Hirshfeld charges
- `adch_charge()` - ADCH charges
- `resp_charge()` - RESP charges
- `chelpg_charge()` - CHELPG charges
- `mulliken_population()` - Mulliken population
- `lowdin_population()` - Löwdin population
- `becke_charge()` - Becke charges
- `cm5_charge()` - CM5 charges
- ... and more

### Bond Analysis
- `mayer_bond_order()` - Mayer bond orders
- `wiberg_bond_order()` - Wiberg bond orders
- `fuzzy_bond_order()` - Fuzzy bond orders
- `laplacian_bond_order()` - Laplacian bond orders
- `ibsi_analysis()` - Intrinsic bond strength index
- `multicenter_bond_order()` - Multi-center bond orders

### Topology
- `topology_analysis_complete()` - Full QTAIM analysis
- `topology_search_cps()` - Find critical points
- `topology_generate_paths()` - Generate bond paths
- `topology_esp_analysis()` - ESP topology

### Weak Interactions
- `nci_analysis()` - NCI analysis
- `igm_analysis()` - IGM analysis
- `igmh_analysis()` - IGMH analysis
- `iri_analysis()` - IRI analysis
- `dori_analysis()` - DORI analysis
- `vdw_potential()` - van der Waals potential

### Spectroscopy
- `plot_ir_spectrum()` - IR spectrum
- `plot_uv_vis_spectrum()` - UV-Vis spectrum
- `plot_nmr_spectrum()` - NMR spectrum
- `plot_ecd_spectrum()` - ECD spectrum
- `plot_raman_spectrum()` - Raman spectrum

### Aromaticity
- `homa_index()` - HOMA index
- `nics_scan()` - NICS scan
- `flu_aromaticity()` - FLU index
- `pdi_aromaticity()` - PDI index
- `icss_analysis()` - ICSS analysis

[See full list of 200+ functions in documentation]

## API Reference

### MultiwfnJob Class

```python
class MultiwfnJob:
    """Main class for Multiwfn job management."""
    
    def __init__(self, input_file, multiwfn_exe=None, working_dir=None):
        """Initialize job with wavefunction file."""
        
    def add_menu_sequence(self, menu_func, **kwargs):
        """Add a menu sequence. Returns self for chaining."""
        
    def add_custom_commands(self, commands):
        """Add custom command list. Returns self for chaining."""
        
    def run(self, verbose=False, timeout=None):
        """Execute the job. Returns self for chaining."""
        
    def get_output(self):
        """Get stdout from execution."""
        
    def parse_charges(self, method="Hirshfeld"):
        """Parse atomic charges from output."""
        
    def parse_bond_orders(self):
        """Parse bond orders from output."""
        
    def parse_critical_points(self):
        """Parse critical points from topology analysis."""
        
    def save_output(self, filename):
        """Save output to file."""
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- [Multiwfn](http://sobereva.com/multiwfn/) by Tian Lu - The powerful wavefunction analysis program this wrapper interfaces with
- The computational chemistry community

## Citation

If you use pyMultiwfn in your research, please cite both this package and Multiwfn:

**Multiwfn:**
> Lu, T.; Chen, F. Multiwfn: A Multifunctional Wavefunction Analyzer. *J. Comput. Chem.* **2012**, 33, 580-592.

## Links

- **Documentation**: https://pymultiwfn.readthedocs.io
- **PyPI**: https://pypi.org/project/pyMultiwfn/
- **GitHub**: https://github.com/yourusername/pyMultiwfn
- **Multiwfn**: http://sobereva.com/multiwfn/