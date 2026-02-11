pyMultiwfn

Python wrapper for Multiwfn wavefunction analysis automation
Features

    🚀 190+ Menu Functions: Complete coverage of Multiwfn's analysis capabilities.

    🔧 Fluent API: Method chaining for building complex workflows.

    📊 Output Parsing: Automatic extraction of charges, bond orders, and critical points.

    🎯 Type Hints: Full typing support for IDE integration.

Installation
Bash

pip install pyMultiwfn

    Note: Requires Multiwfn to be installed separately.

Quick Example
Python

from pyMultiwfn import MultiwfnJob, menu

# Create job and add analyses
job = MultiwfnJob("molecule.wfx")
job.add_menu_sequence(menu.hirshfeld_charge)
job.add_menu_sequence(menu.mayer_bond_order)

# Run and parse results
result = job.run()
charges = result.parse_charges('hirshfeld')
bonds = result.parse_bond_orders()

print(f"Atom 1 charge: {charges[1]:.4f}")

Module Overview
Module	Description	Key Classes/Functions
config.py	Configuration and executable discovery	MultiwfnConfig, MultiwfnError
job.py	Job creation and execution	MultiwfnJob, MultiwfnJobBuilder
result.py	Execution results container	MultiwfnResult
parsers.py	Output parsing utilities	ChargeParser, BondOrderParser
analysis.py	Convenience analysis functions	run_analysis, quick_analysis
menu.py	190+ Multiwfn menu sequences	hirshfeld_charge, nci_analysis
Links

    GitHub Repository

    Official Multiwfn Website
