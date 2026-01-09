"""
Analysis utilities for pyMultiwfn.

Provides simple wrappers to run ANY menu.py function.
"""

from pathlib import Path
from typing import Any, Callable, Dict, List, Union

from .config import MultiwfnConfig
from .job import MultiwfnJob
from .result import MultiwfnResult


def _get_all_menu_functions() -> Dict[str, Callable]:
    """Get all callable functions from menu module."""
    from . import menu
    return {
        name: getattr(menu, name)
        for name in dir(menu)
        if callable(getattr(menu, name)) and not name.startswith('_')
    }


def list_functions() -> List[str]:
    """List all available menu functions."""
    return sorted(_get_all_menu_functions().keys())


def get_function(name: str) -> Callable:
    """Get a menu function by name."""
    funcs = _get_all_menu_functions()
    if name not in funcs:
        raise ValueError(f"Unknown function: {name}. Use list_functions() to see available.")
    return funcs[name]


def _make_config(**kwargs) -> MultiwfnConfig:
    """Create config from kwargs, removing config-related keys."""
    return MultiwfnConfig(
        exe_path=Path(kwargs.pop('exe_path')) if kwargs.get('exe_path') else None,
        working_dir=Path(kwargs.pop('working_dir')) if kwargs.get('working_dir') else Path.cwd(),
        timeout=kwargs.pop('timeout', None),
        verbose=kwargs.pop('verbose', False)
    )


def run_analysis(
    input_file: Union[str, Path],
    analysis_name: str,
    **kwargs
) -> MultiwfnResult:
    """
    Run any analysis by name.
    
    Parameters
    ----------
    input_file : str or Path
        Path to wavefunction file
    analysis_name : str
        Name of the menu function (e.g., 'hirshfeld_charge', 'nci_analysis')
    **kwargs
        - func_kwargs: dict of kwargs to pass to the menu function
        - exe_path, timeout, working_dir, verbose for config
        
    Returns
    -------
    MultiwfnResult
        Result object
        
    Examples
    --------
    >>> result = run_analysis("mol.wfx", "hirshfeld_charge")
    >>> result = run_analysis("mol.wfx", "nci_analysis")
    >>> result = run_analysis("mol.wfx", "cube_orbital", func_kwargs={'orbital_index': 5})
    >>> result = run_analysis("mol.wfx", "hole_electron_analysis", func_kwargs={'state_num': '2'})
    """
    func_kwargs = kwargs.pop('func_kwargs', {})
    config = _make_config(**kwargs)
    menu_func = get_function(analysis_name)
    
    job = MultiwfnJob(input_file, config=config)
    job.add_menu_sequence(menu_func, **func_kwargs)
    return job.run()


def run_menu_sequence(
    input_file: Union[str, Path],
    menu_func: Callable,
    **kwargs
) -> MultiwfnResult:
    """
    Run a menu function directly.
    
    Parameters
    ----------
    input_file : str or Path
        Path to wavefunction file
    menu_func : callable
        Menu function from pyMultiwfn.menu
    **kwargs
        - func_kwargs: dict of kwargs to pass to menu_func
        - exe_path, timeout, working_dir, verbose for config
        
    Returns
    -------
    MultiwfnResult
        Result object
    """
    func_kwargs = kwargs.pop('func_kwargs', {})
    config = _make_config(**kwargs)
    
    job = MultiwfnJob(input_file, config=config)
    job.add_menu_sequence(menu_func, **func_kwargs)
    return job.run()


def run_custom_sequence(
    input_file: Union[str, Path],
    commands: List[str],
    **kwargs
) -> MultiwfnResult:
    """
    Run a custom command sequence.
    
    Parameters
    ----------
    input_file : str or Path
        Path to wavefunction file
    commands : list of str
        List of menu commands to execute
    **kwargs
        exe_path, timeout, working_dir, verbose
        
    Returns
    -------
    MultiwfnResult
        Result object
    """
    config = _make_config(**kwargs)
    job = MultiwfnJob(input_file, config=config)
    job.add_custom_commands(commands)
    return job.run()


# =============================================================================
# Convenience functions for common analyses with parsing
# =============================================================================

def quick_analysis(
    input_file: Union[str, Path],
    analysis_type: str,
    **kwargs
) -> Dict[str, Any]:
    """
    Perform a quick analysis with automatic parsing.
    
    Parameters
    ----------
    input_file : str or Path
        Path to wavefunction file
    analysis_type : str or AnalysisType
        Type: 'charges', 'bonds', 'topology', 'nci'
    **kwargs
        - method: Analysis method for charges/bonds
        - exe_path, timeout, working_dir, verbose
        
    Returns
    -------
    dict
        Parsed results plus raw stdout
    """
    method = kwargs.pop('method', None)
    config = _make_config(**kwargs)
    job = MultiwfnJob(input_file, config=config)
    results = {}
    
    analysis_type_str = str(analysis_type).lower()
    
    if analysis_type_str == 'charges':
        method = method or 'hirshfeld'
        method_str = str(method).lower()
        # Map method names to function names
        func_name = {
            'hirshfeld': 'hirshfeld_charge',
            'mulliken': 'mulliken_population',
            'lowdin': 'lowdin_population',
            'vdd': 'vdd_population',
        }.get(method_str, f"{method_str}_charge")
        
        job.add_menu_sequence(get_function(func_name))
        result = job.run()
        results['charges'] = result.parse_charges(method)
        
    elif analysis_type_str == 'bonds':
        method = method or 'mayer'
        func_name = f"{str(method).lower()}_bond_order"
        job.add_menu_sequence(get_function(func_name))
        result = job.run()
        results['bond_orders'] = result.parse_bond_orders()
        
    elif analysis_type_str == 'topology':
        job.add_menu_sequence(get_function('topology_analysis_complete'))
        result = job.run()
        results['critical_points'] = result.parse_critical_points()
        
    elif analysis_type_str == 'nci':
        job.add_menu_sequence(get_function('nci_analysis'))
        result = job.run()
        results['output'] = result.stdout
        
    else:
        raise ValueError(f"Unknown analysis type: {analysis_type}. Use: charges, bonds, topology, nci")
    
    results['stdout'] = job.result.stdout
    results['success'] = job.result.success
    results['execution_time'] = job.result.execution_time
    return results


def charge_analysis(
    input_file: Union[str, Path],
    method: str = 'hirshfeld',
    **kwargs
) -> Dict[int, float]:
    """Calculate atomic charges."""
    return quick_analysis(input_file, 'charges', method=method, **kwargs)['charges']


def bond_order_analysis(
    input_file: Union[str, Path],
    method: str = 'mayer',
    **kwargs
) -> Dict[tuple, float]:
    """Calculate bond orders."""
    return quick_analysis(input_file, 'bonds', method=method, **kwargs)['bond_orders']


def topology_analysis(input_file: Union[str, Path], **kwargs) -> List[Dict[str, Any]]:
    """Perform QTAIM topology analysis."""
    return quick_analysis(input_file, 'topology', **kwargs)['critical_points']


def nci_analysis(input_file: Union[str, Path], **kwargs) -> str:
    """Perform NCI analysis."""
    return quick_analysis(input_file, 'nci', **kwargs)['output']