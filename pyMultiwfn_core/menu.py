# pyMultiwfn/core.py
"""
Core functionality for pyMultiwfn - A Python wrapper for Multiwfn batch automation.
"""
import subprocess
import tempfile
import os
import re
from pathlib import Path
from typing import List, Dict, Union, Callable, Optional, Tuple, Any


__version__ = "0.1.0"


class MultiwfnError(Exception):
    """Base exception for Multiwfn-related errors."""
    pass


class MultiwfnExecutionError(MultiwfnError):
    """Raised when Multiwfn execution fails."""
    pass


class MultiwfnOutputParser:
    """Parse and extract useful information from Multiwfn output."""
    
    @staticmethod
    def extract_charges(stdout: str, method: str = "Hirshfeld") -> Dict[int, float]:
        """
        Extract atomic charges from Multiwfn output.
        
        Parameters
        ----------
        stdout : str
            Multiwfn standard output
        method : str
            Charge method name (e.g., "Hirshfeld", "ADCH", "RESP")
            
        Returns
        -------
        dict
            Dictionary mapping atom indices to charges
        """
        charges = {}
        pattern = r"^\s*(\d+)\(.*?\)\s+([-+]?\d+\.\d+)"
        
        in_charge_section = False
        for line in stdout.split('\n'):
            if method in line or "Atomic charges" in line:
                in_charge_section = True
                continue
            
            if in_charge_section:
                match = re.match(pattern, line)
                if match:
                    atom_idx = int(match.group(1))
                    charge = float(match.group(2))
                    charges[atom_idx] = charge
                elif line.strip() == "":
                    break
                    
        return charges
    
    @staticmethod
    def extract_bond_orders(stdout: str) -> Dict[Tuple[int, int], float]:
        """
        Extract bond orders from Multiwfn output.
        
        Parameters
        ----------
        stdout : str
            Multiwfn standard output
            
        Returns
        -------
        dict
            Dictionary mapping atom pairs to bond orders
        """
        bond_orders = {}
        pattern = r"^\s*(\d+)\s*-\s*(\d+)\s+([-+]?\d+\.\d+)"
        
        for line in stdout.split('\n'):
            match = re.match(pattern, line)
            if match:
                atom1 = int(match.group(1))
                atom2 = int(match.group(2))
                bo = float(match.group(3))
                bond_orders[(atom1, atom2)] = bo
                
        return bond_orders
    
    @staticmethod
    def extract_critical_points(stdout: str) -> List[Dict[str, Any]]:
        """
        Extract critical point information from topology analysis.
        
        Parameters
        ----------
        stdout : str
            Multiwfn standard output
            
        Returns
        -------
        list
            List of dictionaries containing CP information
        """
        cps = []
        pattern = r"CP\s+(\d+).*?Type:\s*\((\d+),([+-]?\d+)\).*?Position.*?:\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)"
        
        for match in re.finditer(pattern, stdout, re.DOTALL):
            cp = {
                'index': int(match.group(1)),
                'type': (int(match.group(2)), int(match.group(3))),
                'position': (float(match.group(4)), float(match.group(5)), float(match.group(6)))
            }
            cps.append(cp)
            
        return cps
    
    @staticmethod
    def extract_spectrum_data(stdout: str) -> Dict[str, List[float]]:
        """
        Extract spectrum data (frequencies, intensities).
        
        Parameters
        ----------
        stdout : str
            Multiwfn standard output
            
        Returns
        -------
        dict
            Dictionary with 'frequencies' and 'intensities' lists
        """
        spectrum = {'frequencies': [], 'intensities': []}
        
        # Pattern for vibrational data
        pattern = r"(\d+\.\d+)\s+cm\^-1.*?Intensity:\s+([-+]?\d+\.\d+)"
        
        for match in re.finditer(pattern, stdout):
            spectrum['frequencies'].append(float(match.group(1)))
            spectrum['intensities'].append(float(match.group(2)))
            
        return spectrum


class MultiwfnJob:
    """
    Encapsulates a Multiwfn job with input file, menu sequence, and execution.
    """
    
    def __init__(self, input_file: Union[str, Path], 
                 multiwfn_exe: Optional[Union[str, Path]] = None,
                 working_dir: Optional[Union[str, Path]] = None):
        """
        Initialize a Multiwfn job.
        
        Parameters
        ----------
        input_file : str or Path
            Path to wavefunction file
        multiwfn_exe : str or Path, optional
            Path to Multiwfn executable
        working_dir : str or Path, optional
            Working directory for output files
        """
        self.input_file = Path(input_file)
        if not self.input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
            
        self.multiwfn_exe = self._find_executable(multiwfn_exe)
        self.working_dir = Path(working_dir) if working_dir else Path.cwd()
        self.batch_commands = []
        self.stdout = None
        self.stderr = None
        self.returncode = None
        
    def _find_executable(self, exe_path: Optional[Union[str, Path]]) -> Path:
        """Find Multiwfn executable in common locations."""
        if exe_path:
            exe = Path(exe_path)
            if exe.exists():
                return exe
        
        # Check common locations
        search_paths = [
            Path(__file__).parent / "bin" / "Multiwfn.exe",
            Path(__file__).parent / "bin" / "Multiwfn",
            Path.home() / "Multiwfn" / "Multiwfn.exe",
            Path("/usr/local/bin/Multiwfn"),
        ]
        
        # Check PATH
        import shutil
        if shutil.which("Multiwfn"):
            return Path(shutil.which("Multiwfn"))
            
        for path in search_paths:
            if path.exists():
                return path
                
        raise MultiwfnError(
            "Multiwfn executable not found. Please specify path with multiwfn_exe parameter."
        )
    
    def add_menu_sequence(self, menu_func: Callable, **kwargs) -> 'MultiwfnJob':
        """
        Add a menu sequence from a menu function.
        
        Parameters
        ----------
        menu_func : callable
            Menu function from menus.py
        **kwargs
            Arguments to pass to the menu function
            
        Returns
        -------
        self
            For method chaining
        """
        sequence = menu_func(**kwargs) if kwargs else menu_func()
        self.batch_commands.extend(sequence)
        return self
    
    def add_custom_commands(self, commands: List[str]) -> 'MultiwfnJob':
        """
        Add custom command sequence.
        
        Parameters
        ----------
        commands : list of str
            List of commands to add
            
        Returns
        -------
        self
            For method chaining
        """
        self.batch_commands.extend(commands)
        return self
    
    def run(self, verbose: bool = False, timeout: Optional[int] = None) -> 'MultiwfnJob':
        """
        Execute the Multiwfn job.
        
        Parameters
        ----------
        verbose : bool
            Print stdout during execution
        timeout : int, optional
            Timeout in seconds
            
        Returns
        -------
        self
            For method chaining
        """
        # Ensure exit command
        if not self.batch_commands or self.batch_commands[-1] != "q":
            self.batch_commands.append("q")
        
        # Create temporary batch file
        with tempfile.NamedTemporaryFile(
            mode='w', 
            delete=False, 
            suffix=".inp",
            dir=self.working_dir,
            encoding='utf-8'
        ) as batch_file:
            batch_file.write("\n".join(self.batch_commands) + "\n")
            batch_path = batch_file.name
        
        try:
            # Change to working directory
            original_dir = Path.cwd()
            os.chdir(self.working_dir)
            
            # Run Multiwfn
            with open(batch_path, 'r', encoding='utf-8') as batch:
                proc = subprocess.Popen(
                    [str(self.multiwfn_exe), str(self.input_file)],
                    stdin=batch,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
                
                try:
                    self.stdout, self.stderr = proc.communicate(timeout=timeout)
                    self.returncode = proc.returncode
                except subprocess.TimeoutExpired:
                    proc.kill()
                    self.stdout, self.stderr = proc.communicate()
                    raise MultiwfnExecutionError(f"Multiwfn execution timed out after {timeout}s")
            
            if verbose:
                print(self.stdout)
                
            if self.returncode != 0 and self.returncode is not None:
                raise MultiwfnExecutionError(
                    f"Multiwfn failed with return code {self.returncode}\n"
                    f"STDERR: {self.stderr}"
                )
                
        finally:
            os.chdir(original_dir)
            # Clean up batch file
            try:
                os.unlink(batch_path)
            except:
                pass
                
        return self
    
    def get_output(self) -> str:
        """Get standard output from execution."""
        if self.stdout is None:
            raise MultiwfnError("Job has not been executed yet. Call run() first.")
        return self.stdout
    
    def get_errors(self) -> str:
        """Get standard error from execution."""
        if self.stderr is None:
            raise MultiwfnError("Job has not been executed yet. Call run() first.")
        return self.stderr
    
    def parse_charges(self, method: str = "Hirshfeld") -> Dict[int, float]:
        """Parse atomic charges from output."""
        return MultiwfnOutputParser.extract_charges(self.get_output(), method)
    
    def parse_bond_orders(self) -> Dict[Tuple[int, int], float]:
        """Parse bond orders from output."""
        return MultiwfnOutputParser.extract_bond_orders(self.get_output())
    
    def parse_critical_points(self) -> List[Dict[str, Any]]:
        """Parse critical points from output."""
        return MultiwfnOutputParser.extract_critical_points(self.get_output())
    
    def save_output(self, filename: Union[str, Path]) -> None:
        """Save output to file."""
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(self.get_output())


# Legacy functions for backward compatibility
def generate_batch_file(menu_functions: List[Callable], 
                       extra_inputs: Optional[Dict[str, List[str]]] = None,
                       batch_file: Optional[Union[str, Path]] = None) -> str:
    """
    Generate a batch input file for Multiwfn from menu sequences.
    
    **DEPRECATED**: Use MultiwfnJob class for new code.
    
    Parameters
    ----------
    menu_functions : list
        List of menu functions from menus.py
    extra_inputs : dict, optional
        Extra inputs for specific menu functions
    batch_file : str or Path, optional
        Filename to write; if None, a temporary file is created
        
    Returns
    -------
    str
        Path to the batch file
    """
    batch_commands = []
    
    for func in menu_functions:
        seq = func()
        batch_commands.extend(seq)
        
        # Append extra inputs if provided
        if extra_inputs and func.__name__ in extra_inputs:
            batch_commands.extend(extra_inputs[func.__name__])
    
    # Ensure Multiwfn exits at the end
    if not batch_commands or batch_commands[-1] != "q":
        batch_commands.append("q")
    
    # Use temporary file if none provided
    if batch_file is None:
        tmp = tempfile.NamedTemporaryFile(
            mode='w+', 
            delete=False, 
            suffix=".inp",
            encoding='utf-8'
        )
        batch_file = tmp.name
        tmp.write("\n".join(batch_commands) + "\n")
        tmp.close()
    else:
        with open(batch_file, "w", encoding="utf-8") as f:
            f.write("\n".join(batch_commands) + "\n")
    
    return batch_file


def run_multiwfn(input_file: Union[str, Path],
                batch_file: Union[str, Path],
                multiwfn_exe: Optional[Union[str, Path]] = None) -> Tuple[str, str]:
    """
    Run Multiwfn in batch mode with the given input and batch file.
    
    **DEPRECATED**: Use MultiwfnJob class for new code.
    
    Parameters
    ----------
    input_file : str or Path
        Path to the wavefunction file
    batch_file : str or Path
        Path to the batch file
    multiwfn_exe : str or Path, optional
        Path to Multiwfn executable
        
    Returns
    -------
    tuple
        (stdout, stderr) from Multiwfn execution
    """
    job = MultiwfnJob(input_file, multiwfn_exe)
    
    # Read batch file and add commands
    with open(batch_file, 'r', encoding='utf-8') as f:
        commands = [line.strip() for line in f if line.strip()]
    
    job.add_custom_commands(commands)
    job.run()
    
    return job.get_output(), job.get_errors()


# Convenience function for common workflows
def quick_analysis(input_file: Union[str, Path],
                  analysis_type: str,
                  **kwargs) -> Dict[str, Any]:
    """
    Perform a quick analysis with common presets.
    
    Parameters
    ----------
    input_file : str or Path
        Path to wavefunction file
    analysis_type : str
        Type of analysis: 'charges', 'bonds', 'topology', 'nci', etc.
    **kwargs
        Additional arguments passed to specific analysis
        
    Returns
    -------
    dict
        Analysis results
    """
    from . import menus  # Import menu functions
    
    job = MultiwfnJob(input_file, **kwargs)
    results = {}
    
    if analysis_type == 'charges':
        method = kwargs.get('method', 'hirshfeld')
        if method.lower() == 'hirshfeld':
            job.add_menu_sequence(menus.hirshfeld_charge)
        elif method.lower() == 'adch':
            job.add_menu_sequence(menus.adch_charge)
        elif method.lower() == 'resp':
            job.add_menu_sequence(menus.resp_charge)
        elif method.lower() == 'cm5':
            job.add_menu_sequence(menus.cm5_charge)
        else:
            raise ValueError(f"Unknown charge method: {method}")
        
        job.run()
        results['charges'] = job.parse_charges(method)
        
    elif analysis_type == 'bonds':
        method = kwargs.get('method', 'mayer')
        if method.lower() == 'mayer':
            job.add_menu_sequence(menus.mayer_bond_order)
        elif method.lower() == 'wiberg':
            job.add_menu_sequence(menus.wiberg_bond_order)
        elif method.lower() == 'fuzzy':
            job.add_menu_sequence(menus.fuzzy_bond_order)
        else:
            raise ValueError(f"Unknown bond order method: {method}")
        
        job.run()
        results['bond_orders'] = job.parse_bond_orders()
        
    elif analysis_type == 'topology':
        job.add_menu_sequence(menus.topology_analysis_complete)
        job.run()
        results['critical_points'] = job.parse_critical_points()
        
    elif analysis_type == 'nci':
        job.add_menu_sequence(menus.nci_analysis)
        job.run()
        results['output'] = job.get_output()
        
    else:
        raise ValueError(f"Unknown analysis type: {analysis_type}")
    
    results['stdout'] = job.get_output()
    return results


# Export public API
__all__ = [
    'MultiwfnJob',
    'MultiwfnError',
    'MultiwfnExecutionError',
    'MultiwfnOutputParser',
    'generate_batch_file',
    'run_multiwfn',
    'quick_analysis',
]