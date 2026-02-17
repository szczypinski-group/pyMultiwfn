"""
Pytest configuration and fixtures for pyMultiwfn tests.

Test data files should be placed in tests/test_data/:
- M062X_TZVPP_D30.wfx
- M062X_TZVPP_D30.molden
"""

import platform
import sys
from pathlib import Path

import pytest

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


# =============================================================================
# Path Fixtures
# =============================================================================


@pytest.fixture(scope="session")
def project_root() -> Path:
    """Path to project root (PYMULTIWFN/)."""
    return Path(__file__).parent.parent


@pytest.fixture(scope="session")
def test_data_dir(project_root: Path) -> Path:
    """Path to test data directory."""
    return project_root / "tests" / "test_data"


@pytest.fixture(scope="session")
def real_wfx_file(test_data_dir: Path) -> Path:
    """Real .wfx test file."""
    path = test_data_dir / "M062X_TZVPP_D30.wfx"
    if not path.exists():
        pytest.skip(f"Test file not found: {path}")
    return path


@pytest.fixture(scope="session")
def real_molden_file(test_data_dir: Path) -> Path:
    """Real .molden test file."""
    path = test_data_dir / "M062X_TZVPP_D30.molden"
    if not path.exists():
        pytest.skip(f"Test file not found: {path}")
    return path


@pytest.fixture(scope="session")
def real_executable(project_root: Path) -> Path:
    """Real Multiwfn executable from bin/ folder."""
    exe_name = "Multiwfn.exe" if platform.system() == "Windows" else "Multiwfn"
    exe_path = project_root / "bin" / exe_name

    if not exe_path.exists():
        pytest.skip(f"Multiwfn executable not found: {exe_path}")
    return exe_path


# =============================================================================
# Temporary Fixtures
# =============================================================================


@pytest.fixture
def temp_dir(tmp_path: Path) -> Path:
    """Temporary directory for test outputs."""
    return tmp_path


@pytest.fixture
def mock_wfn_file(temp_dir: Path) -> Path:
    """Create a mock wavefunction file."""
    path = temp_dir / "mock.wfn"
    path.write_text("mock wavefunction content")
    return path


@pytest.fixture
def mock_executable(temp_dir: Path) -> Path:
    """Create a mock Multiwfn executable."""
    if platform.system() == "Windows":
        exe_path = temp_dir / "Multiwfn.exe"
        exe_path.write_text("@echo mock")
    else:
        exe_path = temp_dir / "Multiwfn"
        exe_path.write_text("#!/bin/bash\necho 'mock'")
        exe_path.chmod(0o755)
    return exe_path


# =============================================================================
# Sample Output Fixtures
# =============================================================================


@pytest.fixture
def sample_hirshfeld_output() -> str:
    """Sample Multiwfn output for Hirshfeld charge analysis."""
    return """
 Hirshfeld charges:
     1(C )   -0.0523
     2(C )   -0.0892
     3(H )    0.0534
     4(H )    0.0521
     5(O )   -0.3245
     6(H )    0.1823
 Sum of Hirshfeld charges:  -0.1782
"""


@pytest.fixture
def sample_mulliken_output() -> str:
    """Sample Multiwfn output for Mulliken population."""
    return """
 Mulliken atomic charges:
     1(C )   -0.1234
     2(C )   -0.2345
     3(H )    0.1122
     4(H )    0.1133
     5(O )   -0.4567
     6(H )    0.2891
 Sum of Mulliken charges:  -0.3000
"""


@pytest.fixture
def sample_mayer_output() -> str:
    """Sample Multiwfn output for Mayer bond order analysis."""
    return """
 Mayer bond orders:
    1(C ) -    2(C ):  1.4523
    1(C ) -    3(H ):  0.9234
    1(C ) -    4(H ):  0.9123
    2(C ) -    5(O ):  1.8934
    5(O ) -    6(H ):  0.8234
 
 Total Mayer bond order:  6.0048
"""


@pytest.fixture
def sample_wiberg_output() -> str:
    """Sample Multiwfn output for Wiberg bond order analysis."""
    return """
 Wiberg bond orders (NAO basis):
    1(C ) -    2(C ):  1.5234
    1(C ) -    3(H ):  0.9567
    2(C ) -    5(O ):  1.9123
"""


@pytest.fixture
def sample_topology_output() -> str:
    """Sample Multiwfn output for topology analysis."""
    return """
 Searching critical points...
 
 CP   1 (3,-3) Nuclear critical point
   Position (Bohr):    0.000000    0.000000    0.000000
   Density:  2.834523E+02
   
 CP   2 (3,-1) Bond critical point
   Position (Bohr):    1.234567    0.000000    0.000000
   Density:  2.534523E-01
   
 CP   3 (3,+1) Ring critical point
   Position (Bohr):    0.500000    0.866025    0.000000
   Density:  1.234523E-02
   
 CP   4 (3,+3) Cage critical point
   Position (Bohr):    0.333333    0.577350    0.816497
   Density:  5.234523E-03

 Total CPs found: 4
"""


@pytest.fixture
def sample_spectrum_output() -> str:
    """Sample Multiwfn output for spectrum analysis."""
    return """
 IR Spectrum:
 Frequency(cm^-1)    Intensity
      500.00           12.34
     1000.00           45.67
     1500.00           23.45
     2000.00           78.90
     3000.00          156.78
"""


# =============================================================================
# Pytest Markers
# =============================================================================


def pytest_configure(config: pytest.Config) -> None:
    """Register custom markers."""
    config.addinivalue_line("markers", "slow: marks tests as slow")
    config.addinivalue_line("markers", "integration: marks integration tests")
    config.addinivalue_line(
        "markers",
        "requires_multiwfn: marks tests requiring Multiwfn executable",
    )
    config.addinivalue_line(
        "markers", "real_files: marks tests using real wavefunction files"
    )
