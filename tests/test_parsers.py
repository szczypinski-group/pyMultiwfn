"""Tests for parsers.py - Output parsers.

Every test class includes:
  - Positive assertions (happy path, correct values extracted)
  - Negative assertions (wrong keys absent, irrelevant data ignored,
    non-matching patterns rejected)
  - Edge cases (empty input, scientific notation, negative values,
    zero values, single items, large indices, whitespace variation,
    mixed/garbage output, section boundaries, partial data)
"""

import pytest

from pymultiwfn.parsers import (
    AromaticityParser,
    BasinParser,
    BondOrderParser,
    CDFTParser,
    ChargeParser,
    CriticalPointParser,
    CubeParser,
    DOSParser,
    EDAParser,
    ExcitationParser,
    FuzzySpaceParser,
    OrbitalCompositionParser,
    PolarizabilityParser,
    SpectrumParser,
    SurfaceParser,
    UtilityParser,
    WavefunctionParser,
    WeakInteractionParser,
)


# =============================================================================
# Fixtures — shared across multiple test classes
# =============================================================================


@pytest.fixture()
def sample_mayer_output() -> str:
    return (
        "Mayer bond order analysis\n"
        "Bond order list:\n"
        "#    1:         1(C )    2(C )    1.45230000\n"
        "#    2:         1(C )    3(H )    0.92340000\n"
        "#    3:         2(C )    4(H )    0.95670000\n"
        "#    4:         2(C )    5(O )    1.89340000\n"
        "#    5:         5(O )    6(N )    0.82340000\n"
    )


@pytest.fixture()
def sample_wiberg_output() -> str:
    return (
        "Wiberg bond order analysis\n"
        "#    1:         1(C )    2(C )    1.52340000\n"
        "#    2:         1(C )    3(H )    0.89560000\n"
        "#    3:         2(C )    4(H )    0.91230000\n"
    )


@pytest.fixture()
def sample_topology_output() -> str:
    return (
        "Topology analysis\n"
        "CP  1 (3,-3)\n"
        "Position (Bohr):   0.000000   0.000000   0.000000\n"
        "Density of all electrons:   0.29876500\n"
        "Laplacian of electron density:  -1.12345600\n"
        "Ellipticity:   0.00000000\n"
        "\n"
        "CP  2 (3,-1)\n"
        "Position (Bohr):   1.234567   0.567890   0.000000\n"
        "Density of all electrons:   0.25432100\n"
        "Laplacian of electron density:  -0.87654300\n"
        "Ellipticity:   0.04523000\n"
        "\n"
        "CP  3 (3,+1)\n"
        "Position (Bohr):   0.500000   0.500000   0.500000\n"
        "Density of all electrons:   0.01234500\n"
        "Laplacian of electron density:   0.05678900\n"
        "\n"
        "CP  4 (3,+3)\n"
        "Position (Bohr):   1.000000   1.000000   1.000000\n"
        "Density of all electrons:   0.00123400\n"
        "Laplacian of electron density:   0.00987600\n"
    )


@pytest.fixture()
def sample_spectrum_output() -> str:
    return (
        "IR spectrum data\n"
        "  500.0 cm^-1  Intensity:  12.34\n"
        "  800.0 cm^-1  Intensity:  45.67\n"
        " 1200.0 cm^-1  Intensity:  89.01\n"
        " 1500.0 cm^-1  Intensity:  23.45\n"
        " 3000.0 cm^-1  Intensity: 156.78\n"
    )


# =============================================================================
# Tests: ChargeParser
# =============================================================================


class TestChargeParser:
    """Tests for ChargeParser."""

    # ---- Positive ----

    def test_parse_hirshfeld_values(self) -> None:
        """Pattern 1: 'Hirshfeld charge of atom N(X) is VALUE'."""
        output = (
            "Hirshfeld charge analysis\n"
            "Hirshfeld charge of atom     1(C ) is  -0.05230000\n"
            "Hirshfeld charge of atom     2(C ) is   0.12340000\n"
            "Hirshfeld charge of atom     3(H ) is   0.04560000\n"
            "Hirshfeld charge of atom     4(H ) is   0.02340000\n"
            "Hirshfeld charge of atom     5(O ) is  -0.32450000\n"
            "Hirshfeld charge of atom     6(N ) is   0.18230000\n"
        )
        charges = ChargeParser.parse(output, method="Hirshfeld")
        assert len(charges) == 6
        assert charges[1] == pytest.approx(-0.0523)
        assert charges[5] == pytest.approx(-0.3245)
        assert charges[6] == pytest.approx(0.1823)

    def test_parse_mulliken_final_section(self) -> None:
        """Pattern 2 inside 'Final atomic charges' section."""
        output = (
            "Mulliken population analysis\n"
            "Mulliken charge section\n"
            "Final atomic charges:\n"
            "Atom    1(C ):    -0.12340000\n"
            "Atom    2(C ):     0.23450000\n"
            "Atom    5(O ):    -0.45670000\n"
        )
        charges = ChargeParser.parse(output, method="Mulliken")
        assert charges[1] == pytest.approx(-0.1234)
        assert charges[5] == pytest.approx(-0.4567)

    def test_parse_case_insensitive(self) -> None:
        output = "Hirshfeld charge of atom     1(C ) is  -0.05230000\n"
        c1 = ChargeParser.parse(output, method="hirshfeld")
        c2 = ChargeParser.parse(output, method="HIRSHFELD")
        assert c1 == c2
        assert 1 in c1

    def test_parse_population_format(self) -> None:
        """Pattern 5: 'Population of atom N(X) : VALUE'."""
        output = (
            "Mulliken population\n"
            "Population of atom  1(C ) :   5.9679\n"
            "Population of atom  2(N ) :   7.0321\n"
        )
        charges = ChargeParser.parse(output, method="Mulliken")
        assert charges[1] == pytest.approx(5.9679)
        assert charges[2] == pytest.approx(7.0321)

    def test_parse_generic_charge_of_atom(self) -> None:
        """Pattern 6: 'Charge of atom N(X) : VALUE'."""
        output = (
            "ADCH charge\n"
            "Charge of atom  1(C ) :   0.0321\n"
            "Charge of atom  2(N ) :  -0.1500\n"
        )
        charges = ChargeParser.parse(output, method="ADCH")
        assert charges[1] == pytest.approx(0.0321)
        assert charges[2] == pytest.approx(-0.15)

    def test_parse_table_format(self) -> None:
        """Pattern 3: '  1(C )  -0.0523' summary table."""
        output = (
            "Summary of RESP charge\n"
            "  1(C )  -0.05230\n"
            "  2(N )   0.12340\n"
        )
        charges = ChargeParser.parse(output, method="RESP")
        assert charges[1] == pytest.approx(-0.0523)
        assert charges[2] == pytest.approx(0.1234)

    def test_parse_column_format(self) -> None:
        """Pattern 4: '  1  C  -0.052300' column format."""
        output = (
            "EEM charge\n"
            "  1  C  -0.052300\n"
            "  2  N   0.123400\n"
        )
        charges = ChargeParser.parse(output, method="EEM")
        assert charges[1] == pytest.approx(-0.0523)
        assert charges[2] == pytest.approx(0.1234)

    def test_final_charges_override_earlier(self) -> None:
        """'Final atomic charges' clears and replaces previous values."""
        output = (
            "Hirshfeld charge of atom     1(C ) is   0.11111111\n"
            "Final atomic charges:\n"
            "Atom    1(C ):     0.22222222\n"
        )
        charges = ChargeParser.parse(output, method="Hirshfeld")
        assert charges[1] == pytest.approx(0.22222222)

    def test_summary_section_clears_previous(self) -> None:
        """'Summary of … charge' clears prior charges."""
        output = (
            "Hirshfeld charge of atom     1(C ) is   0.11111111\n"
            "Summary of Hirshfeld charge\n"
            "  1(C )   0.22222222\n"
        )
        charges = ChargeParser.parse(output, method="Hirshfeld")
        assert charges[1] == pytest.approx(0.22222222)

    def test_section_end_detection(self) -> None:
        """'calculation took' ends the final section; later lines ignored."""
        output = (
            "Final atomic charges:\n"
            "Atom    1(C ):     0.12340000\n"
            "Atom    2(N ):    -0.56780000\n"
            "calculation took 5.23 seconds\n"
            "Atom    3(O ):     9.99990000\n"
        )
        charges = ChargeParser.parse(output, method="Hirshfeld")
        assert 1 in charges
        assert 2 in charges
        assert 3 not in charges

    def test_parse_dipole(self) -> None:
        output = (
            "Dipole moment (Debye):  "
            "X=   1.2345  Y=  -0.5678  Z=   0.9012  Tot=   1.6234\n"
        )
        dipole = ChargeParser.parse_dipole(output)
        assert dipole is not None
        assert dipole["x"] == pytest.approx(1.2345)
        assert dipole["y"] == pytest.approx(-0.5678)
        assert dipole["z"] == pytest.approx(0.9012)
        assert dipole["total"] == pytest.approx(1.6234)

    # ---- Negative ----

    def test_parse_empty_output(self) -> None:
        assert ChargeParser.parse("", method="Hirshfeld") == {}

    def test_parse_no_charges_section(self) -> None:
        assert ChargeParser.parse("Random text\nNothing here\n", method="Hirshfeld") == {}

    def test_parse_dipole_missing(self) -> None:
        assert ChargeParser.parse_dipole("No dipole here") is None

    def test_wrong_method_no_match(self) -> None:
        """Pattern 1 requires method name match."""
        output = "ADCH charge of atom     1(C ) is   0.12340000\n"
        charges = ChargeParser.parse(output, method="Hirshfeld")
        assert 1 not in charges

    def test_pattern2_requires_section_flag(self) -> None:
        """Pattern 2 lines outside a charge section are skipped."""
        output = "Atom    1(C ):     0.12340000\n"
        assert ChargeParser.parse(output, method="Hirshfeld") == {}

    def test_irrelevant_lines_skipped(self) -> None:
        output = (
            "Hirshfeld charge analysis\n"
            "=== some divider ===\n"
            "Memory: 128 MB\n"
            "Hirshfeld charge of atom     1(C ) is  -0.05230000\n"
            "Running on 4 threads\n"
        )
        charges = ChargeParser.parse(output, method="Hirshfeld")
        assert len(charges) == 1
        assert charges[1] == pytest.approx(-0.0523)

    # ---- Edge cases ----

    def test_scientific_notation(self) -> None:
        output = (
            "Hirshfeld charge analysis\n"
            "Hirshfeld charge of atom     1(C ) is  -5.23E-02\n"
            "Hirshfeld charge of atom     2(H ) is   1.00e+00\n"
        )
        charges = ChargeParser.parse(output, method="Hirshfeld")
        assert charges[1] == pytest.approx(-0.0523)
        assert charges[2] == pytest.approx(1.0)

    def test_zero_charge(self) -> None:
        output = (
            "Hirshfeld charge analysis\n"
            "Hirshfeld charge of atom     1(C ) is   0.00000000\n"
        )
        assert ChargeParser.parse(output, method="Hirshfeld")[1] == pytest.approx(0.0)

    def test_single_atom(self) -> None:
        output = (
            "Hirshfeld charge analysis\n"
            "Hirshfeld charge of atom     1(He) is   0.00000000\n"
        )
        assert len(ChargeParser.parse(output, method="Hirshfeld")) == 1

    def test_large_atom_index(self) -> None:
        output = "Final atomic charges:\nAtom  999(Zn):     0.98760000\n"
        assert ChargeParser.parse(output, method="Hirshfeld")[999] == pytest.approx(0.9876)

    def test_method_with_special_regex_chars(self) -> None:
        """re.escape handles 'Hirshfeld-I'."""
        output = "Hirshfeld-I charge of atom     1(C ) is   0.05000000\n"
        assert ChargeParser.parse(output, method="Hirshfeld-I")[1] == pytest.approx(0.05)

    def test_dipole_all_negative_components(self) -> None:
        output = (
            "Dipole moment (Debye):  "
            "X=  -1.0000  Y=  -2.0000  Z=  -3.0000  Tot=   3.7417\n"
        )
        d = ChargeParser.parse_dipole(output)
        assert d is not None
        assert d["x"] < 0 and d["y"] < 0 and d["z"] < 0
        assert d["total"] > 0

    def test_whitespace_only_input(self) -> None:
        assert ChargeParser.parse("   \n\n   \n", method="Hirshfeld") == {}

    def test_dipole_whitespace_only(self) -> None:
        assert ChargeParser.parse_dipole("   \n") is None


# =============================================================================
# Tests: OrbitalCompositionParser
# =============================================================================


class TestOrbitalCompositionParser:
    """Tests for OrbitalCompositionParser."""

    # ---- Positive ----

    def test_parse_composition(self) -> None:
        output = (
            "Orbital    5  Occ= 2.000000  E= -0.72340 a.u.\n"
            "  C  1    s :   23.45%\n"
            "  O  3    s :   12.34%\n"
        )
        orbs = OrbitalCompositionParser.parse(output)
        assert len(orbs) == 1
        assert orbs[0]["orbital"] == 5
        assert orbs[0]["occupation"] == pytest.approx(2.0)
        assert orbs[0]["energy"] == pytest.approx(-0.7234)
        assert orbs[0]["contributions"]["C1"] == pytest.approx(23.45)
        assert orbs[0]["contributions"]["O3"] == pytest.approx(12.34)

    def test_parse_multiple_orbitals(self) -> None:
        output = (
            "Orbital    5  Occ= 2.000000  E= -0.72340 a.u.\n"
            "  C  1    s :   23.45%\n"
            "\n"
            "Orbital    6  Occ= 2.000000  E= -0.51230 a.u.\n"
            "  N  2    p :   65.43%\n"
        )
        orbs = OrbitalCompositionParser.parse(output)
        assert len(orbs) == 2
        assert orbs[1]["orbital"] == 6
        assert orbs[1]["contributions"]["N2"] == pytest.approx(65.43)

    def test_parse_pattern2_index_format(self) -> None:
        """'  1(C )    45.23%' format."""
        output = (
            "Orbital   10  Occ= 2.000000  E= -0.30000 a.u.\n"
            "  1(C )   45.23%\n"
            "  2(N )   30.00%\n"
        )
        orbs = OrbitalCompositionParser.parse(output)
        assert len(orbs) == 1
        assert orbs[0]["contributions"]["atom_1"] == pytest.approx(45.23)
        assert orbs[0]["contributions"]["atom_2"] == pytest.approx(30.0)

    def test_parse_oxidation_states(self) -> None:
        output = (
            "Atom   1(C ) formal oxidation state:  0.0\n"
            "Atom   2(N ) formal oxidation state: -1.0\n"
            "Atom   3(O ) formal oxidation state: -2.0\n"
            "Atom   4(Fe) formal oxidation state:  3.0\n"
        )
        states = OrbitalCompositionParser.parse_oxidation_states(output)
        assert states == {1: 0, 2: -1, 3: -2, 4: 3}

    # ---- Negative ----

    def test_parse_empty(self) -> None:
        assert OrbitalCompositionParser.parse("") == []

    def test_parse_oxidation_states_empty(self) -> None:
        assert OrbitalCompositionParser.parse_oxidation_states("") == {}

    def test_contribution_lines_without_header_ignored(self) -> None:
        output = "  C  1    s :   23.45%\n  N  2    p :   65.43%\n"
        assert OrbitalCompositionParser.parse(output) == []

    # ---- Edge cases ----

    def test_single_orbital_no_contributions(self) -> None:
        output = "Orbital    1  Occ= 2.000000  E= -20.00000 a.u.\n"
        orbs = OrbitalCompositionParser.parse(output)
        assert len(orbs) == 1
        assert orbs[0]["contributions"] == {}

    def test_fractional_occupation(self) -> None:
        output = "Orbital    5  Occ= 1.500000  E= -0.50000 a.u.\n"
        assert OrbitalCompositionParser.parse(output)[0]["occupation"] == pytest.approx(1.5)

    def test_oxidation_state_rounding(self) -> None:
        """2.9 rounds to 3."""
        output = "Atom   1(Fe) formal oxidation state:  2.9\n"
        assert OrbitalCompositionParser.parse_oxidation_states(output)[1] == 3

    def test_positive_energy_virtual_orbital(self) -> None:
        output = "Orbital   50  Occ= 0.000000  E=  5.12340 a.u.\n"
        orbs = OrbitalCompositionParser.parse(output)
        assert orbs[0]["energy"] > 0
        assert orbs[0]["occupation"] == pytest.approx(0.0)


# =============================================================================
# Tests: BondOrderParser
# =============================================================================


class TestBondOrderParser:
    """Tests for BondOrderParser."""

    # ---- Positive ----

    def test_parse_mayer_pattern1(self, sample_mayer_output: str) -> None:
        bonds = BondOrderParser.parse(sample_mayer_output)
        assert len(bonds) == 5
        assert bonds[(1, 2)] == pytest.approx(1.4523)
        assert bonds[(2, 5)] == pytest.approx(1.8934)
        assert bonds[(5, 6)] == pytest.approx(0.8234)

    def test_parse_wiberg(self, sample_wiberg_output: str) -> None:
        bonds = BondOrderParser.parse(sample_wiberg_output)
        assert len(bonds) == 3
        assert bonds[(1, 2)] == pytest.approx(1.5234)

    def test_parse_pattern2_dash(self) -> None:
        """'1(C ) -    2(N ):  1.4523' format."""
        output = "1(C ) -    2(N ):  1.4523\n3(H ) -    1(C ):  0.9234\n"
        bonds = BondOrderParser.parse(output)
        assert bonds[(1, 2)] == pytest.approx(1.4523)
        assert (1, 3) in bonds

    def test_parse_pattern3_simple_dash(self) -> None:
        """'  1 - 2   1.4523' format."""
        output = "  1 - 2   1.45230\n"
        bonds = BondOrderParser.parse(output)
        assert bonds[(1, 2)] == pytest.approx(1.4523)

    def test_parse_pattern4_column(self) -> None:
        """'  1    2    1.4523' format."""
        output = "  1    2    1.45230\n"
        bonds = BondOrderParser.parse(output)
        assert bonds[(1, 2)] == pytest.approx(1.4523)

    def test_tuple_ordering(self, sample_mayer_output: str) -> None:
        for a, b in BondOrderParser.parse(sample_mayer_output):
            assert a < b

    def test_reverse_atom_order_normalized(self) -> None:
        output = "#    1:         6(N )    1(C )    0.82340000\n"
        bonds = BondOrderParser.parse(output)
        assert (1, 6) in bonds
        assert (6, 1) not in bonds

    def test_parse_valence(self) -> None:
        output = (
            "Total valence of atom    1(C ):   3.94120000\n"
            "Total valence of atom    2(N ):   2.98760000\n"
            "Free valence of atom     1(C ):   0.05880000\n"
            "Free valence of atom     2(N ):   0.01240000\n"
        )
        v = BondOrderParser.parse_valence(output)
        assert v[1]["total_valence"] == pytest.approx(3.9412)
        assert v[1]["free_valence"] == pytest.approx(0.0588)
        assert v[2]["total_valence"] == pytest.approx(2.9876)

    def test_parse_multicenter(self) -> None:
        output = (
            "Multi-center bond order of atoms  1  2  3 :  0.12345\n"
            "Multi-center bond order of atoms  1  2  3  4  5  6 :  0.06789\n"
        )
        results = BondOrderParser.parse_multicenter(output)
        assert len(results) == 2
        assert results[0]["atoms"] == [1, 2, 3]
        assert results[0]["bond_order"] == pytest.approx(0.12345)
        assert len(results[1]["atoms"]) == 6

    def test_parse_decomposition(self) -> None:
        output = (
            "Orbital   5:   0.23456\n"
            "Orbital   6:   0.45678\n"
            "Orbital   7:  -0.01234\n"
        )
        decomp = BondOrderParser.parse_decomposition(output)
        assert len(decomp) == 3
        assert decomp[0] == {"orbital": 5, "contribution": pytest.approx(0.23456)}
        assert decomp[2]["contribution"] == pytest.approx(-0.01234)

    # ---- Negative ----

    def test_parse_empty(self) -> None:
        assert BondOrderParser.parse("") == {}

    def test_parse_valence_empty(self) -> None:
        assert BondOrderParser.parse_valence("") == {}

    def test_parse_multicenter_empty(self) -> None:
        assert BondOrderParser.parse_multicenter("") == []

    def test_parse_decomposition_empty(self) -> None:
        assert BondOrderParser.parse_decomposition("") == []

    def test_irrelevant_text(self) -> None:
        assert BondOrderParser.parse("=== Header ===\nMemory: 256 MB\n") == {}

    # ---- Edge cases ----

    def test_valence_partial_no_free(self) -> None:
        output = "Total valence of atom    1(C ):   3.94120000\n"
        v = BondOrderParser.parse_valence(output)
        assert "total_valence" in v[1]
        assert "free_valence" not in v[1]

    def test_very_small_bond_order(self) -> None:
        output = "#    1:         1(C )    2(C )    0.00001234\n"
        assert BondOrderParser.parse(output)[(1, 2)] == pytest.approx(1.234e-5)

    def test_decomposition_all_negative(self) -> None:
        output = "Orbital  10:  -0.50000\nOrbital  11:  -0.30000\n"
        decomp = BondOrderParser.parse_decomposition(output)
        assert all(d["contribution"] < 0 for d in decomp)

    def test_whitespace_only_input(self) -> None:
        assert BondOrderParser.parse("   \n\n") == {}


# =============================================================================
# Tests: CriticalPointParser
# =============================================================================


class TestCriticalPointParser:
    """Tests for CriticalPointParser."""

    # ---- Positive ----

    def test_parse_all_cp_types(self, sample_topology_output: str) -> None:
        cps = CriticalPointParser.parse(sample_topology_output)
        assert len(cps) == 4

    def test_nuclear_cp(self, sample_topology_output: str) -> None:
        cps = CriticalPointParser.parse(sample_topology_output)
        nuc = [c for c in cps if c["type"] == "(3,-3)"][0]
        assert nuc["cp_type"] == "nuclear"
        assert nuc["position"] == (0.0, 0.0, 0.0)

    def test_bond_cp(self, sample_topology_output: str) -> None:
        bond = [c for c in CriticalPointParser.parse(sample_topology_output) if c["type"] == "(3,-1)"][0]
        assert bond["cp_type"] == "bond"
        assert bond["position"][0] == pytest.approx(1.234567)

    def test_ring_cp(self, sample_topology_output: str) -> None:
        ring = [c for c in CriticalPointParser.parse(sample_topology_output) if c["type"] == "(3,+1)"][0]
        assert ring["cp_type"] == "ring"

    def test_cage_cp(self, sample_topology_output: str) -> None:
        cage = [c for c in CriticalPointParser.parse(sample_topology_output) if c["type"] == "(3,+3)"][0]
        assert cage["cp_type"] == "cage"

    def test_rho_laplacian_ellipticity(self, sample_topology_output: str) -> None:
        cps = CriticalPointParser.parse(sample_topology_output)
        assert cps[0]["rho"] == pytest.approx(0.298765)
        assert cps[0]["laplacian"] == pytest.approx(-1.123456)
        assert cps[1]["ellipticity"] == pytest.approx(0.04523)

    def test_parse_bond_paths(self) -> None:
        output = (
            "Bond path between atom  1(C ) and atom  2(N ), BCP  3, length  2.456\n"
            "Bond path between atom  2(N ) and atom  5(O ), BCP  7, length  3.123\n"
        )
        paths = CriticalPointParser.parse_bond_paths(output)
        assert len(paths) == 2
        assert paths[0] == {"atom1": 1, "atom2": 2, "bcp_index": 3, "path_length": pytest.approx(2.456)}
        assert paths[1]["path_length"] == pytest.approx(3.123)

    def test_summary(self, sample_topology_output: str) -> None:
        cps = CriticalPointParser.parse(sample_topology_output)
        counts = CriticalPointParser.summary(cps)
        assert counts == {"nuclear": 1, "bond": 1, "ring": 1, "cage": 1}

    # ---- Negative ----

    def test_parse_empty(self) -> None:
        assert CriticalPointParser.parse("") == []

    def test_bond_paths_empty(self) -> None:
        assert CriticalPointParser.parse_bond_paths("") == []

    def test_summary_empty(self) -> None:
        assert CriticalPointParser.summary([]) == {}

    def test_cp_without_position(self) -> None:
        """CP header without position line yields position=None."""
        cps = CriticalPointParser.parse("CP  1 (3,-3)\nSome other line\n")
        assert len(cps) == 1
        assert cps[0]["position"] is None

    # ---- Edge cases ----

    def test_ring_cp_no_ellipticity(self, sample_topology_output: str) -> None:
        """Ring CP has no ellipticity line → None."""
        ring = [c for c in CriticalPointParser.parse(sample_topology_output) if c["type"] == "(3,+1)"][0]
        assert ring["ellipticity"] is None

    def test_single_cp(self) -> None:
        output = (
            "CP  1 (3,-3)\n"
            "Position (Bohr):   0.000000   0.000000   0.000000\n"
            "Density of all electrons:   0.10000000\n"
        )
        cps = CriticalPointParser.parse(output)
        assert len(cps) == 1
        assert cps[0]["rho"] == pytest.approx(0.1)
        assert cps[0]["laplacian"] is None

    def test_unknown_cp_type(self) -> None:
        cps = CriticalPointParser.parse("CP  1 (2,-1)\nPosition (Bohr):   0.0   0.0   0.0\n")
        assert cps[0]["cp_type"] == "unknown"

    def test_summary_all_same_type(self) -> None:
        cps = [{"cp_type": "bond"}, {"cp_type": "bond"}, {"cp_type": "bond"}]
        assert CriticalPointParser.summary(cps) == {"bond": 3}

    def test_positive_laplacian(self, sample_topology_output: str) -> None:
        """Ring and cage CPs have positive Laplacian."""
        ring = [c for c in CriticalPointParser.parse(sample_topology_output) if c["type"] == "(3,+1)"][0]
        assert ring["laplacian"] > 0


# =============================================================================
# Tests: DOSParser
# =============================================================================


class TestDOSParser:
    """Tests for DOSParser."""

    # ---- Positive ----

    def test_parse_tdos(self) -> None:
        output = (
            "TDOS data\n"
            "  -15.0000    0.1234\n"
            "  -14.5000    0.2345\n"
            "  -14.0000    0.5678\n"
            "  -13.5000    0.8901\n"
        )
        data = DOSParser.parse(output)
        assert len(data["energies"]) == 4
        assert data["energies"][0] == pytest.approx(-15.0)
        assert data["dos"][0] == pytest.approx(0.1234)

    def test_parse_pdos_multi_column(self) -> None:
        output = (
            "PDOS data\n"
            "  -15.0000    0.1234    0.0567    0.0345\n"
            "  -14.5000    0.2345    0.1234    0.0678\n"
        )
        data = DOSParser.parse(output)
        assert len(data["energies"]) == 2
        assert "pdos_1" in data
        assert "pdos_2" in data
        assert data["pdos_1"][0] == pytest.approx(0.0567)

    def test_parse_orbital_energies(self) -> None:
        output = (
            "   5   -19.234  eV  Occ=  2.000000\n"
            "   6   -15.678  eV  Occ=  2.000000\n"
            "   7    -0.567  eV  Occ=  0.000000\n"
        )
        orbs = DOSParser.parse_orbital_energies(output)
        assert len(orbs) == 3
        assert orbs[0] == {"index": 5, "energy_eV": pytest.approx(-19.234), "occupation": pytest.approx(2.0)}

    # ---- Negative ----

    def test_parse_empty(self) -> None:
        data = DOSParser.parse("")
        assert data["energies"] == []
        assert data["dos"] == []

    def test_orbital_energies_empty(self) -> None:
        assert DOSParser.parse_orbital_energies("") == []

    def test_data_before_header_ignored(self) -> None:
        output = "  -15.0000    0.1234\nTDOS data\n  -14.0000    0.5678\n"
        data = DOSParser.parse(output)
        assert len(data["energies"]) == 1
        assert data["energies"][0] == pytest.approx(-14.0)

    def test_pdos_key_absent_for_two_column(self) -> None:
        output = "TDOS data\n  -10.0000    0.5000\n"
        data = DOSParser.parse(output)
        assert "pdos_1" not in data

    # ---- Edge cases ----

    def test_negative_dos_value_opdos(self) -> None:
        output = "OPDOS data\n  -10.0000   -0.5000\n"
        assert DOSParser.parse(output)["dos"][0] == pytest.approx(-0.5)

    def test_single_data_line(self) -> None:
        output = "TDOS data\n  -10.0000    0.5000\n"
        data = DOSParser.parse(output)
        assert len(data["energies"]) == 1

    def test_orbital_energies_virtual(self) -> None:
        output = "  50     5.123  eV  Occ=  0.000000\n"
        orbs = DOSParser.parse_orbital_energies(output)
        assert orbs[0]["energy_eV"] > 0
        assert orbs[0]["occupation"] == pytest.approx(0.0)


# =============================================================================
# Tests: SpectrumParser
# =============================================================================


class TestSpectrumParser:
    """Tests for SpectrumParser."""

    # ---- Positive ----

    def test_parse_ir(self, sample_spectrum_output: str) -> None:
        sp = SpectrumParser.parse(sample_spectrum_output)
        assert len(sp["frequencies"]) == 5
        assert sp["frequencies"][0] == pytest.approx(500.0)
        assert sp["intensities"][-1] == pytest.approx(156.78)

    def test_parse_uv_vis(self) -> None:
        output = "  345.67 nm  f= 0.1234\n  280.12 nm  f= 0.5678\n"
        sp = SpectrumParser.parse(output)
        assert "wavelengths" in sp
        assert sp["wavelengths"][0] == pytest.approx(345.67)
        assert sp["intensities"][0] == pytest.approx(0.1234)

    def test_parse_nmr(self) -> None:
        output = (
            "Atom  1(C ) shift:  123.45 ppm\n"
            "Atom  2(C ) shift:   45.67 ppm\n"
            "Atom  3(H ) shift:    7.89 ppm\n"
        )
        sp = SpectrumParser.parse(output)
        assert "chemical_shifts" in sp
        assert len(sp["chemical_shifts"]) == 3
        assert sp["chemical_shifts"][0] == pytest.approx(123.45)

    def test_parse_transitions(self) -> None:
        output = (
            "Excited state   1:  E= 3.4567 eV  lam= 358.7 nm  f= 0.0123\n"
            "Excited state   2:  E= 4.1234 eV  lam= 300.5 nm  f= 0.5678\n"
        )
        t = SpectrumParser.parse_transitions(output)
        assert len(t) == 2
        assert t[0]["state"] == 1
        assert t[0]["energy_eV"] == pytest.approx(3.4567)
        assert t[0]["wavelength_nm"] == pytest.approx(358.7)
        assert t[0]["osc_strength"] == pytest.approx(0.0123)

    def test_parse_color_full(self) -> None:
        output = "CIE: X=  0.3456  Y=  0.3210  Z=  0.2890  RGB: R= 180  G= 120  B=  90\n"
        c = SpectrumParser.parse_color(output)
        assert c is not None
        assert c["X"] == pytest.approx(0.3456)
        assert c["R"] == 180
        assert c["B"] == 90

    # ---- Negative ----

    def test_parse_empty(self) -> None:
        assert SpectrumParser.parse("") == {"frequencies": [], "intensities": []}

    def test_transitions_empty(self) -> None:
        assert SpectrumParser.parse_transitions("") == []

    def test_color_missing(self) -> None:
        assert SpectrumParser.parse_color("No color data") is None

    def test_uv_keys_absent_for_ir(self, sample_spectrum_output: str) -> None:
        sp = SpectrumParser.parse(sample_spectrum_output)
        assert "wavelengths" not in sp
        assert "chemical_shifts" not in sp

    def test_nmr_keys_absent_for_ir(self, sample_spectrum_output: str) -> None:
        sp = SpectrumParser.parse(sample_spectrum_output)
        assert "atom_indices" not in sp

    def test_ir_keys_absent_for_uv(self) -> None:
        sp = SpectrumParser.parse("  345.67 nm  f= 0.1234\n")
        assert sp.get("frequencies") is None or sp.get("frequencies") == []

    # ---- Edge cases ----

    def test_ir_zero_intensity(self) -> None:
        sp = SpectrumParser.parse("  1500.0 cm^-1  Intensity:  0.00\n")
        assert sp["intensities"][0] == pytest.approx(0.0)

    def test_single_transition(self) -> None:
        output = "Excited state   1:  E= 3.0000 eV  lam= 413.3 nm  f= 0.0001\n"
        assert len(SpectrumParser.parse_transitions(output)) == 1

    def test_color_rgb_only(self) -> None:
        c = SpectrumParser.parse_color("R= 255  G= 0  B= 128\n")
        assert c is not None
        assert c["R"] == 255
        assert "X" not in c

    def test_generic_two_column_fallback(self) -> None:
        sp = SpectrumParser.parse("  100.00    5.67\n  200.00   10.23\n")
        assert len(sp["frequencies"]) == 2

    def test_dark_transition_low_f(self) -> None:
        output = "Excited state   1:  E= 5.0000 eV  lam= 248.0 nm  f= 0.0000\n"
        t = SpectrumParser.parse_transitions(output)
        assert t[0]["osc_strength"] == pytest.approx(0.0)


# =============================================================================
# Tests: SurfaceParser
# =============================================================================


class TestSurfaceParser:
    """Tests for SurfaceParser."""

    # ---- Positive ----

    def test_parse_all_descriptors(self) -> None:
        output = (
            "Overall surface area:  234.5678\n"
            "Enclosed volume:  345.6789\n"
            "V_S+ :   15.4567\n"
            "V_S- :  -23.4567\n"
            "sigma^2_tot :   0.01234\n"
            "nu :   0.23456\n"
            "Pi :   1.23456\n"
            "Balance :   0.34567\n"
            "Global surface maximum:   56.7890\n"
            "Global surface minimum:  -45.6789\n"
        )
        r = SurfaceParser.parse(output)
        assert r["area"] == pytest.approx(234.5678)
        assert r["volume"] == pytest.approx(345.6789)
        assert r["V_S_plus"] == pytest.approx(15.4567)
        assert r["V_S_minus"] == pytest.approx(-23.4567)
        assert r["sigma2_total"] == pytest.approx(0.01234)
        assert r["nu"] == pytest.approx(0.23456)
        assert r["V_S_max"] == pytest.approx(56.789)
        assert r["V_S_min"] == pytest.approx(-45.6789)

    def test_parse_extrema(self) -> None:
        output = (
            "Local minimum  1:   -45.678  at   1.234   2.345   3.456\n"
            "Local maximum  1:    56.789  at   7.890   8.901   9.012\n"
        )
        ext = SurfaceParser.parse_extrema(output)
        mins = [e for e in ext if e["type"] == "min"]
        maxs = [e for e in ext if e["type"] == "max"]
        assert len(mins) == 1 and len(maxs) == 1
        assert mins[0]["value"] == pytest.approx(-45.678)
        assert mins[0]["position"] == pytest.approx((1.234, 2.345, 3.456), abs=1e-3)

    # ---- Negative ----

    def test_parse_empty(self) -> None:
        assert SurfaceParser.parse("") == {}

    def test_extrema_empty(self) -> None:
        assert SurfaceParser.parse_extrema("") == []

    def test_unrelated_text(self) -> None:
        assert SurfaceParser.parse("Running calculation...\n") == {}

    # ---- Edge cases ----

    def test_partial_descriptors(self) -> None:
        r = SurfaceParser.parse("Overall surface area:  100.000\n")
        assert "area" in r
        assert "volume" not in r

    def test_extrema_negative_position(self) -> None:
        output = "Local minimum  1:   -10.000  at  -1.000  -2.000  -3.000\n"
        ext = SurfaceParser.parse_extrema(output)
        assert ext[0]["position"] == pytest.approx((-1.0, -2.0, -3.0), abs=1e-3)

    def test_extrema_single(self) -> None:
        output = "Local maximum  1:    10.000  at   0.0   0.0   0.0\n"
        assert len(SurfaceParser.parse_extrema(output)) == 1


# =============================================================================
# Tests: FuzzySpaceParser
# =============================================================================


class TestFuzzySpaceParser:
    """Tests for FuzzySpaceParser."""

    # ---- Positive ----

    def test_parse_atomic_properties(self) -> None:
        output = (
            "Atom   1(C ): population=  5.96780\n"
            "Atom   2(N ): population=  7.12340\n"
            "Atomic dipole of atom  1(C ):  X= 0.1234  Y= -0.5678  Z= 0.9012\n"
            "Atom  1  volume:  23.456\n"
        )
        a = FuzzySpaceParser.parse_atomic_properties(output)
        assert a[1]["population"] == pytest.approx(5.9678)
        assert a[1]["dipole_x"] == pytest.approx(0.1234)
        assert a[1]["dipole_y"] == pytest.approx(-0.5678)
        assert a[1]["volume"] == pytest.approx(23.456)
        assert a[2]["population"] == pytest.approx(7.1234)

    def test_parse_delocalization_indices(self) -> None:
        output = (
            "Localization index of atom  1(C ):  3.45670\n"
            "Localization index of atom  2(N ):  5.67890\n"
            "Delocalization index of atom  1(C ) and atom  2(N ):  1.23456\n"
            "Delocalization index of atom  1(C ) and atom  3(O ):  0.45678\n"
        )
        r = FuzzySpaceParser.parse_delocalization_indices(output)
        assert r["localization"][1] == pytest.approx(3.4567)
        assert r["delocalization"][(1, 2)] == pytest.approx(1.23456)
        assert r["delocalization"][(1, 3)] == pytest.approx(0.45678)

    def test_parse_aromaticity_index(self) -> None:
        output = "PDI= 0.05678\nFLU= 0.00123\nMCI= 0.04567\nIring= 0.03456\n"
        r = FuzzySpaceParser.parse_aromaticity_index(output)
        assert r["PDI"] == pytest.approx(0.05678)
        assert r["FLU"] == pytest.approx(0.00123)
        assert r["Iring"] == pytest.approx(0.03456)

    # ---- Negative ----

    def test_atomic_properties_empty(self) -> None:
        assert FuzzySpaceParser.parse_atomic_properties("") == {}

    def test_delocalization_empty(self) -> None:
        r = FuzzySpaceParser.parse_delocalization_indices("")
        assert r == {"localization": {}, "delocalization": {}}

    def test_aromaticity_index_empty(self) -> None:
        assert FuzzySpaceParser.parse_aromaticity_index("") == {}

    def test_unrelated_text(self) -> None:
        assert FuzzySpaceParser.parse_atomic_properties("memory=1024\n") == {}

    # ---- Edge cases ----

    def test_di_pair_reordering(self) -> None:
        """DI pair (3,1) stored as (1,3)."""
        output = "Delocalization index of atom  3(O ) and atom  1(C ):  0.50000\n"
        r = FuzzySpaceParser.parse_delocalization_indices(output)
        assert (1, 3) in r["delocalization"]
        assert (3, 1) not in r["delocalization"]

    def test_population_only_no_dipole(self) -> None:
        a = FuzzySpaceParser.parse_atomic_properties("Atom   1(C ): population=  6.0\n")
        assert "population" in a[1]
        assert "dipole_x" not in a[1]

    def test_aromaticity_partial(self) -> None:
        r = FuzzySpaceParser.parse_aromaticity_index("HOMA=  0.98760\n")
        assert "PDI" not in r and "FLU" not in r


# =============================================================================
# Tests: BasinParser
# =============================================================================


class TestBasinParser:
    """Tests for BasinParser."""

    # ---- Positive ----

    def test_parse_basins(self) -> None:
        output = (
            "Basin   1  attractor at atom  1(C )  population:  5.96780\n"
            "Basin   2  attractor at atom  2(N )  population:  7.12340\n"
        )
        b = BasinParser.parse(output)
        assert len(b) == 2
        assert b[0]["basin"] == 1
        assert b[0]["attractor_atom"] == 1
        assert b[0]["attractor_element"] == "C"
        assert b[0]["population"] == pytest.approx(5.9678)

    def test_parse_simple_format(self) -> None:
        output = "basin analysis\nbasin  population\n  1   C    6.01230\n  2   N    7.12340\n"
        b = BasinParser.parse(output)
        assert len(b) == 2
        assert b[0]["attractor_element"] == "C"
        assert b[0]["population"] == pytest.approx(6.0123)

    def test_parse_charges_aim(self) -> None:
        output = "AIM charge of atom  1(C ):  0.03220\nAIM charge of atom  2(N ): -0.12340\n"
        c = BasinParser.parse_charges(output)
        assert c[1] == pytest.approx(0.0322)
        assert c[2] == pytest.approx(-0.1234)

    def test_parse_charges_bader(self) -> None:
        output = "Bader charge of atom  1(C ):  0.50000\n"
        assert BasinParser.parse_charges(output)[1] == pytest.approx(0.5)

    # ---- Negative ----

    def test_parse_empty(self) -> None:
        assert BasinParser.parse("") == []

    def test_charges_empty(self) -> None:
        assert BasinParser.parse_charges("") == {}

    def test_irrelevant_output(self) -> None:
        assert BasinParser.parse("Done.\nMemory freed.\n") == []

    # ---- Edge cases ----

    def test_single_basin(self) -> None:
        output = "Basin   1  attractor at atom  1(He)  population:  2.00000\n"
        assert len(BasinParser.parse(output)) == 1

    def test_two_letter_element(self) -> None:
        output = "Basin   1  attractor at atom  1(Fe)  population:  23.50000\n"
        b = BasinParser.parse(output)
        assert b[0]["attractor_element"] == "Fe"


# =============================================================================
# Tests: ExcitationParser
# =============================================================================


class TestExcitationParser:
    """Tests for ExcitationParser."""

    # ---- Positive ----

    def test_parse_hole_electron(self) -> None:
        output = (
            "D index:  2.34560\n"
            "Sr:  0.56780\n"
            "t index:  1.23450\n"
            "H index:  3.45670\n"
            "E index:  4.56780\n"
            "HDI:  0.12340\n"
            "EDI:  0.23450\n"
            "hole centroid (Angstrom):  1.234  2.345  3.456\n"
            "electron centroid (Angstrom):  4.567  5.678  6.789\n"
        )
        r = ExcitationParser.parse_hole_electron(output)
        assert r["D_index"] == pytest.approx(2.3456)
        assert r["Sr"] == pytest.approx(0.5678)
        assert r["t_index"] == pytest.approx(1.2345)
        assert r["H_index"] == pytest.approx(3.4567)
        assert r["E_index"] == pytest.approx(4.5678)
        assert r["HDI"] == pytest.approx(0.1234)
        assert r["EDI"] == pytest.approx(0.2345)
        assert r["hole_centroid"] == pytest.approx((1.234, 2.345, 3.456), abs=1e-3)
        assert r["electron_centroid"] == pytest.approx((4.567, 5.678, 6.789), abs=1e-3)

    def test_parse_charge_transfer(self) -> None:
        output = (
            "CT distance:  3.45670\n"
            "CT amount:  0.56780\n"
            "Fragment  1:  hole= 0.8500  electron= 0.1500\n"
            "Fragment  2:  hole= 0.1500  electron= 0.8500\n"
        )
        r = ExcitationParser.parse_charge_transfer(output)
        assert r["CT_distance"] == pytest.approx(3.4567)
        assert r["CT_amount"] == pytest.approx(0.5678)
        assert len(r["fragments"]) == 2
        assert r["fragments"][0]["hole"] == pytest.approx(0.85)

    def test_parse_delta_r(self) -> None:
        output = "State   1  Delta_r:  1.23450\nState   2  Delta_r:  4.56780\n"
        r = ExcitationParser.parse_delta_r(output)
        assert len(r) == 2
        assert r[0]["delta_r"] == pytest.approx(1.2345)

    def test_parse_lambda_index(self) -> None:
        output = "State   1  Lambda:  0.78900\nState   2  Lambda:  0.23400\n"
        r = ExcitationParser.parse_lambda_index(output)
        assert len(r) == 2
        assert r[0]["lambda_index"] == pytest.approx(0.789)

    # ---- Negative ----

    def test_hole_electron_empty(self) -> None:
        assert ExcitationParser.parse_hole_electron("") == {}

    def test_charge_transfer_empty(self) -> None:
        assert ExcitationParser.parse_charge_transfer("") == {}

    def test_delta_r_empty(self) -> None:
        assert ExcitationParser.parse_delta_r("") == []

    def test_lambda_index_empty(self) -> None:
        assert ExcitationParser.parse_lambda_index("") == []

    def test_charge_transfer_no_fragments_key(self) -> None:
        """Only distance, no fragment data."""
        r = ExcitationParser.parse_charge_transfer("CT distance:  1.00000\n")
        assert "CT_distance" in r
        assert "fragments" not in r

    # ---- Edge cases ----

    def test_hole_electron_partial(self) -> None:
        r = ExcitationParser.parse_hole_electron("D index:  5.00000\n")
        assert r["D_index"] == pytest.approx(5.0)
        assert "Sr" not in r
        assert "hole_centroid" not in r

    def test_lambda_below_ct_threshold(self) -> None:
        r = ExcitationParser.parse_lambda_index("State   1  Lambda:  0.10000\n")
        assert r[0]["lambda_index"] < 0.3  # CT state

    def test_single_delta_r(self) -> None:
        r = ExcitationParser.parse_delta_r("State   1  Delta_r:  0.00100\n")
        assert len(r) == 1
        assert r[0]["delta_r"] == pytest.approx(0.001)


# =============================================================================
# Tests: WeakInteractionParser
# =============================================================================


class TestWeakInteractionParser:
    """Tests for WeakInteractionParser."""

    # ---- Positive ----

    def test_parse_full(self) -> None:
        output = (
            "delta_g_inter:  0.12340\n"
            "delta_g_intra:  0.56780\n"
            "Integral over isosurface:  0.03456\n"
            "RDG.cube has been generated\n"
            "signlambda2.cube has been generated\n"
        )
        r = WeakInteractionParser.parse(output)
        assert r["delta_g_inter"] == pytest.approx(0.1234)
        assert r["delta_g_intra"] == pytest.approx(0.5678)
        assert r["isosurface_integral"] == pytest.approx(0.03456)
        assert set(r["cube_files"]) == {"RDG.cube", "signlambda2.cube"}

    # ---- Negative ----

    def test_parse_empty(self) -> None:
        assert WeakInteractionParser.parse("") == {}

    def test_no_cube_files_key_absent(self) -> None:
        r = WeakInteractionParser.parse("delta_g_inter:  0.12340\n")
        assert "cube_files" not in r

    # ---- Edge cases ----

    def test_cube_files_only(self) -> None:
        r = WeakInteractionParser.parse("IRI.cube has been generated\n")
        assert r["cube_files"] == ["IRI.cube"]
        assert "delta_g_inter" not in r

    def test_many_cube_files(self) -> None:
        output = "a.cube has been generated\nb.cube has been generated\nc.cube has been generated\n"
        assert len(WeakInteractionParser.parse(output)["cube_files"]) == 3


# =============================================================================
# Tests: EDAParser
# =============================================================================


class TestEDAParser:
    """Tests for EDAParser."""

    # ---- Positive ----

    def test_parse_all_components(self) -> None:
        output = (
            "Electrostatic:  -45.6789\n"
            "Exchange:  -23.4567\n"
            "Pauli repulsion:   67.8901\n"
            "Polarization:  -12.3456\n"
            "Dispersion:   -5.6789\n"
            "Orbital interaction:  -34.5678\n"
            "Total interaction:  -19.3456\n"
        )
        r = EDAParser.parse(output)
        assert r["electrostatic"] == pytest.approx(-45.6789)
        assert r["exchange"] == pytest.approx(-23.4567)
        assert r["repulsion"] == pytest.approx(67.8901)
        assert r["polarization"] == pytest.approx(-12.3456)
        assert r["dispersion"] == pytest.approx(-5.6789)
        assert r["orbital_interaction"] == pytest.approx(-34.5678)
        assert r["total_interaction"] == pytest.approx(-19.3456)

    def test_parse_dispersion_contributions(self) -> None:
        output = "Atom  1(C )  D3 dispersion:  -2.34560\nAtom  2(N )  D3 dispersion:  -1.56780\n"
        c = EDAParser.parse_dispersion_contributions(output)
        assert c[1] == pytest.approx(-2.3456)
        assert c[2] == pytest.approx(-1.5678)

    # ---- Negative ----

    def test_parse_empty(self) -> None:
        assert EDAParser.parse("") == {}

    def test_dispersion_contributions_empty(self) -> None:
        assert EDAParser.parse_dispersion_contributions("") == {}

    def test_unrelated_output(self) -> None:
        assert EDAParser.parse("Running job...\nDone.\n") == {}

    # ---- Edge cases ----

    def test_partial_components(self) -> None:
        r = EDAParser.parse("Electrostatic:  -10.0000\nDispersion:  -1.0000\n")
        assert "electrostatic" in r and "dispersion" in r
        assert "exchange" not in r and "polarization" not in r

    def test_positive_repulsion(self) -> None:
        assert EDAParser.parse("Pauli repulsion:   100.0000\n")["repulsion"] > 0


# =============================================================================
# Tests: CDFTParser
# =============================================================================


class TestCDFTParser:
    """Tests for CDFTParser."""

    # ---- Positive ----

    def test_parse_global_indices(self) -> None:
        output = (
            "Chemical potential:  -0.15234\n"
            "Hardness:   0.21345\n"
            "Softness:   4.68520\n"
            "Electrophilicity:   0.05432\n"
            "Nucleophilicity:   2.34560\n"
            "IP (Ionization potential):   0.28456\n"
            "EA (Electron affinity):   0.07123\n"
        )
        r = CDFTParser.parse_global_indices(output)
        assert r["chemical_potential"] == pytest.approx(-0.15234)
        assert r["hardness"] == pytest.approx(0.21345)
        assert r["softness"] == pytest.approx(4.6852)
        assert r["IP"] == pytest.approx(0.28456)
        assert r["EA"] == pytest.approx(0.07123)

    def test_parse_condensed_fukui(self) -> None:
        output = (
            "Atom   1(C ):  f+= 0.12340  f-= 0.05670  f0= 0.09010\n"
            "Atom   2(N ):  f+= 0.23450  f-= 0.15670  f0= 0.19560\n"
            "Atom   3(O ):  f+= 0.04560  f-= 0.34560  f0= 0.19560\n"
        )
        f = CDFTParser.parse_condensed_fukui(output)
        assert len(f) == 3
        assert f[1]["f_plus"] == pytest.approx(0.1234)
        assert f[1]["f_minus"] == pytest.approx(0.0567)
        assert f[1]["f_zero"] == pytest.approx(0.0901)
        assert f[3]["f_minus"] == pytest.approx(0.3456)

    def test_parse_dual_descriptor(self) -> None:
        output = (
            "Atom   1(C ):  dual= 0.06670\n"
            "Atom   2(N ):  dual= 0.07780\n"
            "Atom   3(O ):  dual=-0.30000\n"
        )
        dd = CDFTParser.parse_dual_descriptor(output)
        assert dd[1] == pytest.approx(0.0667)
        assert dd[3] == pytest.approx(-0.3)

    # ---- Negative ----

    def test_global_indices_empty(self) -> None:
        assert CDFTParser.parse_global_indices("") == {}

    def test_condensed_fukui_empty(self) -> None:
        assert CDFTParser.parse_condensed_fukui("") == {}

    def test_dual_descriptor_empty(self) -> None:
        assert CDFTParser.parse_dual_descriptor("") == {}

    def test_unrelated_text(self) -> None:
        assert CDFTParser.parse_global_indices("Done in 5 sec\n") == {}

    # ---- Edge cases ----

    def test_global_indices_partial(self) -> None:
        r = CDFTParser.parse_global_indices("Chemical potential:  -0.1\nHardness:   0.2\n")
        assert "chemical_potential" in r and "hardness" in r
        assert "softness" not in r and "IP" not in r

    def test_fukui_nucleophilic_site(self) -> None:
        """f+ > f- indicates nucleophilic attack site."""
        f = CDFTParser.parse_condensed_fukui("Atom   1(C ):  f+= 0.5  f-= 0.1  f0= 0.3\n")
        assert f[1]["f_plus"] > f[1]["f_minus"]

    def test_dual_descriptor_sign(self) -> None:
        dd = CDFTParser.parse_dual_descriptor("Atom   1(C ):  dual= 0.50000\nAtom   2(O ):  dual=-0.40000\n")
        assert dd[1] > 0  # nucleophilic
        assert dd[2] < 0  # electrophilic

    def test_condensed_fukui_separate_patterns(self) -> None:
        """Fallback: separate lines for f+, f-, f0."""
        output = (
            "Atom   1(C ):  f+= 0.12340\n"
            "Atom   1(C ):  f-= 0.05670\n"
            "Atom   1(C ):  f0= 0.09010\n"
        )
        f = CDFTParser.parse_condensed_fukui(output)
        assert f[1]["f_plus"] == pytest.approx(0.1234)
        assert f[1]["f_minus"] == pytest.approx(0.0567)
        assert f[1]["f_zero"] == pytest.approx(0.0901)


# =============================================================================
# Tests: PolarizabilityParser
# =============================================================================


class TestPolarizabilityParser:
    """Tests for PolarizabilityParser."""

    # ---- Positive ----

    def test_parse_full(self) -> None:
        output = (
            "Isotropic polarizability:  45.67890\n"
            "Anisotropy:  12.34560\n"
            "alpha_xx:  56.78900\nalpha_xy:   1.23400\nalpha_xz:   0.56780\n"
            "alpha_yy:  42.34500\nalpha_yz:  -0.89010\nalpha_zz:  37.90100\n"
            "Beta total:  123.45600\nGamma average:  456.78900\n"
        )
        r = PolarizabilityParser.parse(output)
        assert r["alpha_iso"] == pytest.approx(45.6789)
        assert r["alpha_aniso"] == pytest.approx(12.3456)
        assert r["alpha_xx"] == pytest.approx(56.789)
        assert r["alpha_yz"] == pytest.approx(-0.8901)
        assert r["beta_total"] == pytest.approx(123.456)
        assert r["gamma_total"] == pytest.approx(456.789)

    # ---- Negative ----

    def test_parse_empty(self) -> None:
        assert PolarizabilityParser.parse("") == {}

    def test_unrelated_text(self) -> None:
        assert PolarizabilityParser.parse("Energy: -76.0 Hartree\n") == {}

    # ---- Edge cases ----

    def test_partial_iso_only(self) -> None:
        r = PolarizabilityParser.parse("Isotropic polarizability:  30.00000\n")
        assert r["alpha_iso"] == pytest.approx(30.0)
        assert "beta_total" not in r and "alpha_xx" not in r

    def test_negative_anisotropy(self) -> None:
        assert PolarizabilityParser.parse("Anisotropy:  -5.0\n")["alpha_aniso"] == pytest.approx(-5.0)

    def test_all_six_tensor_components(self) -> None:
        output = "alpha_xx: 1\nalpha_xy: 2\nalpha_xz: 3\nalpha_yy: 4\nalpha_yz: 5\nalpha_zz: 6\n"
        r = PolarizabilityParser.parse(output)
        assert len([k for k in r if k.startswith("alpha_")]) == 6


# =============================================================================
# Tests: AromaticityParser
# =============================================================================


class TestAromaticityParser:
    """Tests for AromaticityParser."""

    # ---- Positive ----

    def test_parse_indices(self) -> None:
        output = "NICS(0):  -8.12340\nNICS(1):  -10.5678\nNICS_ZZ: -25.67890\nHOMA:   0.98760\nBird index:  92.345\n"
        r = AromaticityParser.parse(output)
        assert r["NICS"] == pytest.approx(-8.1234)
        assert r["NICS_1"] == pytest.approx(-10.5678)
        assert r["NICS_ZZ"] == pytest.approx(-25.6789)
        assert r["HOMA"] == pytest.approx(0.9876)
        assert r["Bird"] == pytest.approx(92.345)

    def test_parse_nics_scan(self) -> None:
        output = (
            "NICS scan data\n"
            "  0.0000   -8.1234\n  0.5000   -9.5678\n"
            "  1.0000  -10.5678\n  2.0000   -6.5678\n"
        )
        d = AromaticityParser.parse_nics_scan(output)
        assert len(d["distances"]) == 4
        assert d["nics_values"][2] == pytest.approx(-10.5678)

    def test_stanger_indices(self) -> None:
        r = AromaticityParser.parse("EN_GEO:  0.12340\nEN_BLA:  0.05670\n")
        assert r["EN_GEO"] == pytest.approx(0.1234)
        assert r["EN_BLA"] == pytest.approx(0.0567)

    # ---- Negative ----

    def test_parse_empty(self) -> None:
        assert AromaticityParser.parse("") == {}

    def test_nics_scan_empty(self) -> None:
        d = AromaticityParser.parse_nics_scan("")
        assert d == {"distances": [], "nics_values": []}

    def test_nics_scan_no_header(self) -> None:
        """Data without 'nics scan' header ignored."""
        d = AromaticityParser.parse_nics_scan("  0.0000   -8.1234\n  1.0000  -10.5678\n")
        assert d["distances"] == []

    # ---- Edge cases ----

    def test_positive_nics_antiaromatic(self) -> None:
        assert AromaticityParser.parse("NICS(0):   15.00000\n")["NICS"] > 0

    def test_negative_homa_antiaromatic(self) -> None:
        assert AromaticityParser.parse("HOMA:  -0.50000\n")["HOMA"] < 0

    def test_partial_indices(self) -> None:
        r = AromaticityParser.parse("HOMA:   0.95000\n")
        assert "HOMA" in r and "NICS" not in r and "Bird" not in r


# =============================================================================
# Tests: WavefunctionParser
# =============================================================================


class TestWavefunctionParser:
    """Tests for WavefunctionParser."""

    # ---- Positive ----

    def test_parse_orbital_info(self) -> None:
        output = (
            "   5   Alpha   Occ= 2.000000   E=  -0.72340 a.u.  -19.684 eV\n"
            "   6   Alpha   Occ= 2.000000   E=  -0.51230 a.u.  -13.941 eV\n"
            "   7   Beta    Occ= 0.000000   E=   0.12340 a.u.    3.358 eV\n"
        )
        o = WavefunctionParser.parse_orbital_info(output)
        assert len(o) == 3
        assert o[0]["index"] == 5
        assert o[0]["spin"] == "alpha"
        assert o[0]["occupation"] == pytest.approx(2.0)
        assert o[0]["energy_au"] == pytest.approx(-0.7234)
        assert o[0]["energy_eV"] == pytest.approx(-19.684)
        assert o[2]["spin"] == "beta"
        assert o[2]["occupation"] == pytest.approx(0.0)

    # ---- Negative ----

    def test_parse_empty(self) -> None:
        assert WavefunctionParser.parse_orbital_info("") == []

    def test_unrelated_text(self) -> None:
        assert WavefunctionParser.parse_orbital_info("Total energy: -76.0\n") == []

    # ---- Edge cases ----

    def test_single_orbital(self) -> None:
        output = "   1   Alpha   Occ= 2.000000   E=  -20.50000 a.u.  -557.715 eV\n"
        o = WavefunctionParser.parse_orbital_info(output)
        assert len(o) == 1
        assert o[0]["energy_eV"] == pytest.approx(-557.715)

    def test_virtual_orbital_positive_energy(self) -> None:
        output = "  10   Alpha   Occ= 0.000000   E=   5.12340 a.u.   139.415 eV\n"
        o = WavefunctionParser.parse_orbital_info(output)
        assert o[0]["energy_au"] > 0


# =============================================================================
# Tests: CubeParser
# =============================================================================


class TestCubeParser:
    """Tests for CubeParser."""

    # ---- Positive ----

    def test_parse_full(self) -> None:
        output = (
            "Grid dimensions: 80 x 80 x 80\n"
            "Minimum value:  0.00001234\n"
            "Maximum value:  0.29876500\n"
            "Mean value:  0.00345670\n"
            "Integral of the grid:  9.99876500\n"
            "Std dev:  0.01234567\n"
            "density.cube has been generated\n"
        )
        r = CubeParser.parse(output)
        assert r["grid_points"] == (80, 80, 80)
        assert r["min"] == pytest.approx(0.00001234)
        assert r["max"] == pytest.approx(0.298765)
        assert r["mean"] == pytest.approx(0.0034567)
        assert r["integral"] == pytest.approx(9.998765)
        assert r["std_dev"] == pytest.approx(0.01234567)
        assert "density.cube" in r["cube_files"]

    def test_multiple_cubes(self) -> None:
        output = "density.cube has been generated\nELF.cube has been generated\n"
        assert len(CubeParser.parse(output)["cube_files"]) == 2

    # ---- Negative ----

    def test_parse_empty(self) -> None:
        assert CubeParser.parse("") == {}

    def test_no_cubes_key_absent(self) -> None:
        r = CubeParser.parse("Grid dimensions: 50 x 50 x 50\n")
        assert "cube_files" not in r

    # ---- Edge cases ----

    def test_cubes_only_no_stats(self) -> None:
        r = CubeParser.parse("ESP.cube has been generated\n")
        assert r["cube_files"] == ["ESP.cube"]
        assert "min" not in r and "grid_points" not in r

    def test_stats_only_no_cubes(self) -> None:
        r = CubeParser.parse("Minimum value:  -1.0\nMaximum value:  1.0\n")
        assert r["min"] == pytest.approx(-1.0)
        assert "cube_files" not in r

    def test_asymmetric_grid(self) -> None:
        assert CubeParser.parse("Grid dimensions: 100 x 80 x 60\n")["grid_points"] == (100, 80, 60)

    def test_negative_minimum(self) -> None:
        assert CubeParser.parse("Minimum value:  -0.50000\n")["min"] == pytest.approx(-0.5)


# =============================================================================
# Tests: UtilityParser
# =============================================================================


class TestUtilityParser:
    """Tests for UtilityParser."""

    # ---- Positive ----

    def test_parse_geometry(self) -> None:
        output = (
            "Bond length between atom  1(C ) and atom  2(N ):  1.3456\n"
            "Bond length between atom  2(N ) and atom  3(O ):  1.4567\n"
            "Angle  1-2-3 :  120.345 degree\n"
            "Dihedral  1-2-3-4 :  -45.678 degree\n"
        )
        r = UtilityParser.parse_geometry(output)
        assert len(r["bond_lengths"]) == 2
        assert r["bond_lengths"][0] == {"atom1": 1, "atom2": 2, "length": pytest.approx(1.3456)}
        assert r["angles"][0]["angle"] == pytest.approx(120.345)
        assert r["angles"][0]["atoms"] == (1, 2, 3)
        assert r["dihedrals"][0]["dihedral"] == pytest.approx(-45.678)
        assert r["dihedrals"][0]["atoms"] == (1, 2, 3, 4)

    def test_parse_multipole_moments(self) -> None:
        output = (
            "Dipole moment:  X=  1.2345  Y= -0.5678  Z=  0.9012\n"
            "Quadrupole moment:  XX= -5.6789  XY=  0.1234\n"
        )
        r = UtilityParser.parse_multipole_moments(output)
        assert r["dipole"]["x"] == pytest.approx(1.2345)
        assert r["dipole"]["y"] == pytest.approx(-0.5678)
        assert r["quadrupole"]["XX"] == pytest.approx(-5.6789)

    def test_parse_coordination_numbers(self) -> None:
        output = (
            "Atom  1(C ) coordination number:  3.00000\n"
            "Atom  2(N ) coordination number:  2.00000\n"
        )
        c = UtilityParser.parse_coordination_numbers(output)
        assert c[1] == pytest.approx(3.0)
        assert c[2] == pytest.approx(2.0)

    def test_parse_bla_boa(self) -> None:
        r = UtilityParser.parse_bla_boa("BLA=  0.05670\nBOA=  0.12340\n")
        assert r["BLA"] == pytest.approx(0.0567)
        assert r["BOA"] == pytest.approx(0.1234)

    def test_parse_generated_files(self) -> None:
        output = (
            "density.cube has been generated\n"
            "new.wfn has been generated\n"
            "output.xyz has been generated\n"
            "spectrum.txt has been generated\n"
        )
        f = UtilityParser.parse_generated_files(output)
        assert set(f) == {"density.cube", "new.wfn", "output.xyz", "spectrum.txt"}

    # ---- Negative ----

    def test_geometry_empty(self) -> None:
        assert UtilityParser.parse_geometry("") == {"bond_lengths": [], "angles": [], "dihedrals": []}

    def test_multipole_empty(self) -> None:
        assert UtilityParser.parse_multipole_moments("") == {}

    def test_coordination_empty(self) -> None:
        assert UtilityParser.parse_coordination_numbers("") == {}

    def test_bla_boa_empty(self) -> None:
        assert UtilityParser.parse_bla_boa("") == {}

    def test_generated_files_empty(self) -> None:
        assert UtilityParser.parse_generated_files("") == []

    def test_unrecognized_extension_ignored(self) -> None:
        assert UtilityParser.parse_generated_files("output.dat has been generated\n") == []

    # ---- Edge cases ----

    def test_geometry_only_bonds(self) -> None:
        r = UtilityParser.parse_geometry("Bond length between atom  1(C ) and atom  2(N ):  1.5\n")
        assert len(r["bond_lengths"]) == 1
        assert r["angles"] == [] and r["dihedrals"] == []

    def test_negative_dihedral(self) -> None:
        r = UtilityParser.parse_geometry("Dihedral  1-2-3-4 :  -179.999 degree\n")
        assert r["dihedrals"][0]["dihedral"] == pytest.approx(-179.999)

    def test_linear_angle(self) -> None:
        r = UtilityParser.parse_geometry("Angle  1-2-3 :  180.000 degree\n")
        assert r["angles"][0]["angle"] == pytest.approx(180.0)

    def test_bla_only(self) -> None:
        r = UtilityParser.parse_bla_boa("BLA=  0.10000\n")
        assert "BLA" in r and "BOA" not in r

    def test_dipole_only_no_quadrupole(self) -> None:
        r = UtilityParser.parse_multipole_moments("Dipole moment:  X=  0.0  Y=  0.0  Z=  1.5\n")
        assert "dipole" in r and "quadrupole" not in r

    def test_all_file_extensions(self) -> None:
        output = (
            "a.cube has been generated\nb.wfn has been generated\n"
            "c.molden has been generated\nd.fch has been generated\n"
            "e.xyz has been generated\nf.pdb has been generated\n"
            "g.txt has been generated\nh.gjf has been generated\n"
        )
        assert len(UtilityParser.parse_generated_files(output)) == 8

    # def test_fractional_coordination_number(self) -> None:
    #     c = UtilityParser.parse_coordination_numbers("Atom  1(Pt) coordination number:  5.50000\n")
    #     assert c[1] == pytest.approx(5.5)