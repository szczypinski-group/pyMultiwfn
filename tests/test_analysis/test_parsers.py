"""Tests for pymultiwfn.analysis.parsers — output parsers.

Covers positive, negative, and edge cases for all parser classes.
"""

import pytest

from pymultiwfn.analysis.parsers import (
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


# =========================================================================
# ChargeParser
# =========================================================================


class TestChargeParser:
    def test_parse_final_charges(self) -> None:
        output = (
            "Final atomic charges:\n"
            "Atom    1(C ):    -0.05230000\n"
            "Atom    2(N ):     0.12340000\n"
        )
        charges = ChargeParser.parse(output)
        assert len(charges) == 2
        assert charges[1] == pytest.approx(-0.0523)
        assert charges[2] == pytest.approx(0.1234)

    def test_parse_pattern1_method_name(self) -> None:
        output = "Hirshfeld charge of atom     1(C ) is  0.03208687\n"
        charges = ChargeParser.parse(output)
        assert charges[1] == pytest.approx(0.03208687)

    def test_final_overrides_earlier(self) -> None:
        output = (
            "Hirshfeld charge of atom     1(C ) is  0.111\n"
            "Final atomic charges:\n"
            "Atom    1(C ):     0.222\n"
        )
        charges = ChargeParser.parse(output)
        assert charges[1] == pytest.approx(0.222)

    def test_empty_output(self) -> None:
        assert ChargeParser.parse("") == {}

    def test_parse_dipole(self) -> None:
        output = (
            "Dipole moment:  X=   1.234  Y=  -0.567  Z=   0.901  Tot=   1.623\n"
        )
        d = ChargeParser.parse_dipole(output)
        assert d is not None
        assert d["total"] == pytest.approx(1.623)

    def test_dipole_not_found(self) -> None:
        assert ChargeParser.parse_dipole("no dipole here") is None

    def test_scientific_notation(self) -> None:
        output = "Final atomic charges:\nAtom    1(C ):    -5.23E-02\n"
        assert ChargeParser.parse(output)[1] == pytest.approx(-0.0523)


# =========================================================================
# BondOrderParser
# =========================================================================


class TestBondOrderParser:
    def test_parse_pattern1(self) -> None:
        output = "#    1:         1(C )    2(N )    1.45230000\n"
        bonds = BondOrderParser.parse(output)
        assert bonds[(1, 2)] == pytest.approx(1.4523)

    def test_tuple_ordering(self) -> None:
        output = "#    1:         6(N )    1(C )    0.82340000\n"
        bonds = BondOrderParser.parse(output)
        assert (1, 6) in bonds
        assert (6, 1) not in bonds

    def test_parse_valence(self) -> None:
        output = (
            "Total valence of atom    1(C ):   3.94120000\n"
            "Free valence of atom     1(C ):   0.05880000\n"
        )
        v = BondOrderParser.parse_valence(output)
        assert v[1]["total_valence"] == pytest.approx(3.9412)
        assert v[1]["free_valence"] == pytest.approx(0.0588)

    def test_parse_multicenter(self) -> None:
        output = "Multi-center bond order of atoms  1  2  3 :  0.12345\n"
        results = BondOrderParser.parse_multicenter(output)
        assert results[0]["atoms"] == [1, 2, 3]

    def test_parse_decomposition(self) -> None:
        output = "Orbital   5:   0.23456\nOrbital   6:  -0.01234\n"
        decomp = BondOrderParser.parse_decomposition(output)
        assert len(decomp) == 2

    def test_empty(self) -> None:
        assert BondOrderParser.parse("") == {}


# =========================================================================
# CriticalPointParser
# =========================================================================


class TestCriticalPointParser:
    def test_parse_all_types(self, sample_topology_output: str) -> None:
        cps = CriticalPointParser.parse(sample_topology_output)
        assert len(cps) == 4
        types = {cp["cp_type"] for cp in cps}
        assert types == {"nuclear", "bond", "ring", "cage"}

    def test_rho_laplacian(self, sample_topology_output: str) -> None:
        cps = CriticalPointParser.parse(sample_topology_output)
        assert cps[0]["rho"] == pytest.approx(0.298765)
        assert cps[1]["ellipticity"] == pytest.approx(0.04523)

    def test_bond_paths(self) -> None:
        output = (
            "Bond path between atom  1(C ) and atom  2(N ), "
            "BCP  3, length  2.456\n"
        )
        paths = CriticalPointParser.parse_bond_paths(output)
        assert paths[0]["atom1"] == 1
        assert paths[0]["path_length"] == pytest.approx(2.456)

    def test_summary(self) -> None:
        cps = [{"cp_type": "bond"}, {"cp_type": "bond"}, {"cp_type": "ring"}]
        assert CriticalPointParser.summary(cps) == {"bond": 2, "ring": 1}

    def test_empty(self) -> None:
        assert CriticalPointParser.parse("") == []

    def test_unknown_type(self) -> None:
        cps = CriticalPointParser.parse("CP  1 (2,-1)\n")
        assert cps[0]["cp_type"] == "unknown"


# =========================================================================
# SpectrumParser
# =========================================================================


class TestSpectrumParser:
    def test_parse_ir(self, sample_spectrum_output: str) -> None:
        sp = SpectrumParser.parse(sample_spectrum_output)
        assert len(sp["frequencies"]) == 5
        assert sp["intensities"][-1] == pytest.approx(156.78)

    def test_parse_uv_vis(self) -> None:
        output = "  345.67 nm  f= 0.1234\n"
        sp = SpectrumParser.parse(output)
        assert "wavelengths" in sp

    def test_parse_nmr(self) -> None:
        output = "Atom  1(C ) shift:  123.45 ppm\n"
        sp = SpectrumParser.parse(output)
        assert "chemical_shifts" in sp

    def test_parse_transitions(self) -> None:
        output = (
            "Excited state   1:  E= 3.4567 eV  lam= 358.7 nm  f= 0.0123\n"
        )
        t = SpectrumParser.parse_transitions(output)
        assert t[0]["state"] == 1
        assert t[0]["energy_eV"] == pytest.approx(3.4567)

    def test_parse_color(self) -> None:
        output = "X=  0.3456  Y=  0.3210  Z=  0.2890\nR= 180  G= 120  B=  90\n"
        c = SpectrumParser.parse_color(output)
        assert c is not None
        assert c["R"] == 180

    def test_empty(self) -> None:
        sp = SpectrumParser.parse("")
        assert sp["frequencies"] == []


# =========================================================================
# DOSParser
# =========================================================================


class TestDOSParser:
    def test_parse_tdos(self) -> None:
        output = "TDOS data\n  -15.0000    0.1234\n  -14.0000    0.5678\n"
        data = DOSParser.parse(output)
        assert len(data["energies"]) == 2

    def test_parse_orbital_energies(self) -> None:
        output = "   5   -19.234  eV  Occ=  2.000000\n"
        orbs = DOSParser.parse_orbital_energies(output)
        assert orbs[0]["index"] == 5

    def test_empty(self) -> None:
        assert DOSParser.parse("")["energies"] == []


# =========================================================================
# SurfaceParser
# =========================================================================


class TestSurfaceParser:
    def test_parse_descriptors(self) -> None:
        output = "Overall surface area:  234.5678\nEnclosed volume:  345.6789\n"
        r = SurfaceParser.parse(output)
        assert r["area"] == pytest.approx(234.5678)

    def test_parse_extrema(self) -> None:
        output = "Local minimum  1:   -45.678  at   1.234   2.345   3.456\n"
        ext = SurfaceParser.parse_extrema(output)
        assert ext[0]["type"] == "min"

    def test_empty(self) -> None:
        assert SurfaceParser.parse("") == {}


# =========================================================================
# Remaining parsers — focused coverage
# =========================================================================


class TestOrbitalCompositionParser:
    def test_parse_composition(self) -> None:
        output = (
            "Orbital    5  Occ= 2.000000  E= -0.72340 a.u.\n"
            "  C  1    s :   23.45%\n"
        )
        orbs = OrbitalCompositionParser.parse(output)
        assert orbs[0]["orbital"] == 5
        assert orbs[0]["contributions"]["C1"] == pytest.approx(23.45)

    def test_parse_oxidation(self) -> None:
        output = "Atom   1(Fe) formal oxidation state:  3.0\n"
        assert OrbitalCompositionParser.parse_oxidation_states(output)[1] == 3

    def test_empty(self) -> None:
        assert OrbitalCompositionParser.parse("") == []


class TestFuzzySpaceParser:
    def test_atomic_properties(self) -> None:
        output = "Atom   1(C ): population=  5.96780\n"
        a = FuzzySpaceParser.parse_atomic_properties(output)
        assert a[1]["population"] == pytest.approx(5.9678)

    def test_delocalization(self) -> None:
        output = "Localization index of atom  1(C ):  3.45670\n"
        r = FuzzySpaceParser.parse_delocalization_indices(output)
        assert r["localization"][1] == pytest.approx(3.4567)

    def test_aromaticity_index(self) -> None:
        output = "PDI= 0.05678\nFLU= 0.00123\n"
        r = FuzzySpaceParser.parse_aromaticity_index(output)
        assert r["PDI"] == pytest.approx(0.05678)


class TestBasinParser:
    def test_parse(self) -> None:
        output = "Basin   1  attractor at atom  1(C )  population:  5.96780\n"
        b = BasinParser.parse(output)
        assert b[0]["population"] == pytest.approx(5.9678)

    def test_parse_charges(self) -> None:
        output = "AIM charge of atom  1(C ):  0.03220\n"
        assert BasinParser.parse_charges(output)[1] == pytest.approx(0.0322)

    def test_empty(self) -> None:
        assert BasinParser.parse("") == []


class TestExcitationParser:
    def test_hole_electron(self) -> None:
        output = "D index:  2.34560\nSr:  0.56780\n"
        r = ExcitationParser.parse_hole_electron(output)
        assert r["D_index"] == pytest.approx(2.3456)

    def test_delta_r(self) -> None:
        output = "State   1  Delta_r:  1.23450\n"
        assert ExcitationParser.parse_delta_r(output)[0]["delta_r"] == pytest.approx(1.2345)

    def test_lambda_index(self) -> None:
        output = "State   1  Lambda:  0.78900\n"
        assert ExcitationParser.parse_lambda_index(output)[0]["lambda_index"] == pytest.approx(0.789)


class TestWeakInteractionParser:
    def test_parse(self) -> None:
        output = "delta_g_inter:  0.12340\nRDG.cube has been generated\n"
        r = WeakInteractionParser.parse(output)
        assert r["delta_g_inter"] == pytest.approx(0.1234)
        assert "RDG.cube" in r["cube_files"]

    def test_empty(self) -> None:
        assert WeakInteractionParser.parse("") == {}


class TestEDAParser:
    def test_parse(self) -> None:
        output = "Electrostatic:  -45.6789\nDispersion:   -5.6789\n"
        r = EDAParser.parse(output)
        assert r["electrostatic"] == pytest.approx(-45.6789)

    def test_empty(self) -> None:
        assert EDAParser.parse("") == {}


class TestCDFTParser:
    def test_global_indices(self) -> None:
        output = "Chemical potential:  -0.15234\nHardness:   0.21345\n"
        r = CDFTParser.parse_global_indices(output)
        assert r["chemical_potential"] == pytest.approx(-0.15234)

    def test_condensed_fukui(self) -> None:
        output = "Atom   1(C ):  f+= 0.1234  f-= 0.0567  f0= 0.0901\n"
        f = CDFTParser.parse_condensed_fukui(output)
        assert f[1]["f_plus"] == pytest.approx(0.1234)

    def test_dual_descriptor(self) -> None:
        output = "Atom   1(C ):  dual= 0.06670\n"
        assert CDFTParser.parse_dual_descriptor(output)[1] == pytest.approx(0.0667)


class TestPolarizabilityParser:
    def test_parse(self) -> None:
        output = "Isotropic polarizability:  45.67890\nalpha_xx:  56.78900\n"
        r = PolarizabilityParser.parse(output)
        assert r["alpha_iso"] == pytest.approx(45.6789)
        assert r["alpha_xx"] == pytest.approx(56.789)

    def test_empty(self) -> None:
        assert PolarizabilityParser.parse("") == {}


class TestAromaticityParser:
    def test_parse(self) -> None:
        output = "NICS(0):  -8.12340\nHOMA:   0.98760\n"
        r = AromaticityParser.parse(output)
        assert r["NICS"] == pytest.approx(-8.1234)
        assert r["HOMA"] == pytest.approx(0.9876)

    def test_nics_scan(self) -> None:
        output = "NICS scan data\n  0.0000   -8.1234\n  1.0000  -10.5678\n"
        d = AromaticityParser.parse_nics_scan(output)
        assert len(d["distances"]) == 2

    def test_empty(self) -> None:
        assert AromaticityParser.parse("") == {}


class TestWavefunctionParser:
    def test_orbital_info(self) -> None:
        output = "   5   Alpha   Occ= 2.000000   E=  -0.72340 a.u.  -19.684 eV\n"
        o = WavefunctionParser.parse_orbital_info(output)
        assert o[0]["index"] == 5
        assert o[0]["spin"] == "alpha"

    def test_empty(self) -> None:
        assert WavefunctionParser.parse_orbital_info("") == []


class TestCubeParser:
    def test_parse(self) -> None:
        output = (
            "Grid dimensions: 80 x 80 x 80\n"
            "density.cube has been generated\n"
        )
        r = CubeParser.parse(output)
        assert r["grid_points"] == (80, 80, 80)
        assert "density.cube" in r["cube_files"]

    def test_empty(self) -> None:
        assert CubeParser.parse("") == {}


class TestUtilityParser:
    def test_geometry(self) -> None:
        output = "Bond length between atom  1(C ) and atom  2(N ):  1.3456\n"
        r = UtilityParser.parse_geometry(output)
        assert r["bond_lengths"][0]["length"] == pytest.approx(1.3456)

    def test_bla_boa(self) -> None:
        r = UtilityParser.parse_bla_boa("BLA=  0.05670\nBOA=  0.12340\n")
        assert r["BLA"] == pytest.approx(0.0567)

    def test_generated_files(self) -> None:
        output = "density.cube has been generated\nnew.wfn has been generated\n"
        assert len(UtilityParser.parse_generated_files(output)) == 2

    def test_empty_geometry(self) -> None:
        r = UtilityParser.parse_geometry("")
        assert r == {"bond_lengths": [], "angles": [], "dihedrals": []}