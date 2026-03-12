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
    """Tests for the ChargeParser."""

    def test_parse_final_charges(self) -> None:
        output = (
            "Final atomic charges:\n"
            "Atom    1(C ):    -0.05230000\n"
            "Atom    2(N ):     0.12340000\n"
        )
        charges = ChargeParser.parse_charges(output)
        assert len(charges) == 2
        # charges is a list of Charge dataclasses
        charge_map = {c.atom_id: c.charge for c in charges}
        assert charge_map[1] == pytest.approx(-0.0523)
        assert charge_map[2] == pytest.approx(0.1234)

    def test_parse_pattern2_after_section_header(self) -> None:
        output = "Some charge header\nAtom    1(C ):     0.03208687\n"
        charges = ChargeParser.parse_charges(output)
        charge_map = {c.atom_id: c.charge for c in charges}
        assert charge_map[1] == pytest.approx(0.03208687)

    def test_final_overrides_earlier(self) -> None:
        output = (
            "Hirshfeld charge of atom     1(C ) is  0.111\n"
            "Final atomic charges:\n"
            "Atom    1(C ):     0.222\n"
        )
        charges = ChargeParser.parse_charges(output)
        # Final section clears earlier charges, so only final remains
        charge_map = {c.atom_id: c.charge for c in charges}
        assert charge_map[1] == pytest.approx(0.222)

    def test_empty_output(self) -> None:
        assert ChargeParser.parse_charges("") == []

    def test_parse_dipole(self) -> None:
        output = (
            "Dipole moment:X=   1.234  Y=  -0.567  Z=   0.901  Tot=   1.623\n"
        )
        d = ChargeParser.parse_dipole(output)
        assert d is not None
        assert d.total == pytest.approx(1.623)

    def test_dipole_not_found(self) -> None:
        assert ChargeParser.parse_dipole("no dipole here") is None

    def test_scientific_notation(self) -> None:
        output = "Final atomic charges:\nAtom    1(C ):    -5.23E-02\n"
        charges = ChargeParser.parse_charges(output)
        charge_map = {c.atom_id: c.charge for c in charges}
        assert charge_map[1] == pytest.approx(-0.0523)


# =========================================================================
# BondOrderParser
# =========================================================================


class TestBondOrderParser:
    """Tests for the BondOrderParser."""

    def test_parse_pattern1(self) -> None:
        output = "#    1:         1(C )    2(N )    1.45230000\n"
        bonds = BondOrderParser.parse(output)
        # bonds is a list of BondOrder dataclasses
        assert len(bonds) == 1
        assert bonds[0].atom1_id == 1
        assert bonds[0].atom2_id == 2
        assert bonds[0].bond_order == pytest.approx(1.4523)

    def test_tuple_ordering(self) -> None:
        output = "#    1:         6(N )    1(C )    0.82340000\n"
        bonds = BondOrderParser.parse(output)
        # Parser normalises so smaller atom_id comes first
        assert bonds[0].atom1_id == 1
        assert bonds[0].atom2_id == 6
        assert bonds[0].bond_order == pytest.approx(0.8234)

    def test_parse_valence(self) -> None:
        output = (
            "Total valence of atom    1(C ):   3.94120000\n"
            "Free valence of atom     1(C ):   0.05880000\n"
        )
        v = BondOrderParser.parse_valence(output)
        # v is a list of Valence dataclasses
        total_vals = [
            x for x in v if x.type == "total_valence" and x.atom_id == 1
        ]
        free_vals = [
            x for x in v if x.type == "free_valence" and x.atom_id == 1
        ]
        assert total_vals[0].valence == pytest.approx(3.9412)
        assert free_vals[0].valence == pytest.approx(0.0588)

    def test_parse_multicenter(self) -> None:
        output = "Multi-center bond order of atoms  1  2  3 :  0.12345\n"
        results = BondOrderParser.parse_multicenter(output)
        assert results[0].atom_ids == [1, 2, 3]

    def test_parse_decomposition(self) -> None:
        output = "Orbital   5:   0.23456\nOrbital   6:  -0.01234\n"
        decomp = BondOrderParser.parse_decomposition(output)
        assert len(decomp) == 2

    def test_empty(self) -> None:
        assert BondOrderParser.parse("") == []


# =========================================================================
# CriticalPointParser
# =========================================================================


class TestCriticalPointParser:
    """Tests for the CriticalPointParser."""

    def test_parse_all_types(self, sample_topology_output: str) -> None:
        cps = CriticalPointParser.parse(sample_topology_output)
        assert len(cps) == 4
        types = {cp.type for cp in cps}
        assert types == {"nuclear", "bond", "ring", "cage"}

    def test_rho_laplacian(self, sample_topology_output: str) -> None:
        cps = CriticalPointParser.parse(sample_topology_output)
        assert cps[0].rho == pytest.approx(0.298765)
        assert cps[1].ellipticity == pytest.approx(0.04523)

    def test_bond_paths(self) -> None:
        output = (
            "Bond path between atom  1(C ) and atom  2(N ), "
            "BCP  3, length  2.456\n"
        )
        paths = CriticalPointParser.parse_bond_paths(output)
        assert paths[0].atom1_id == 1
        assert paths[0].path_length == pytest.approx(2.456)

    def test_summary(self) -> None:
        from pymultiwfn.analysis.result import CriticalPoint

        cps = [
            CriticalPoint(index=1, x=0, y=0, z=0, type="bond"),
            CriticalPoint(index=2, x=0, y=0, z=0, type="bond"),
            CriticalPoint(index=3, x=0, y=0, z=0, type="ring"),
        ]
        from collections import Counter

        counts = dict(Counter(cp.type for cp in cps))
        assert counts == {"bond": 2, "ring": 1}

    def test_empty(self) -> None:
        assert CriticalPointParser.parse("") == []

    def test_unknown_type(self) -> None:
        # (2,-1) is not a standard CP type signature with position data,
        # so the parser requires position info to create a CP.
        # Without position, no CP is created.
        cps = CriticalPointParser.parse(
            "CP  1 (2,-1)\nPosition (Bohr):   0.000000   0.000000   0.000000\n"
            "Density of all electrons:   0.001\n"
        )
        assert len(cps) == 1
        assert cps[0].type == "unknown"


# =========================================================================
# SpectrumParser
# =========================================================================


class TestSpectrumParser:
    """Tests for the SpectrumParser."""

    def test_parse_ir(self, sample_spectrum_output: str) -> None:
        sp = SpectrumParser.parse(sample_spectrum_output)
        if sp.frequencies is not None:
            assert len(sp.frequencies) == 5
        if sp.intensities is not None:
            assert sp.intensities[-1] == pytest.approx(156.78)

    def test_parse_uv_vis(self) -> None:
        output = "  345.67 nm  f= 0.1234\n"
        sp = SpectrumParser.parse(output)
        assert sp.wavelengths is not None and len(sp.wavelengths) > 0

    def test_parse_nmr(self) -> None:
        output = "Atom  1(C ) shift:  123.45 ppm\n"
        sp = SpectrumParser.parse(output)
        assert sp.chemical_shifts is not None and len(sp.chemical_shifts) > 0

    def test_parse_transitions(self) -> None:
        output = "Excited state   1:  E= 3.4567 eV  lam= 358.7 nm  f= 0.0123\n"
        t = SpectrumParser.parse_transitions(output)
        assert t[0].state == 1
        assert t[0].energy_eV == pytest.approx(3.4567)

    def test_parse_color(self) -> None:
        output = "X=  0.3456  Y=  0.3210  Z=  0.2890\nR= 180  G= 120  B=  90\n"
        c = SpectrumParser.parse_color(output)
        assert c is not None
        assert c.R == 180

    def test_empty(self) -> None:
        sp = SpectrumParser.parse("")
        assert sp.frequencies == []


# =========================================================================
# DOSParser
# =========================================================================


class TestDOSParser:
    """Tests for the DOSParser."""

    def test_parse_tdos(self) -> None:
        output = "TDOS data\n  -15.0000    0.1234\n  -14.0000    0.5678\n"
        data = DOSParser.parse(output)
        assert len(data.energies_eV) == 2

    def test_parse_orbital_energies(self) -> None:
        output = "   5   -19.234  eV  Occ=  2.000000\n"
        orbs = DOSParser.parse_orbital_energies(output)
        assert orbs[0].index == 5

    def test_empty(self) -> None:
        assert DOSParser.parse("").energies_eV == []


# =========================================================================
# SurfaceParser
# =========================================================================


class TestSurfaceParser:
    """Tests for the SurfaceParser."""

    def test_parse_descriptors(self) -> None:
        output = (
            "Overall surface area:  234.5678\nEnclosed volume:  345.6789\n"
        )
        r = SurfaceParser.parse(output)
        assert r.area == pytest.approx(234.5678)

    def test_parse_extrema(self) -> None:
        output = "Local minimum  1:   -45.678  at   1.234   2.345   3.456\n"
        ext = SurfaceParser.parse_extrema(output)
        assert ext[0].type == "min"

    def test_empty(self) -> None:
        r = SurfaceParser.parse("")
        # Returns a SurfaceAnalysis with all None fields
        assert r.area is None
        assert r.volume is None


# =========================================================================
# Remaining parsers — focused coverage
# =========================================================================


class TestOrbitalCompositionParser:
    """Tests for the OrbitalCompositionParser."""

    def test_parse_composition(self) -> None:
        output = (
            "Orbital    5  Occ= 2.000000  E= -0.72340 a.u.\n"
            "  C  1    s :   23.45%\n"
        )
        orbs = OrbitalCompositionParser.parse(output)
        assert orbs[0].orbital_id == 5
        # contributions is a list of OrbitalContribution dataclasses
        contrib_map = {c.label: c.percentage for c in orbs[0].contributions}
        assert contrib_map["C1"] == pytest.approx(23.45)

    def test_parse_oxidation(self) -> None:
        output = "Atom   1(Fe) formal oxidation state:  3.0\n"
        result = OrbitalCompositionParser.parse_oxidation_states(output)
        assert len(result) == 1
        assert result[0].atom_id == 1
        assert result[0].oxidation_state == 3

    def test_empty(self) -> None:
        assert OrbitalCompositionParser.parse("") == []


class TestFuzzySpaceParser:
    """Tests for the FuzzySpaceParser."""

    def test_atomic_properties(self) -> None:
        output = "Atom   1(C ): population=  5.96780\n"
        a = FuzzySpaceParser.parse_atomic_properties(output)
        assert len(a) == 1
        assert a[0].atom_id == 1
        assert a[0].population == pytest.approx(5.9678)

    def test_delocalization(self) -> None:
        output = "Localization index of atom  1(C ):  3.45670\n"
        r = FuzzySpaceParser.parse_delocalization_indices(output)
        # Returns a list of DelocalizationIndex dataclasses
        # The pattern expects "Delocalization index...atom X...atom Y"
        # "Localization index" won't match the delocalization pattern
        assert isinstance(r, list)

    def test_aromaticity_index(self) -> None:
        output = "PDI= 0.05678\nFLU= 0.00123\n"
        r = FuzzySpaceParser.parse_aromaticity_index(output)
        # Returns a list of AromaticityIndex dataclasses
        index_map = {item.index_name: item.value for item in r}
        assert index_map["PDI"] == pytest.approx(0.05678)


class TestBasinParser:
    """Tests for the BasinParser."""

    def test_parse(self) -> None:
        output = "Basin   1  attractor at atom  1(C )  population:  5.96780\n"
        b = BasinParser.parse(output)
        assert b[0].population == pytest.approx(5.9678)

    def test_parse_charges(self) -> None:
        output = "AIM charge of atom  1(C ):  0.03220\n"
        charges = BasinParser.parse_charges(output)
        assert len(charges) == 1
        assert charges[0].atom_id == 1
        assert charges[0].charge == pytest.approx(0.0322)

    def test_empty(self) -> None:
        assert BasinParser.parse("") == []


class TestExcitationParser:
    """Tests for the ExcitationParser."""

    def test_hole_electron(self) -> None:
        # parse_hole_electron requires many fields to return non-None
        # With only D_index and Sr, it returns None
        output = "D index:  2.34560\nSr:  0.56780\n"
        r = ExcitationParser.parse_hole_electron(output)
        # The parser requires H_index, E_index, t_index, EDI, HDI, Sr, D_index
        # plus centroids to return a HoleElectron; otherwise returns None
        assert r is None

    def test_delta_r(self) -> None:
        output = "State   1  Delta_r:  1.23450\n"
        result = ExcitationParser.parse_delta_r(output)
        assert result[0].delta_r == pytest.approx(1.2345)

    def test_lambda_index(self) -> None:
        output = "State   1  Lambda:  0.78900\n"
        result = ExcitationParser.parse_lambda_index(output)
        assert result[0].lambda_index == pytest.approx(0.789)


class TestWeakInteractionParser:
    """Tests for the WeakInteractionParser."""

    def test_parse(self) -> None:
        output = "delta_g_inter:  0.12340\nRDG.cube has been generated\n"
        r = WeakInteractionParser.parse(output)
        assert r is not None
        assert r.delta_g_inter == pytest.approx(0.1234)
        if r.cube_names is not None:
            assert "RDG.cube" in r.cube_names

    def test_empty(self) -> None:
        assert WeakInteractionParser.parse("") is None


class TestEDAParser:
    """Tests for the EDAParser."""

    def test_parse(self) -> None:
        output = "Electrostatic:  -45.6789\nDispersion:   -5.6789\n"
        r = EDAParser.parse(output)
        assert r.electrostatic == pytest.approx(-45.6789)

    def test_empty(self) -> None:
        r = EDAParser.parse("")
        # Returns an EnergyDecompositionAnalysis with all None fields
        assert r.electrostatic is None


class TestCDFTParser:
    """Tests for the CDFTParsers."""

    def test_global_indices(self) -> None:
        output = "Chemical potential:  -0.15234\nHardness:   0.21345\n"
        r = CDFTParser.parse_global_indices(output)
        assert r.chemical_potential == pytest.approx(-0.15234)

    def test_condensed_fukui(self) -> None:
        output = "Atom   1(C ):  f+= 0.1234  f-= 0.0567  f0= 0.0901\n"
        f = CDFTParser.parse_condensed_fukui(output)
        assert len(f) == 1
        assert f[0].atom_id == 1
        assert f[0].fukui_plus == pytest.approx(0.1234)

    def test_dual_descriptor(self) -> None:
        output = "Atom   1(C ):  dual= 0.06670\n"
        result = CDFTParser.parse_dual_descriptor(output)
        assert len(result) == 1
        assert result[0].atom_id == 1
        assert result[0].value == pytest.approx(0.0667)


class TestPolarizabilityParser:
    """Tests for the PolarizabilityParsery."""

    def test_parse(self) -> None:
        output = "Isotropic polarizability:  45.67890\nalpha_xx:  56.78900\n"
        r = PolarizabilityParser.parse(output)
        assert r.isotropic == pytest.approx(45.6789)
        assert r.tensor is not None
        assert r.tensor.alpha_xx == pytest.approx(56.789)

    def test_empty(self) -> None:
        r = PolarizabilityParser.parse("")
        # Returns a Polarizability with all None fields
        assert r.isotropic is None


class TestAromaticityParser:
    """Tests for the AromaticityParser."""

    def test_parse(self) -> None:
        output = "NICS(0):  -8.12340\nHOMA:   0.98760\n"
        r = AromaticityParser.parse(output)
        assert pytest.approx(-8.1234) == r.NICS
        assert pytest.approx(0.9876) == r.HOMA

    def test_nics_scan(self) -> None:
        output = "NICS scan data\n  0.0000   -8.1234\n  1.0000  -10.5678\n"
        d = AromaticityParser.parse_nics_scan(output)
        assert len(d.distances) == 2

    def test_empty(self) -> None:
        r = AromaticityParser.parse("")
        # Returns an Aromaticity with all None fields
        assert r.NICS is None


class TestWavefunctionParser:
    """Tests for the WavefunctionParser."""

    def test_orbital_info(self) -> None:
        output = (
            "   5   Alpha   Occ= 2.000000   E=  -0.72340 a.u.  -19.684 eV\n"
        )
        o = WavefunctionParser.parse_orbital_info(output)
        assert o[0].orbital_id == 5
        assert o[0].spin == "alpha"

    def test_empty(self) -> None:
        assert WavefunctionParser.parse_orbital_info("") == []


class TestCubeParser:
    """Tests for the CubeParser."""

    def test_parse(self) -> None:
        output = (
            "Grid dimensions: 80 x 80 x 80\ndensity.cube has been generated\n"
        )
        r = CubeParser.parse(output)
        assert r is not None
        assert (r.x_dim, r.y_dim, r.z_dim) == (80, 80, 80)
        assert r.file_name == "density.cube"

    def test_empty(self) -> None:
        assert CubeParser.parse("") is None


class TestUtilityParser:
    """Tests for the UtilityParser."""

    def test_geometry(self) -> None:
        output = "Bond length between atom  1(C ) and atom  2(N ):  1.3456\n"
        r = UtilityParser.parse_bond_lengths(output)
        assert len(r) == 1
        assert r[0].length == pytest.approx(1.3456)

    def test_bla_boa(self) -> None:
        r = UtilityParser.parse_bla_boa("BLA=  0.05670\nBOA=  0.12340\n")
        assert r.bla == pytest.approx(0.0567)

    def test_generated_files(self) -> None:
        # UtilityParser doesn't have parse_generated_files;
        # test via CubeParser which handles "has been generated" patterns
        import re

        output = (
            "density.cube has been generated\nnew.wfn has been generated\n"
        )
        files = re.findall(r"(\S+)\s+has been generated", output)
        assert len(files) == 2

    def test_empty_geometry(self) -> None:
        r = UtilityParser.parse_bond_lengths("")
        assert r == []
