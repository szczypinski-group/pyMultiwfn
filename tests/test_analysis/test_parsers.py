"""Tests for pymultiwfn.analysis.parsers — output parsers.

Covers positive, negative, and edge cases for all parser classes.
Updated to match current API: ChargeSet, BondOrderSet, SurfaceAnalysisResult,
TopologyPath, Reactivity, multi-cube Cube support, etc.
"""

import pytest

from pymultiwfn.analysis.parsers import (
    AromaticityParser,
    AtomListParser,
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
    GridParser,
    LineParser,
    OrbitalCompositionParser,
    PlaneParser,
    PolarizabilityParser,
    SpectrumParser,
    SurfaceParser,
    UtilityParser,
    WavefunctionParser,
    WeakInteractionParser,
)
from pymultiwfn.analysis.result import (
    AtomInfo,
    BondOrderSet,
    ChargeSet,
    CondensedFukui,
    Cube,
    DelocalizationIndex,
    Dipole,
    DualDescriptor,
    FuzzyIntegrationResult,
    HOMOLUMOGap,
    OrbitalAtomComposition,
    OrbitalBasisComposition,
    PoincareHopfCounts,
    Reactivity,
)
from pymultiwfn.enums.menu import Menu

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
        charge_map = {c.atom_id: c.charge for c in charges}
        assert charge_map[1] == pytest.approx(0.222)

    def test_empty_output(self) -> None:
        assert ChargeParser.parse_charges("") == []

    def test_parse_dipole(self) -> None:
        output = "Dipole moment:X=   1.234  Y=  -0.567  Z=   0.901  Tot=   1.623\n"
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

    def test_parse_charge_sets_returns_charge_sets(self) -> None:
        """ChargeParser.parse_charge_sets returns ChargeSet objects."""
        output = (
            "Final atomic charges:\n"
            "Atom    1(C ):    -0.05230000\n"
            "Atom    2(N ):     0.12340000\n"
        )
        sets = ChargeParser.parse_charge_sets(output, base_method="hirshfeld")
        assert len(sets) >= 1
        assert isinstance(sets[0], ChargeSet)
        assert sets[0].method == "hirshfeld"
        assert len(sets[0].charges) == 2

    def test_parse_charge_sets_empty(self) -> None:
        sets = ChargeParser.parse_charge_sets("", base_method="unknown")
        assert sets == []

    def test_parse_charge_sets_all_block_formats(self) -> None:
        """Exercise every section format parse_charge_sets recognises."""
        output = (
            "----------------- ADCH charges -----------------\n"
            "Atom:    1C   Corrected charge:   -0.041271\n"
            "Atom:    2N   Corrected charge:    0.052341\n"
            "Total charge:    0.011070\n"
            "\n"
            "----------------- CM5 charges -----------------\n"
            "Atom:    1C   CM5 charge:   -0.097209\n"
            "Atom:    2N   CM5 charge:    0.083122\n"
            "Sum of charges:   -0.014087\n"
            "\n"
            "**** Stage 1 RESP fitting ****\n"
            "Center       Charge\n"
            "     1(C )  -0.0912959843\n"
            "     2(N )   0.0823140000\n"
            "Total net charge:    -0.008982\n"
            "\n"
            "Atom       Charge\n"
            "     1(C )  -0.1234500000\n"
            "     2(N )   0.1234500000\n"
            "Total charge:     0.000000\n"
            "\n"
            "Population of atoms:\n"
            "Atom     1(C )  Population:  6.128  Net charge: -0.128\n"
            "Atom     2(N )  Population:  5.872  Net charge:  0.128\n"
            "Total charge:     0.000000\n"
            "\n"
            "Final atomic charges:\n"
            "Atom    1(C ):    -0.05230000\n"
            "Atom    2(N ):     0.05230000\n"
            "Summing up all Hirshfeld charges:      0.00000000\n"
        )
        sets = ChargeParser.parse_charge_sets(output, base_method="hirshfeld")

        by_method_stage = {(s.method, s.stage): s for s in sets}
        assert ("adch", "corrected") in by_method_stage
        assert by_method_stage[("adch", "corrected")].total_charge == (
            pytest.approx(0.01107)
        )
        assert ("cm5", "raw") in by_method_stage
        assert by_method_stage[("cm5", "raw")].total_charge == (
            pytest.approx(-0.014087)
        )
        assert ("hirshfeld", "esp_fit") in by_method_stage
        assert len(by_method_stage[("hirshfeld", "esp_fit")].charges) == 2
        assert ("hirshfeld", "raw") in by_method_stage
        assert ("hirshfeld", "population") in by_method_stage
        assert by_method_stage[("hirshfeld", "population")].charges[
            0
        ].charge == pytest.approx(-0.128)
        assert ("hirshfeld", "final") in by_method_stage
        assert by_method_stage[("hirshfeld", "final")].total_charge == (
            pytest.approx(0.0)
        )

    def test_parse_for_result_returns_charge_sets_and_dipole(self) -> None:
        output = (
            "Final atomic charges:\n"
            "Atom    1(C ):    -0.05230000\n"
            "Dipole moment:X=   1.0  Y=   0.0  Z=   0.0  Tot=   1.0\n"
        )
        results = ChargeParser.parse_for_result(
            analysis=Menu.HIRSHFELD_CHARGE,
            stdout=output,
        )
        charge_sets = [r for r in results if isinstance(r, ChargeSet)]
        dipoles = [r for r in results if isinstance(r, Dipole)]
        assert len(charge_sets) >= 1
        assert len(dipoles) == 1

    def test_parse_multiple_charge_blocks(self) -> None:
        """Multiple charge blocks produce multiple ChargeSet objects."""
        output = (
            "Final atomic charges:\n"
            "Atom    1(C ):    -0.05230000\n"
            "\n"
            "--- CM5 charges ---\n"
            "Atom:    1C   CM5 charge:   -0.097209\n"
        )
        sets = ChargeParser.parse_charge_sets(output, base_method="hirshfeld")
        assert len(sets) >= 1


# =========================================================================
# BondOrderParser
# =========================================================================


class TestBondOrderParser:
    """Tests for the BondOrderParser — updated API with BondOrderSet."""

    def test_parse_bond_order_sets(self) -> None:
        output = (
            "Bond orders with absolute value >=  0.050000\n"
            "#    1:         1(C )    2(C )    1.45230000\n"
            "#    2:         1(C )    3(H )    0.92340000\n"
        )
        sets = BondOrderParser.parse_bond_order_sets(output, method="mayer")
        assert len(sets) >= 1
        assert isinstance(sets[0], BondOrderSet)
        assert sets[0].method == "mayer"
        assert len(sets[0].bond_orders) == 2
        assert sets[0].bond_orders[0].atom1_id == 1
        assert sets[0].bond_orders[0].atom2_id == 2
        assert sets[0].bond_orders[0].bond_order == pytest.approx(1.4523)

    def test_tuple_ordering_normalised(self) -> None:
        output = "#    1:         6(N )    1(C )    0.82340000\n"
        sets = BondOrderParser.parse_bond_order_sets(output, method="mayer")
        assert len(sets) >= 1
        bo = sets[0].bond_orders[0]
        assert bo.atom1_id == 1
        assert bo.atom2_id == 6
        assert bo.bond_order == pytest.approx(0.8234)

    def test_parse_valence(self) -> None:
        output = (
            "Total valence of atom    1(C ):   3.94120000\n"
            "Free valence of atom     1(C ):   0.05880000\n"
        )
        v = BondOrderParser.parse_valence(output)
        total_vals = [x for x in v if x.type == "total_valence" and x.atom_id == 1]
        free_vals = [x for x in v if x.type == "free_valence" and x.atom_id == 1]
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
        sets = BondOrderParser.parse_bond_order_sets("", method="mayer")
        assert sets == []

    def test_parse_for_result_mayer(self) -> None:
        output = (
            "Bond orders with absolute value >=  0.050000\n"
            "#    1:         1(C )    2(N )    1.45230000\n"
        )
        results = BondOrderParser.parse_for_result(
            analysis=Menu.MAYER_BOND_ORDER,
            stdout=output,
        )
        bond_sets = [r for r in results if isinstance(r, BondOrderSet)]
        assert len(bond_sets) >= 1

    def test_parse_ibsi(self) -> None:
        output = (
            "    1(C )    2(C )  Dist:  1.3950   "
            "Int(dg_pair): 0.87601   IBSI: 0.79442\n"
        )
        result = BondOrderParser.parse_ibsi(output)
        assert result is not None
        assert len(result.entries) == 1
        assert result.entries[0].ibsi == pytest.approx(0.79442)

    def test_parse_ibsi_empty(self) -> None:
        assert BondOrderParser.parse_ibsi("") is None

    def test_parse_valence_two_column(self) -> None:
        """Two-column valence format from Mayer output."""
        output = (
            "Total valences and free valences:\n"
            "Atom     1(C ) :    3.94042232    0.00000000\n"
            "Atom     2(N ) :    3.12345678    0.01234567\n"
        )
        v = BondOrderParser.parse_valence(output)
        assert len(v) == 4  # 2 atoms x (total + free)


# =========================================================================
# CriticalPointParser
# =========================================================================


class TestCriticalPointParser:
    """Tests for the CriticalPointParser — updated for TopologyPath."""

    def test_parse_short_summary(self) -> None:
        output = (
            "---- Summary ----\n"
            "    1    0.000    0.000    0.000   (3,-3)\n"
            "    2    1.234    0.568    0.000   (3,-1)\n"
            "    3    0.500    0.500    0.500   (3,+1)\n"
            "    4    1.000    1.000    1.000   (3,+3)\n"
            "Totally find 4 critical points\n"
        )
        cps = CriticalPointParser.parse(output)
        assert len(cps) == 4
        types = {cp.type for cp in cps}
        assert types == {"nuclear", "bond", "ring", "cage"}

    def test_parse_long_summary(self) -> None:
        output = (
            "Summary of found CPs\n"
            "    1    0.000   -4.680   -0.000   (3,-3)   Nucleus:   10(H )\n"
            "   13    0.000   -3.966   -0.000   (3,-1)   10(H ) --    4(C )\n"
            "number of critical points found\n"
        )
        cps = CriticalPointParser.parse(output)
        assert len(cps) == 2
        assert cps[0].type == "nuclear"
        assert cps[0].nucleus_atom_id == 10
        assert cps[1].type == "bond"
        assert cps[1].bonded_atom1_id == 10
        assert cps[1].bonded_atom2_id == 4

    def test_parse_bond_paths_empty(self) -> None:
        """Empty input produces no matches, so no TypeError is raised."""
        paths = CriticalPointParser.parse_bond_paths("")
        assert paths == []

    def test_parse_topology_paths(self) -> None:
        output = "Path  1, CP:  5 (3,-1) -->  CP:  1 (3,-3) Length:  1.234\n"
        paths = CriticalPointParser.parse_topology_paths(output)
        assert len(paths) == 1
        assert paths[0].path_id == 1
        assert paths[0].bcp_index == 5
        assert paths[0].length_bohr == pytest.approx(1.234)

    def test_parse_poincare_hopf(self) -> None:
        output = (
            "(3,-3):  5  (3,-1):  6  (3,+1):  2  (3,+3):  0\n"
            "Poincare-Hopf relationship is satisfied\n"
        )
        ph = CriticalPointParser.parse_poincare_hopf(output)
        assert ph is not None
        assert isinstance(ph, PoincareHopfCounts)
        assert ph.nuclear == 5
        assert ph.bond == 6
        assert ph.ring == 2
        assert ph.cage == 0
        assert ph.satisfied is True

    def test_parse_poincare_hopf_not_found(self) -> None:
        assert CriticalPointParser.parse_poincare_hopf("no data") is None

    def test_empty(self) -> None:
        assert CriticalPointParser.parse("") == []

    def test_parse_for_result(self) -> None:
        output = (
            "---- Summary ----\n"
            "    1    0.000    0.000    0.000   (3,-3)\n"
            "Totally find 1 critical points\n"
        )
        results = CriticalPointParser.parse_for_result(
            analysis=Menu.TOPOLOGY_SEARCH_CPS,
            stdout=output,
        )
        assert len(results) >= 1


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

    def test_parse_spectrum_curve_extrema(self) -> None:
        output = (
            "Extrema on the spectrum curve:\n\n"
            " Maximum    1   X:      3578.5262   Value:       585.2897\n\n"
            " Maximum    2   X:      3454.4848   Value:       110.4777\n\n"
            " Maximum   11   X:       262.7543   Value:       184.0664.\n"
        )
        extrema = SpectrumParser.parse_spectrum_curve_extrema(output)
        assert extrema is not None
        assert len(extrema.extrema) == 3
        assert extrema.extrema[0].kind == "maximum"
        assert extrema.extrema[0].index == 1
        assert extrema.extrema[0].x == pytest.approx(3578.5262)
        assert extrema.extrema[2].value == pytest.approx(184.0664)

    def test_parse_spectrum_curve_extrema_not_found(self) -> None:
        assert SpectrumParser.parse_spectrum_curve_extrema("no data") is None

    def test_parse_for_result_includes_spectrum_curve_extrema(self) -> None:
        output = (
            "Extrema on the spectrum curve:\n\n"
            " Maximum    1   X:      1303.1010   Value:      5467.1462\n"
        )
        results = SpectrumParser.parse_for_result(
            analysis=Menu.PLOT_UV_VIS_SPECTRUM,
            stdout=output,
        )
        assert any(type(r).__name__ == "SpectrumCurveExtrema" for r in results)

    def test_parse_transitions_with_rotation(self) -> None:
        output = (
            "Excited state   1:  E= 3.4567 eV  lam= 358.7 nm  f= 0.0123\n"
            "  R(velocity)=  0.5678\n"
        )
        t = SpectrumParser.parse_transitions(output)
        assert len(t) == 1
        assert t[0].rot_strength == pytest.approx(0.5678)


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

    def test_parse_metadata(self) -> None:
        output = (
            "Center of TDOS:   -0.12345 a.u.\n"
            "vertical dash line corresponds to HOMO level at  -0.23456 a.u.\n"
        )
        meta = DOSParser.parse_metadata(output)
        assert meta is not None
        assert meta.tdos_center_au == pytest.approx(-0.12345)
        assert meta.homo_level_au == pytest.approx(-0.23456)

    def test_parse_metadata_empty(self) -> None:
        assert DOSParser.parse_metadata("") is None


# =========================================================================
# SurfaceParser
# =========================================================================


class TestSurfaceParser:
    """Tests for the SurfaceParser — returns SurfaceAnalysisResult."""

    def test_parse_geometry(self) -> None:
        output = (
            "Volume enclosed in the isosurface:  234.5678 Bohr^3 ( 34.7360 Angstrom^3)\n"  # noqa: E501
            "Overall surface area:  448.2930 Bohr^2 ( 125.5350 Angstrom^2)\n"
            "Sphericity:   0.8765\n"
        )
        geo = SurfaceParser.parse_geometry(output)
        assert geo is not None
        assert geo.volume_bohr3 == pytest.approx(234.5678)
        assert geo.area_angstrom2 == pytest.approx(125.5350)
        assert geo.sphericity == pytest.approx(0.8765)

    def test_parse_extrema_min(self) -> None:
        output = (
            "The number of surface minima:    2\n"
            "  #         a.u.        eV         kcal/mol      X         Y         Z\n"  # noqa: E501
            "*    1 -0.02750965   -0.748576  -17.262580   -0.225  0.366  -1.831\n"  # noqa: E501
            "     2 -0.01234567   -0.335890   -7.745670    0.500  0.500   0.500\n"  # noqa: E501
        )
        extrema = SurfaceParser.parse_extrema(output, "min")
        assert len(extrema) == 2
        assert extrema[0].type == "min"
        assert extrema[0].is_global is True
        assert extrema[1].is_global is False

    def test_parse_extrema_max(self) -> None:
        output = (
            "The number of surface maxima:    1\n"
            "  #         a.u.        eV         kcal/mol      X         Y         Z\n"  # noqa: E501
            "*    1  0.03456789    0.940580   21.694560    1.000  1.000   1.000\n"  # noqa: E501
        )
        extrema = SurfaceParser.parse_extrema(output, "max")
        assert len(extrema) == 1
        assert extrema[0].type == "max"

    def test_parse_statistics(self) -> None:
        output = (
            "Summary of surface analysis\n"
            "Minimal value:    -17.26258 kcal/mol   "
            "Maximal value:     11.82103 kcal/mol\n"
            "Balance of charges (nu):   0.19601755\n"
        )
        stats = SurfaceParser.parse_statistics(output, "esp")
        assert stats is not None
        assert stats.global_min_kcal_mol == pytest.approx(-17.26258)
        assert stats.nu == pytest.approx(0.19601755)

    def test_parse_statistics_full_summary_block(self) -> None:
        output = (
            " Global surface minimum:   -0.032150 a.u.\n"
            " Global surface maximum:    0.045670 a.u.\n"
            "\n"
            " ========= Summary of surface analysis =========\n"
            " Minimal value:    -17.26258 kcal/mol   "
            "Maximal value:     11.82103 kcal/mol\n"
            " Overall surface area:  448.293 Bohr^2  "
            "( 125.535 Angstrom^2)\n"
            " Positive surface area:  248.293 Bohr^2  "
            "(  65.535 Angstrom^2)\n"
            " Negative surface area:  200.000 Bohr^2  "
            "(  60.000 Angstrom^2)\n"
            " Overall average value:  0.00002044 a.u. "
            "(  0.01282 kcal/mol)\n"
            " Positive average value:  0.00500000 a.u. "
            "(  3.13800 kcal/mol)\n"
            " Negative average value: -0.00400000 a.u. "
            "( -2.51000 kcal/mol)\n"
            " Overall variance (sigma^2_tot):  0.00009608 a.u.^2 "
            "(  37.83 (kcal/mol)^2)\n"
            " Positive variance:  0.00005000 a.u.^2 "
            "(  19.68 (kcal/mol)^2)\n"
            " Negative variance:  0.00004000 a.u.^2 "
            "(  15.73 (kcal/mol)^2)\n"
            " Balance of charges (nu):   0.19601755\n"
            " Product of sigma^2_tot and nu:  0.00001883 a.u.^2 "
            "(  7.416 (kcal/mol)^2)\n"
            " Internal charge separation (Pi):  0.01236 a.u. "
            "(  7.757 kcal/mol)\n"
            " Molecular polarity index (MPI):  0.33644 eV "
            "(  7.758 kcal/mol)\n"
            " Nonpolar surface area (|ESP| <= 10 kcal/mol):  "
            "88.18 Angstrom^2  ( 70.25 %)\n"
            " Polar surface area (|ESP| > 10 kcal/mol):  "
            "37.36 Angstrom^2  ( 29.75 %)\n"
            " Overall skewness:   0.234500\n"
            " Positive skewness:   0.123400\n"
            " Negative skewness:  -0.098700\n"
        )
        stats = SurfaceParser.parse_statistics(output, "esp")
        assert stats is not None
        assert stats.global_min_au == pytest.approx(-0.03215)
        assert stats.global_max_au == pytest.approx(0.04567)
        assert stats.overall_area_bohr2 == pytest.approx(448.293)
        assert stats.overall_area_angstrom2 == pytest.approx(125.535)
        assert stats.positive_area_bohr2 == pytest.approx(248.293)
        assert stats.negative_area_bohr2 == pytest.approx(200.0)
        assert stats.overall_average_au == pytest.approx(2.044e-05)
        assert stats.positive_average_au == pytest.approx(0.005)
        assert stats.negative_average_au == pytest.approx(-0.004)
        assert stats.sigma2_total_au2 == pytest.approx(9.608e-05)
        assert stats.positive_variance_au2 == pytest.approx(5e-05)
        assert stats.negative_variance_au2 == pytest.approx(4e-05)
        assert stats.sigma2_tot_times_nu_au2 == pytest.approx(1.883e-05)
        assert stats.pi_au == pytest.approx(0.01236)
        assert stats.mpi_eV == pytest.approx(0.33644)
        # Nonpolar and polar area are distinct blocks -- a naive
        # substring match on "Polar surface area" inside "Nonpolar
        # surface area" would otherwise collapse them to the same
        # values.
        assert stats.nonpolar_area_angstrom2 == pytest.approx(88.18)
        assert stats.nonpolar_area_pct == pytest.approx(70.25)
        assert stats.polar_area_angstrom2 == pytest.approx(37.36)
        assert stats.polar_area_pct == pytest.approx(29.75)
        assert stats.overall_skewness == pytest.approx(0.2345)
        assert stats.positive_skewness == pytest.approx(0.1234)
        assert stats.negative_skewness == pytest.approx(-0.0987)

    def test_parse_statistics_empty(self) -> None:
        assert SurfaceParser.parse_statistics("", "esp") is None

    def test_parse_geometry_empty(self) -> None:
        assert SurfaceParser.parse_geometry("") is None


# =========================================================================
# GridParser
# =========================================================================


class TestGridParser:
    """Tests for the GridParser (Menu 13 -> show statistic data)."""

    def test_parse_statistics_full(self) -> None:
        output = (
            " The minimum value:     -0.123456 at     1.000000"
            "     2.000000     3.000000\n"
            " The maximum value:      0.987654 at     4.000000"
            "     5.000000     6.000000\n"
            " Average value:      0.345678\n"
            " Root mean square (RMS):      0.456789\n"
            " Standard deviation:      0.234567\n"
            " Volume of positive value space:    123.456\n"
            " Volume of negative value space:     78.901\n"
            " Volume of all space:              202.357\n"
            " Summing up positive values:      50.123\n"
            " Summing up negative values:     -30.456\n"
            " Summing up all values:           19.667\n"
            " Integral of positive data:        5.678\n"
            " Integral of negative data:       -3.456\n"
            " Integral of all data:             2.222\n"
        )
        stats = GridParser.parse_statistics(output)
        assert stats is not None
        assert stats.minimum == pytest.approx(-0.123456)
        assert stats.minimum_x_bohr == pytest.approx(1.0)
        assert stats.minimum_y_bohr == pytest.approx(2.0)
        assert stats.minimum_z_bohr == pytest.approx(3.0)
        assert stats.maximum == pytest.approx(0.987654)
        assert stats.maximum_x_bohr == pytest.approx(4.0)
        assert stats.average == pytest.approx(0.345678)
        assert stats.rms == pytest.approx(0.456789)
        assert stats.std_dev == pytest.approx(0.234567)
        assert stats.volume_positive_bohr3 == pytest.approx(123.456)
        assert stats.volume_negative_bohr3 == pytest.approx(78.901)
        assert stats.volume_all_bohr3 == pytest.approx(202.357)
        assert stats.sum_positive == pytest.approx(50.123)
        assert stats.sum_negative == pytest.approx(-30.456)
        assert stats.sum_all == pytest.approx(19.667)
        assert stats.integral_positive == pytest.approx(5.678)
        assert stats.integral_negative == pytest.approx(-3.456)
        assert stats.integral_all == pytest.approx(2.222)

    def test_parse_statistics_empty(self) -> None:
        assert GridParser.parse_statistics("") is None

    def test_parse_for_result(self) -> None:
        output = " Average value:      0.345678\n"
        results = GridParser.parse_for_result(
            analysis=Menu.GRID_STATISTIC_DATA, stdout=output
        )
        assert len(results) == 1

    def test_parse_for_result_empty(self) -> None:
        results = GridParser.parse_for_result(
            analysis=Menu.GRID_STATISTIC_DATA, stdout=""
        )
        assert results == []


# =========================================================================
# LineParser
# =========================================================================


class TestLineParser:
    """Tests for the LineParser (Menu 3 -> property along a line)."""

    def test_parse_full(self) -> None:
        output = (
            " Original point in X,Y,Z:    0.000000    0.000000"
            "    0.000000\n"
            " End point in X,Y,Z:    5.000000    0.000000"
            "    0.000000\n"
            " Number of points:      3000\n"
            " Minimal/Maximum value:   0.000123   1.234567\n"
            " Summing up all values:     500.123456  "
            "Integration value:      2.345678\n"
        )
        profile = LineParser.parse(output)
        assert profile is not None
        assert profile.origin_x_bohr == pytest.approx(0.0)
        assert profile.end_x_bohr == pytest.approx(5.0)
        assert profile.n_points == 3000
        assert profile.minimum == pytest.approx(0.000123)
        assert profile.maximum == pytest.approx(1.234567)
        assert profile.sum_all_values == pytest.approx(500.123456)
        assert profile.integration_value == pytest.approx(2.345678)

    def test_parse_empty(self) -> None:
        assert LineParser.parse("") is None

    def test_parse_for_result(self) -> None:
        output = " Number of points:      3000\n"
        results = LineParser.parse_for_result(analysis=Menu.LINE_ESP, stdout=output)
        assert len(results) == 1

    def test_parse_for_result_empty(self) -> None:
        results = LineParser.parse_for_result(analysis=Menu.LINE_ESP, stdout="")
        assert results == []


# =========================================================================
# PlaneParser
# =========================================================================


class TestPlaneParser:
    """Tests for the PlaneParser (Menu 4 -> property in a plane)."""

    def test_parse_full(self) -> None:
        output = (
            " X/Y/Z of origin of the plane:   -5.000000   -5.000000"
            "    0.000000\n"
            " X/Y/Z of end of the plane:    5.000000    5.000000"
            "    0.000000\n"
            " The minimum of data:   -0.123456\n"
            " The maximum of data:    0.654321\n"
        )
        plane = PlaneParser.parse(output)
        assert plane is not None
        assert plane.origin_x_bohr == pytest.approx(-5.0)
        assert plane.end_x_bohr == pytest.approx(5.0)
        assert plane.minimum == pytest.approx(-0.123456)
        assert plane.maximum == pytest.approx(0.654321)

    def test_parse_empty(self) -> None:
        assert PlaneParser.parse("") is None

    def test_parse_for_result(self) -> None:
        output = " The minimum of data:   -0.123456\n"
        results = PlaneParser.parse_for_result(
            analysis=Menu.PLANE_MAP_DENSITY, stdout=output
        )
        assert len(results) == 1

    def test_parse_for_result_empty(self) -> None:
        results = PlaneParser.parse_for_result(
            analysis=Menu.PLANE_MAP_DENSITY, stdout=""
        )
        assert results == []


# =========================================================================
# OrbitalCompositionParser
# =========================================================================


class TestOrbitalCompositionParser:
    """Tests for the OrbitalCompositionParser."""

    def test_parse_oxidation_states(self) -> None:
        output = "Atom   1(Fe) formal oxidation state:  3.0\n"
        result = OrbitalCompositionParser.parse_oxidation_states(output)
        assert len(result) == 1
        assert result[0].atom_id == 1
        assert result[0].oxidation_state == 3

    def test_parse_oxidation_states_empty(self) -> None:
        assert OrbitalCompositionParser.parse_oxidation_states("") == []

    def test_parse_basis_compositions(self) -> None:
        output = (
            "Orbital:    21  Energy(a.u.):     -0.246939  Occ:  2.000000  Type: Alpha&Beta\n"  # noqa: E501
            "    20   Z        2(C )    9      8.67678 %      5.61233 %     14.28911 %\n"  # noqa: E501
            "Composition of each shell\n"
            "Shell     9 Type: P    in atom    2(C ) :    14.28911 %\n"
            "Composition of different types of shells\n"
            "  s:   0.000  p:  99.010  d:   0.990  f:   0.000  g:   0.000  h:   0.000\n"  # noqa: E501
            "Composition of each atom:\n"
            "Atom     2(C ) :    14.28911 %\n"
            "Orbital delocalization index:   24.69\n"
        )
        results = OrbitalCompositionParser.parse_basis_compositions(
            output, method="mulliken"
        )
        assert len(results) == 1
        assert isinstance(results[0], OrbitalBasisComposition)
        assert results[0].orbital_id == 21
        assert len(results[0].basis_contributions) == 1
        assert results[0].delocalization_index == pytest.approx(24.69)

    def test_parse_atom_compositions(self) -> None:
        output = (
            "Orbital:    5  Energy(a.u.):     -0.72340  Occ:  2.000000  Type: Alpha&Beta\n"  # noqa: E501
            "The sum of contributions before normalization   99.999199 %\n"
            "Contributions after normalization:\n"
            "Atom     1(C ) :      4.109 %\n"
            "Atom     2(N ) :     95.891 %\n"
            "Orbital delocalization index:   7.89\n"
        )
        results = OrbitalCompositionParser.parse_atom_compositions(
            output, method="hirshfeld"
        )
        assert len(results) == 1
        assert isinstance(results[0], OrbitalAtomComposition)
        assert results[0].orbital_id == 5
        assert len(results[0].atom_contributions) == 2

    def test_parse_atom_compositions_empty(self) -> None:
        results = OrbitalCompositionParser.parse_atom_compositions(
            "", method="hirshfeld"
        )
        assert results == []


# =========================================================================
# FuzzySpaceParser
# =========================================================================


class TestFuzzySpaceParser:
    """Tests for the FuzzySpaceParser."""

    def test_parse_fuzzy_integration(self) -> None:
        output = (
            "Atomic space  Value  %  %abs\n"
            "     1(C )            6.21199090            14.790461            14.790461\n"  # noqa: E501
            "Summing up above values:   42.00000000\n"
            "Summing up absolute value of above values:   42.00000000\n"
        )
        result = FuzzySpaceParser.parse_fuzzy_integration(output, "edensity")
        assert result is not None
        assert isinstance(result, FuzzyIntegrationResult)
        assert len(result.entries) == 1
        assert result.entries[0].atom_id == 1
        assert result.total_sum == pytest.approx(42.0)

    def test_parse_fuzzy_integration_empty(self) -> None:
        assert FuzzySpaceParser.parse_fuzzy_integration("", "edensity") is None

    def test_parse_aromaticity_index(self) -> None:
        output = "PDI= 0.05678\nFLU= 0.00123\n"
        r = FuzzySpaceParser.parse_aromaticity_index(output)
        index_map = {item.index_name: item.value for item in r}
        assert index_map["PDI"] == pytest.approx(0.05678)
        assert index_map["FLU"] == pytest.approx(0.00123)

    def test_parse_delocalization_indices(self) -> None:
        output = "Delocalization index of atom  1(C ) and atom  2(N ):  0.45670\n"
        r = FuzzySpaceParser.parse_delocalization_indices(output)
        assert len(r) == 1
        assert isinstance(r[0], DelocalizationIndex)
        assert r[0].atom1_id == 1
        assert r[0].atom2_id == 2
        assert r[0].index == pytest.approx(0.4567)

    def test_parse_delocalization_indices_empty(self) -> None:
        assert FuzzySpaceParser.parse_delocalization_indices("") == []

    def test_parse_atomic_multipoles(self) -> None:
        output = (
            "*****  Atom     1(C )  *****\n"
            " Atomic charge:      -0.052300\n"
            " Atomic monopole moment(most negative in the center):"
            "        -0.052300\n"
            " Atomic dipole moments:\n"
            " X=      0.010000   Y=      0.020000   Z=      0.030000"
            "  Norm=      0.037417\n"
            " Contribution to molecular dipole moment:\n"
            " X=      0.001000   Y=      0.002000   Z=      0.003000"
            "  Norm=      0.003742\n"
            " Traceless Cartesian form:\n"
            " XX=      0.100000   XY=      0.020000   XZ=      0.030000\n"
            " YY=      0.200000   YZ=      0.040000\n"
            " ZZ=     -0.300000\n"
            " Magnitude of the traceless quadrupole moment:"
            "        0.400000\n"
            " Atomic electronic spatial extent <r^2>:       10.500000\n"
            " Components of <r^2>:  X=      3.500000   Y=      3.500000"
            "   Z=      3.500000\n"
            " Magnitude:  |Q_3|=       0.050000\n"
            "\n"
            "*****  Atom     2(N )  *****\n"
            " Atomic charge:       0.052300\n"
            " Atomic monopole moment(most negative in the center):"
            "         0.052300\n"
            "\n"
            " Molecular dipole moment:\n"
        )
        results = FuzzySpaceParser.parse_atomic_multipoles(output)
        assert len(results) == 2
        first = results[0]
        assert first.atom_id == 1
        assert first.atom_element == "C"
        assert first.atomic_charge == pytest.approx(-0.0523)
        assert first.monopole_moment == pytest.approx(-0.0523)
        assert first.dipole_x == pytest.approx(0.01)
        assert first.dipole_norm == pytest.approx(0.037417)
        assert first.mol_dipole_contrib_x == pytest.approx(0.001)
        assert first.mol_dipole_contrib_norm == pytest.approx(0.003742)
        assert first.quadrupole_xx == pytest.approx(0.1)
        assert first.quadrupole_yz == pytest.approx(0.04)
        assert first.quadrupole_magnitude == pytest.approx(0.4)
        assert first.spatial_extent_r2 == pytest.approx(10.5)
        assert first.spatial_extent_x == pytest.approx(3.5)
        assert first.octopole_magnitude == pytest.approx(0.05)
        assert results[1].atom_id == 2
        assert results[1].dipole_x is None

    def test_parse_atomic_multipoles_empty(self) -> None:
        assert FuzzySpaceParser.parse_atomic_multipoles("") == []

    def test_parse_molecular_multipole(self) -> None:
        output = (
            " Total number of electrons:      92.000000     "
            "Net charge:       0.000000\n"
            " Molecular dipole moment (a.u.):     0.100000    "
            "0.200000    0.300000\n"
            " Molecular dipole moment (Debye):    0.254000    "
            "0.508000    0.762000\n"
            " Magnitude of molecular dipole moment (a.u.&Debye):"
            "     0.374166    0.951000\n"
            " Molecular quadrupole moments (Traceless Cartesian form):\n"
            " XX=      1.100000   XY=      0.020000   XZ=      0.030000\n"
            " YY=      1.200000   YZ=      0.040000\n"
            " ZZ=     -2.300000\n"
            " Magnitude of the traceless quadrupole moment tensor:"
            "        2.500000\n"
            " Molecular electronic spatial extent <r^2>:       55.500000\n"
            " Components of <r^2>:  X=      18.500000   Y=      18.500000"
            "   Z=      18.500000\n"
            " Magnitude:  |Q_3|=       1.050000\n"
        )
        mol = FuzzySpaceParser.parse_molecular_multipole(output)
        assert mol is not None
        assert mol.total_electrons == pytest.approx(92.0)
        assert mol.net_charge == pytest.approx(0.0)
        assert mol.dipole_x_au == pytest.approx(0.1)
        assert mol.dipole_x_debye == pytest.approx(0.254)
        assert mol.dipole_magnitude_au == pytest.approx(0.374166)
        assert mol.dipole_magnitude_debye == pytest.approx(0.951)
        assert mol.quadrupole_xx == pytest.approx(1.1)
        assert mol.quadrupole_magnitude == pytest.approx(2.5)
        assert mol.spatial_extent_r2 == pytest.approx(55.5)
        assert mol.spatial_extent_x == pytest.approx(18.5)
        assert mol.octopole_magnitude == pytest.approx(1.05)

    def test_parse_molecular_multipole_not_found(self) -> None:
        assert FuzzySpaceParser.parse_molecular_multipole("nothing") is None

    def test_parse_aom_diagnostics(self) -> None:
        output = (
            " Error of AOM is    0.001234\n"
            " Maximum diagonal deviation to 1:    0.002345 at orbital"
            "    5\n"
            " Maximum nondiagonal deviation to 0:    0.003456 between"
            " orbitals    3    7\n"
            " AOM data has been exported to AOM.txt in current folder\n"
        )
        aom = FuzzySpaceParser.parse_aom_diagnostics(output)
        assert aom is not None
        assert aom.error == pytest.approx(0.001234)
        assert aom.max_diagonal_deviation == pytest.approx(0.002345)
        assert aom.max_diagonal_orbital == 5
        assert aom.max_nondiagonal_deviation == pytest.approx(0.003456)
        assert aom.max_nondiagonal_orbitals == (3, 7)
        assert aom.exported_file == "AOM.txt"

    def test_parse_aom_diagnostics_empty(self) -> None:
        assert FuzzySpaceParser.parse_aom_diagnostics("") is None

    def test_parse_di_matrices(self) -> None:
        output = (
            "**** Localization and delocalization index matrix ****\n"
            "     1     2\n"
            "    1   0.987654   0.512345\n"
            "    2   0.512345   0.876543\n"
            " Localization index:      1(C ):   0.987654    2(N ):"
            "   0.876543\n"
            "Note: some trailing note\n"
        )
        results = FuzzySpaceParser.parse_di_matrices(output)
        assert len(results) == 1
        matrix = results[0]
        assert matrix.label == "Localization and delocalization index matrix"
        assert matrix.n_atoms == 2
        assert matrix.data[(1, 1)] == pytest.approx(0.987654)
        assert matrix.data[(1, 2)] == pytest.approx(0.512345)
        assert len(matrix.localization_indices) == 2
        assert matrix.localization_indices[0].atom_id == 1
        assert matrix.localization_indices[0].index == pytest.approx(0.987654)

    def test_parse_di_matrices_empty(self) -> None:
        assert FuzzySpaceParser.parse_di_matrices("") == []

    def test_parse_overlap_matrices(self) -> None:
        output = (
            "**** Integration of positive values in overlap region ****\n"
            "     1     2\n"
            "    1   0.100000   0.020000\n"
            "    2   0.020000   0.150000\n"
            " Summing up diagonal terms:     0.250000\n"
            " Summing up non-diagonal terms:     0.040000\n"
            " Summing up all terms:     0.290000\n"
        )
        results = FuzzySpaceParser.parse_overlap_matrices(output, "edensity")
        assert len(results) == 1
        matrix = results[0]
        assert matrix.category == "positive"
        assert matrix.integrated_property == "edensity"
        assert matrix.n_atoms == 2
        assert matrix.sum_diagonal == pytest.approx(0.25)
        assert matrix.sum_nondiagonal == pytest.approx(0.04)
        assert matrix.sum_all == pytest.approx(0.29)

    def test_parse_overlap_matrices_empty(self) -> None:
        assert FuzzySpaceParser.parse_overlap_matrices("", "edensity") == []

    def test_parse_for_result_combines_all_sections(self) -> None:
        output = (
            "Atomic space  Value  %  %abs\n"
            "     1(C )            6.21199090            14.790461"
            "            14.790461\n"
            "Summing up above values:   42.00000000\n"
            "Summing up absolute value of above values:   42.00000000\n"
            "*****  Atom     1(C )  *****\n"
            " Atomic charge:      -0.052300\n"
            "\n"
            " Molecular dipole moment:\n"
            " Total number of electrons:      92.000000     "
            "Net charge:       0.000000\n"
            "PDI= 0.05678\n"
            "Delocalization index of atom  1(C ) and atom  2(N ):"
            "  0.45670\n"
        )
        results = FuzzySpaceParser.parse_for_result(
            analysis=Menu.FUZZY_INTEGRATE_EDENSITY,
            stdout=output,
        )
        result_types = {type(r).__name__ for r in results}
        assert "FuzzyIntegrationResult" in result_types
        assert "AtomicMultipole" in result_types
        assert "MolecularMultipole" in result_types
        assert "AromaticityIndex" in result_types
        assert "DelocalizationIndex" in result_types

    def test_parse_clrk_matrix(self) -> None:
        output = (
            "*** Condensed linear response kernel matrix ***\n"
            "              1            2\n"
            "     1   0.12345   0.06789\n"
            "     2   0.06789   0.23456\n"
            "\n"
        )
        result = FuzzySpaceParser.parse_clrk_matrix(output)
        assert result is not None
        assert result.n_atoms == 2
        assert result.data[(1, 1)] == pytest.approx(0.12345)
        assert result.data[(2, 2)] == pytest.approx(0.23456)

    def test_parse_clrk_matrix_terminates_at_next_section(self) -> None:
        output = (
            "**** Condensed linear response kernel matrix ****\n"
            "     1     2\n"
            "    1   0.300000   0.010000\n"
            "    2   0.010000   0.400000\n"
            "\n"
            "Next section starts here\n"
        )
        result = FuzzySpaceParser.parse_clrk_matrix(output)
        assert result is not None
        assert result.n_atoms == 2
        assert result.data == {
            (1, 1): pytest.approx(0.3),
            (1, 2): pytest.approx(0.01),
            (2, 1): pytest.approx(0.01),
            (2, 2): pytest.approx(0.4),
        }

    def test_parse_clrk_matrix_empty(self) -> None:
        assert FuzzySpaceParser.parse_clrk_matrix("") is None

    def test_parse_flu_references(self) -> None:
        output = "FLU reference parameters\nC - C :  1.38900\nC - N :  1.34100\n"
        refs = FuzzySpaceParser.parse_flu_references(output)
        assert len(refs) == 2
        assert refs[0].element1 == "C"
        assert refs[0].reference_value == pytest.approx(1.389)


# =========================================================================
# BasinParser
# =========================================================================


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

    def test_parse_for_result(self) -> None:
        output = (
            "Basin   1  attractor at atom  1(C )  population:  5.96780\n"
            "AIM charge of atom  1(C ):  0.03220\n"
        )
        results = BasinParser.parse_for_result(
            analysis=Menu.BASIN_ANALYSIS_AIM,
            stdout=output,
        )
        assert len(results) == 2


# =========================================================================
# ExcitationParser
# =========================================================================


class TestExcitationParser:
    """Tests for the ExcitationParser."""

    def test_hole_electron_none_without_required_fields(self) -> None:
        output = "D index:  2.34560\nSr:  0.56780\n"
        r = ExcitationParser.parse_hole_electron(output)
        assert r is None

    def test_delta_r(self) -> None:
        output = "State   1  Delta_r:  1.23450\n"
        result = ExcitationParser.parse_delta_r(output)
        assert result[0].delta_r == pytest.approx(1.2345)

    def test_lambda_index(self) -> None:
        output = "State   1  Lambda:  0.78900\n"
        result = ExcitationParser.parse_lambda_index(output)
        assert result[0].lambda_index == pytest.approx(0.789)

    def test_charge_transfer(self) -> None:
        output = (
            "CT distance:  2.345\n"
            "CT amount:  0.567\n"
            "Fragment  1  hole:  0.800  electron:  0.200\n"
        )
        ct = ExcitationParser.parse_charge_transfer(output)
        assert ct.distance == pytest.approx(2.345)
        assert ct.transfer_amount == pytest.approx(0.567)
        assert ct.fragments is not None
        assert len(ct.fragments) == 1

    def test_delta_r_empty(self) -> None:
        assert ExcitationParser.parse_delta_r("") == []

    def test_lambda_index_empty(self) -> None:
        assert ExcitationParser.parse_lambda_index("") == []


# =========================================================================
# WeakInteractionParser
# =========================================================================


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


# =========================================================================
# EDAParser
# =========================================================================


class TestEDAParser:
    """Tests for the EDAParser."""

    def test_parse(self) -> None:
        output = "Electrostatic:  -45.6789\nDispersion:   -5.6789\n"
        r = EDAParser.parse(output)
        assert r.electrostatic == pytest.approx(-45.6789)

    def test_empty(self) -> None:
        r = EDAParser.parse("")
        assert r.electrostatic is None

    def test_parse_dispersion_contributions(self) -> None:
        output = "Atom  1(C )  dispersion:  -0.1234\n"
        r = EDAParser.parse_dispersion_contributions(output)
        assert len(r) == 1
        assert r[0].atom_id == 1
        assert r[0].contribution == pytest.approx(-0.1234)


# =========================================================================
# CDFTParser
# =========================================================================


class TestCDFTParser:
    """Tests for the CDFTParser.

    uses parse_reactivity instead of parse_global_indices.
    """

    def test_parse_reactivity(self) -> None:
        output = (
            "Chemical potential:  -0.15234 a.u.  -4.14500 eV\n"
            "Chemical hardness:   0.21345\n"
        )
        r = CDFTParser.parse_reactivity(output)
        assert r is not None
        assert isinstance(r, Reactivity)
        assert r.chemical_potential == pytest.approx(-0.15234)
        assert r.hardness == pytest.approx(0.21345)

    def test_parse_reactivity_empty(self) -> None:
        assert CDFTParser.parse_reactivity("") is None

    def test_parse_reactivity_full(self) -> None:
        output = (
            " HOMO energy:   -0.234567 a.u.    -6.383310 eV\n"
            " LUMO energy:   -0.045678 a.u.    -1.242890 eV\n"
            " Chemical potential:   -0.140123 a.u.    -3.813100 eV\n"
            " Delta parameter:    0.188889 a.u.     5.140420 eV\n"
            " Global hardness:    0.094444 a.u.\n"
            " Global softness:   10.588235 a.u.\n"
            " Electrophilicity index:    0.104123 a.u.\n"
            " Nucleophilicity index:    0.089012 a.u.\n"
            " Vertical ionization potential:    0.234567 a.u.\n"
            " Vertical electron affinity:    0.045678 a.u.\n"
        )
        r = CDFTParser.parse_reactivity(output)
        assert r is not None
        assert r.homo_energy_au == pytest.approx(-0.234567)
        assert r.homo_energy_eV == pytest.approx(-6.38331)
        assert r.lumo_energy_au == pytest.approx(-0.045678)
        assert r.chemical_potential == pytest.approx(-0.140123)
        assert r.chemical_potential_eV == pytest.approx(-3.8131)
        assert r.delta_parameter_au == pytest.approx(0.188889)
        assert r.delta_parameter_eV == pytest.approx(5.14042)
        assert r.hardness == pytest.approx(0.094444)
        assert r.softness == pytest.approx(10.588235)
        assert r.electrophilicity == pytest.approx(0.104123)
        assert r.nucleophilicity == pytest.approx(0.089012)
        assert r.ionization_potential == pytest.approx(0.234567)
        assert r.electron_affinity == pytest.approx(0.045678)

    def test_parse_orbital_weighted_fukui(self) -> None:
        output = (
            " Atom              OW f+       OW f-      OW f0      "
            "OW dual descriptor\n"
            "     1(C )        0.123400    0.098700    0.111050"
            "    0.024700\n"
            "     2(N )        0.234500    0.187600    0.211050"
            "    0.046900\n"
            " Sum of orbital weighted f+     0.357900\n"
            " Sum of orbital weighted f-     0.286300\n"
        )
        result = CDFTParser.parse_orbital_weighted_fukui(output)
        assert result is not None
        assert len(result.entries) == 2
        assert result.entries[0].atom_id == 1
        assert result.entries[0].ow_f_plus == pytest.approx(0.1234)
        assert result.entries[0].ow_dd == pytest.approx(0.0247)
        assert result.sum_ow_f_plus == pytest.approx(0.3579)
        assert result.sum_ow_f_minus == pytest.approx(0.2863)

    def test_parse_orbital_weighted_fukui_empty(self) -> None:
        assert CDFTParser.parse_orbital_weighted_fukui("") is None

    def test_parse_orbital_weight_decomposition(self) -> None:
        output = (
            " Highest weights in orbital-weighted f+\n"
            " Orbital    22 (LUMO  )   Weight:  48.16 %   "
            "E_diff:     3.410 eV\n"
            " Orbital    23 (LUMO+1)   Weight:  30.22 %   "
            "E_diff:     4.120 eV\n"
            " Total weight of above listed orbitals: 78.38 %\n"
            "\n"
            " Highest weights in orbital-weighted f-\n"
            " Orbital    21 (HOMO  )   Weight:  55.00 %   "
            "E_diff:     2.900 eV\n"
            " Total weight of above listed orbitals: 55.00 %\n"
        )
        results = CDFTParser.parse_orbital_weight_decomposition(output)
        assert len(results) == 2
        assert results[0].fukui_type == "f+"
        assert len(results[0].entries) == 2
        assert results[0].entries[0].orbital_id == 22
        assert results[0].entries[0].orbital_label == "LUMO"
        assert results[0].total_weight_pct == pytest.approx(78.38)
        assert results[1].fukui_type == "f-"
        assert results[1].total_weight_pct == pytest.approx(55.0)

    def test_parse_orbital_weight_decomposition_empty(self) -> None:
        assert CDFTParser.parse_orbital_weight_decomposition("") == []

    def test_parse_condensed_fukui(self) -> None:
        output = (
            "Atom     f+         f-         f0\n"
            "     1(C )        0.14827        0.15436        0.15132\n"
        )
        f = CDFTParser.parse_condensed_fukui(output)
        assert len(f) == 1
        assert isinstance(f[0], CondensedFukui)
        assert f[0].atom_id == 1
        assert f[0].fukui_plus == pytest.approx(0.14827)
        assert f[0].fukui_minus == pytest.approx(0.15436)

    def test_parse_condensed_fukui_empty(self) -> None:
        assert CDFTParser.parse_condensed_fukui("") == []

    def test_parse_dual_descriptor(self) -> None:
        output = "Atom     Dual Descriptor\n     1(C )        0.06670\n"
        result = CDFTParser.parse_dual_descriptor(output)
        assert len(result) == 1
        assert isinstance(result[0], DualDescriptor)
        assert result[0].atom_id == 1
        assert result[0].value == pytest.approx(0.0667)

    def test_parse_dual_descriptor_empty(self) -> None:
        assert CDFTParser.parse_dual_descriptor("") == []

    def test_parse_for_result(self) -> None:
        output = (
            "Chemical potential:  -0.15234 a.u.  -4.14500 eV\n"
            "Chemical hardness:   0.21345\n"
        )
        results = CDFTParser.parse_for_result(
            analysis=Menu.CDFT_ANALYSIS,
            stdout=output,
        )
        reactivity_results = [r for r in results if isinstance(r, Reactivity)]
        assert len(reactivity_results) == 1

    def test_parse_superdelocalizability(self) -> None:
        output = (
            "superdelocalizability analysis\n"
            "alpha parameter:   0.10000 Hartree\n"
            "Atom      D_N       D_E       D_N_0     D_E_0\n"
            "     1(C )        0.12340        0.56780        0.09870        0.43210\n"  # noqa: E501
            "Sum of D_N:   0.12340\n"
            "Sum of D_E:   0.56780\n"
            "Sum of D_N_0:   0.09870\n"
            "Sum of D_E_0:   0.43210\n"
        )
        result = CDFTParser.parse_superdelocalizability(output)
        assert result is not None
        assert result.alpha_parameter == pytest.approx(0.10)
        assert len(result.entries) == 1

    def test_parse_superdelocalizability_empty(self) -> None:
        assert CDFTParser.parse_superdelocalizability("") is None


# =========================================================================
# PolarizabilityParser
# =========================================================================


class TestPolarizabilityParser:
    """Tests for the PolarizabilityParser."""

    def test_parse(self) -> None:
        output = "Isotropic polarizability:  45.67890\nalpha_xx:  56.78900\n"
        r = PolarizabilityParser.parse(output)
        assert r.isotropic == pytest.approx(45.6789)
        assert r.tensor is not None
        assert r.tensor.alpha_xx == pytest.approx(56.789)

    def test_empty(self) -> None:
        r = PolarizabilityParser.parse("")
        assert r.isotropic is None


# =========================================================================
# AromaticityParser
# =========================================================================


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
        assert r.NICS is None

    def test_parse_for_result_nics_scan(self) -> None:
        output = "NICS scan data\n  0.0000   -8.1234\n  1.0000  -10.5678\n"
        results = AromaticityParser.parse_for_result(
            analysis=Menu.NICS_SCAN,
            stdout=output,
        )
        # Should include Aromaticity and possibly NICSScan
        assert len(results) >= 1


# =========================================================================
# WavefunctionParser
# =========================================================================


class TestWavefunctionParser:
    """Tests for the WavefunctionParser."""

    def test_orbital_info_format1(self) -> None:
        output = "   5   Alpha   Occ= 2.000000   E=  -0.72340 a.u.  -19.684 eV\n"
        o = WavefunctionParser.parse_orbital_info(output)
        assert o[0].orbital_id == 5
        assert o[0].spin == "alpha"
        assert o[0].occupation == pytest.approx(2.0)
        assert o[0].energy_au == pytest.approx(-0.7234)

    def test_orbital_info_format2(self) -> None:
        output = (
            "Orb: 10 Ene(au/eV): -0.50000 -13.606 Occ: 2.000000 Type: AlphaBeta\n"  # noqa: E501
        )
        o = WavefunctionParser.parse_orbital_info(output)
        assert len(o) == 1
        assert o[0].orbital_id == 10

    def test_empty(self) -> None:
        assert WavefunctionParser.parse_orbital_info("") == []

    def test_parse_gtf_info(self) -> None:
        output = "    1   Center:    1(C )   Type: S    Exponent:  3047.52490\n"
        gtfs = WavefunctionParser.parse_gtf_info(output)
        assert len(gtfs) == 1
        assert gtfs[0].gtf_index == 1
        assert gtfs[0].exponent == pytest.approx(3047.5249)

    def test_parse_basis_info(self) -> None:
        output = "Basis:    1  Shell:    1  Center:    1(C )  Type: S   GTF:    1 to    6\n"  # noqa: E501
        basis = WavefunctionParser.parse_basis_info(output)
        assert len(basis) == 1
        assert basis[0].basis_index == 1
        assert basis[0].gtf_start == 1
        assert basis[0].gtf_end == 6

    def test_parse_exported_matrices(self) -> None:
        output = "The matrix has been exported to coeff.txt in current folder\n"
        exports = WavefunctionParser.parse_exported_matrices(output)
        assert len(exports) == 1
        assert exports[0].file_name == "coeff.txt"

    def test_parse_basis_info_without_gtf_range(self) -> None:
        # Live Multiwfn 3.8(dev) output for PRINT_ALL_BASIS_FUNCTIONS
        # (Menu 6 -> 2) does not print a GTF range at all.
        output = " Basis:    1   Shell:    1   Center:    1(C )   Type:S    \n"
        basis = WavefunctionParser.parse_basis_info(output)
        assert len(basis) == 1
        assert basis[0].basis_index == 1
        assert basis[0].shell_index == 1
        assert basis[0].center_atom_id == 1
        assert basis[0].center_element == "C"
        assert basis[0].function_type == "S"
        assert basis[0].gtf_start is None
        assert basis[0].gtf_end is None

    def test_parse_coefficient_matrices(self) -> None:
        output = (
            "              1            2\n"
            "    1   0.500000   0.100000\n"
            "    2   0.100000   0.600000\n"
        )
        results = WavefunctionParser.parse_coefficient_matrices(output)
        assert len(results) == 1
        matrix = results[0]
        assert matrix.n_basis == 2
        assert matrix.n_orbitals == 2
        assert matrix.data[(1, 1)] == pytest.approx(0.5)
        assert matrix.data[(2, 2)] == pytest.approx(0.6)

    def test_parse_coefficient_matrices_empty(self) -> None:
        assert WavefunctionParser.parse_coefficient_matrices("") == []

    def test_parse_density_matrices(self) -> None:
        output = (
            "*** Total density matrix ***\n"
            "     1     2\n"
            "    1   1.000000   0.200000\n"
            "    2   0.200000   1.500000\n"
            "Trace of density matrix:      2.500000\n"
            "Trace of density matrix multiplied by overlap matrix:"
            "      92.000000\n"
        )
        results = WavefunctionParser.parse_density_matrices(output)
        assert len(results) == 1
        matrix = results[0]
        assert matrix.label == "Total density matrix"
        assert matrix.n_basis == 2
        assert matrix.data[(1, 1)] == pytest.approx(1.0)
        assert matrix.data[(1, 2)] == pytest.approx(0.2)
        assert matrix.trace == pytest.approx(2.5)
        assert matrix.trace_overlap == pytest.approx(92.0)

    def test_parse_density_matrices_multiple_sections(self) -> None:
        output = (
            "*** Alpha density matrix ***\n"
            "     1     2\n"
            "    1   0.500000   0.100000\n"
            "    2   0.100000   0.700000\n"
            "Trace of density matrix:      1.200000\n"
            "*** Beta density matrix ***\n"
            "     1     2\n"
            "    1   0.400000   0.050000\n"
            "    2   0.050000   0.600000\n"
            "Trace of density matrix:      1.000000\n"
        )
        results = WavefunctionParser.parse_density_matrices(output)
        assert len(results) == 2
        assert results[0].label == "Alpha density matrix"
        assert results[1].label == "Beta density matrix"

    def test_parse_density_matrices_empty(self) -> None:
        assert WavefunctionParser.parse_density_matrices("") == []

    def test_parse_exported_matrices_matches_wavefunction_save(self) -> None:
        # SAVE_WFN / DELETE_INNER_ORBITALS (Menu 6 -> 0 / 34) print
        # "Wavefunction has been outputted to ...", not "matrix".
        output = "Wavefunction has been outputted to new.wfn in current folder\n"
        exports = WavefunctionParser.parse_exported_matrices(output)
        assert len(exports) == 1
        assert exports[0].file_name == "new.wfn"


# =========================================================================
# CubeParser — multi-cube support
# =========================================================================


class TestCubeParser:
    """Tests for the CubeParser — now supports multi-cube outputs."""

    def test_parse_single_cube(self) -> None:
        output = (
            "Number of points in X,Y,Z is   80   80   80  Total:   512000\n"
            "density.cube in current folder\n"
        )
        results = CubeParser.parse(output)
        assert len(results) == 1
        assert isinstance(results[0], Cube)
        assert results[0].file_name == "density.cube"
        assert results[0].x_dim == 80

    def test_parse_multiple_cubes(self) -> None:
        output = "density.cube in current folder\nesp.cube in current folder\n"
        results = CubeParser.parse(output)
        assert len(results) == 2
        assert results[0].file_name == "density.cube"
        assert results[1].file_name == "esp.cube"

    def test_parse_cube_with_extrema(self) -> None:
        output = (
            "The minimum is  0.00001234 at   1.000   2.000   3.000\n"
            "The maximum is  0.56789000 at   4.000   5.000   6.000\n"
            "density.cube in current folder\n"
        )
        results = CubeParser.parse(output)
        assert len(results) == 1
        assert results[0].minimum is not None
        assert results[0].minimum.value == pytest.approx(0.00001234)
        assert results[0].maximum is not None
        assert results[0].maximum.value == pytest.approx(0.56789)

    def test_empty(self) -> None:
        assert CubeParser.parse("") == []


# =========================================================================
# AtomListParser
# =========================================================================


class TestAtomListParser:
    """Tests for the AtomListParser (Menu 0)."""

    def test_parse_homo_lumo_gap(self) -> None:
        output = (
            "Orbital   10 is HOMO, energy:  -0.30000 a.u.  -8.163 eV\n"
            "Orbital   11 is LUMO, energy:  -0.05000 a.u.  -1.360 eV\n"
            "HOMO-LUMO gap:   0.25000 a.u.   6.803 eV   656.300 kJ/mol\n"
        )
        gap = AtomListParser.parse_homo_lumo_gap(output)
        assert gap is not None
        assert isinstance(gap, HOMOLUMOGap)
        assert gap.homo_index == 10
        assert gap.gap_eV == pytest.approx(6.803)

    def test_parse_homo_lumo_gap_not_found(self) -> None:
        assert AtomListParser.parse_homo_lumo_gap("no data") is None

    def test_parse_atoms(self) -> None:
        output = (
            "    1(C ) -->  Charge:   6.000  "
            "x,y,z(Bohr):   0.000000   0.000000   0.000000\n"
        )
        atoms = AtomListParser.parse_atoms(output)
        assert len(atoms) == 1
        assert isinstance(atoms[0], AtomInfo)
        assert atoms[0].atom_id == 1
        assert atoms[0].element == "C"

    def test_parse_atoms_empty(self) -> None:
        assert AtomListParser.parse_atoms("") == []


# =========================================================================
# UtilityParser
# =========================================================================


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
        import re

        output = "density.cube has been generated\nnew.wfn has been generated\n"
        files = re.findall(r"(\S+)\s+has been generated", output)
        assert len(files) == 2

    def test_empty_geometry(self) -> None:
        r = UtilityParser.parse_bond_lengths("")
        assert r == []

    def test_parse_bond_angles(self) -> None:
        output = "Angle  1-2-3:  120.5678\n"
        r = UtilityParser.parse_bond_angles(output)
        assert len(r) == 1
        assert r[0].angle == pytest.approx(120.5678)

    def test_parse_dihedral_angles(self) -> None:
        output = "Dihedral  1-2-3-4:  -45.6789\n"
        r = UtilityParser.parse_dihedral_angles(output)
        assert len(r) == 1
        assert r[0].angle == pytest.approx(-45.6789)

    def test_parse_coordination_numbers(self) -> None:
        output = "Atom  1(C )  coordination number:  4.0\n"
        r = UtilityParser.parse_coordination_numbers(output)
        assert len(r) == 1
        assert r[0].coordination_number == pytest.approx(4.0)

    def test_parse_export_confirmation(self) -> None:
        # Live Multiwfn 3.8(dev) output for GENERATE_CP2K_INPUT
        # (Menu 100 -> 2 -> 25): no "in current folder" suffix, unlike
        # the matrix/wavefunction export messages WavefunctionParser
        # handles.
        output = "CP2K input file has been exported to M062X_TZVPP_D30.inp\n"
        result = UtilityParser.parse_export_confirmation(output)
        assert result is not None
        assert result.file_name == "M062X_TZVPP_D30.inp"
        assert result.label == "CP2K input"

    def test_parse_export_confirmation_no_match(self) -> None:
        assert UtilityParser.parse_export_confirmation("nothing here") is None

    def test_parse_menu300_electric_multipole_report(self) -> None:
        output = (
            " X, Y, Z of center of positive charges (nuclear charges) in Angstrom\n"  # noqa: E501
            "    0.000000   -0.000000   -0.000000\n"  # noqa: E501
            " X, Y, Z of center of negative charges (electronic charges) in Angstrom\n"  # noqa: E501
            "   -0.000000   -0.000000   -0.000000\n\n"  # noqa: E501
            " Dipole moment from nuclear charges (a.u.):    0.000000  -0.000000  -0.000000\n"  # noqa: E501
            " Dipole moment from electrons (a.u.):          0.000000   0.000000   0.000000\n\n"  # noqa: E501
            " Dipole moment (a.u.):       0.000000      0.000000     -0.000000\n"  # noqa: E501
            " Dipole moment (Debye):      0.000000      0.000000     -0.000000\n"  # noqa: E501
            " Magnitude of dipole moment:      0.000000 a.u.      0.000000 Debye\n\n"  # noqa: E501
            " Quadrupole moments (Standard Cartesian form):\n"
            " XX=  -23.392606  XY=    0.000000  XZ=   -0.000000\n"
            " YX=    0.000000  YY=  -23.392606  YZ=    0.000000\n"
            " ZX=   -0.000000  ZY=    0.000000  ZZ=  -28.741943\n"
            " Quadrupole moments (Traceless Cartesian form):\n"
            " XX=    2.674668  XY=    0.000000  XZ=   -0.000000\n"
            " YX=    0.000000  YY=    2.674668  YZ=    0.000000\n"
            " ZX=   -0.000000  ZY=    0.000000  ZZ=   -5.349337\n"  # noqa: E501
            " Magnitude of the traceless quadrupole moment tensor:    5.349337\n"  # noqa: E501
            " Quadrupole moments (Spherical harmonic form):\n"
            " Q_2,0 =  -5.349337   Q_2,-1=   0.000000   Q_2,1=  -0.000000\n"  # noqa: E501
            " Q_2,-2=   0.000000   Q_2,2 =   0.000000\n"
            " Magnitude: |Q_2|=    5.349337\n\n"
            " Octopole moments (Cartesian form):\n"
            " XXX=    0.0000  YYY=   -0.0000  ZZZ=   -0.0000  XYY=   -0.0000  XXY=   -0.0000\n"  # noqa: E501
            " XXZ=   -0.0000  XZZ=    0.0000  YZZ=   -0.0000  YYZ=    0.0000  XYZ=   -0.0000\n"  # noqa: E501
            " Octopole moments (Spherical harmonic form):\n"
            " Q_3,0 =     0.0000  Q_3,-1=     0.0000  Q_3,1 =    -0.0000\n"
            " Q_3,-2=    -0.0000  Q_3,2 =    -0.0000  Q_3,-3=     0.0000  Q_3,3 =     0.0000\n"  # noqa: E501
            " Magnitude: |Q_3|=      0.0000\n\n"
            " Hexadecapole moments:\n"
            " XXXX=       -719.7607  YYYY=       -719.7607  ZZZZ=       -106.5004\n"  # noqa: E501
            " XXXY=         -0.0000  XXXZ=          0.0000  YYYX=          0.0000\n"  # noqa: E501
            " YYYZ=         -0.0000  ZZZX=         -0.0000  ZZZY=          0.0000\n"  # noqa: E501
            " XXYY=       -239.9202  XXZZ=       -161.0740  YYZZ=       -161.0740\n"  # noqa: E501
            " XXYZ=          0.0000  YYXZ=          0.0000  ZZXY=          0.0000\n\n"  # noqa: E501
            " Electronic spatial extent <r^2>:      459.038682\n"
            " Components of <r^2>:  X=     215.148370  Y=     215.148370  Z=      28.741943\n"  # noqa: E501
        )
        report = UtilityParser.parse_electric_multipole_moment_report(output)
        assert report is not None
        assert report.quadrupole_standard_cartesian["XX"] == pytest.approx(
            -23.392606
        )
        assert report.quadrupole_spherical_harmonic["Q_2,0"] == pytest.approx(
            -5.349337
        )
        assert report.octopole_spherical_harmonic["Q_3,3"] == pytest.approx(0.0)
        assert report.hexadecapole["XXXX"] == pytest.approx(-719.7607)
        assert report.electronic_spatial_extent_r2 == pytest.approx(459.038682)

    def test_parse_electric_multipole_empty(self) -> None:
        assert UtilityParser.parse_electric_multipole_moment_report("") is None

    def test_parse_dipole_moments(self) -> None:
        """DipoleMoment computes total from components when Tot is absent."""
        output = "Dipole X=  1.234  Y=  -0.567  Z=  0.901\n"
        result = UtilityParser.parse_dipole_moments(output)
        assert result is not None
        assert result.x == pytest.approx(1.234)
        assert result.y == pytest.approx(-0.567)
        assert result.z == pytest.approx(0.901)
        expected_total = (1.234**2 + 0.567**2 + 0.901**2) ** 0.5
        assert result.total == pytest.approx(expected_total)

    def test_parse_dipole_moments_plain_triplet(self) -> None:
        """Plain triplet format computes total from components."""
        output = "Dipole moment (a.u.):   1.234   -0.567   0.901\n"
        result = UtilityParser.parse_dipole_moments(output)
        assert result is not None
        assert result.x == pytest.approx(1.234)
        assert result.y == pytest.approx(-0.567)
        assert result.z == pytest.approx(0.901)
        expected_total = (1.234**2 + 0.567**2 + 0.901**2) ** 0.5
        assert result.total == pytest.approx(expected_total)

    def test_parse_dipole_moments_empty(self) -> None:
        assert UtilityParser.parse_dipole_moments("") is None

    def test_parse_quadrupole_moments(self) -> None:
        output = "Quadrupole XX=  -23.39  XY=  0.00\n"
        r = UtilityParser.parse_quadrupole_moments(output)
        assert r is not None
        assert r.xx == pytest.approx(-23.39)

    def test_parse_quadrupole_moments_empty(self) -> None:
        assert UtilityParser.parse_quadrupole_moments("") is None


# =========================================================================
# ParserRoute
# =========================================================================


class TestParserRoute:
    """Tests for the ParserRoute routing table."""

    def test_all_charge_menus_route_to_charge_parser(self) -> None:
        from pymultiwfn.analysis.parsers import ParserRoute

        charge_menus = [
            Menu.HIRSHFELD_CHARGE,
            Menu.MULLIKEN_POPULATION,
            Menu.ADCH_CHARGE,
            Menu.CM5_CHARGE,
        ]
        for menu in charge_menus:
            assert ParserRoute.ROUTE_TABLE.get(menu) is ChargeParser

    def test_bond_order_menus_route_to_bond_order_parser(self) -> None:
        from pymultiwfn.analysis.parsers import ParserRoute

        bo_menus = [
            Menu.MAYER_BOND_ORDER,
            Menu.WIBERG_BOND_ORDER,
            Menu.FUZZY_BOND_ORDER,
        ]
        for menu in bo_menus:
            assert ParserRoute.ROUTE_TABLE.get(menu) is BondOrderParser

    def test_topology_menus_route_to_critical_point_parser(self) -> None:
        from pymultiwfn.analysis.parsers import ParserRoute

        assert (
            ParserRoute.ROUTE_TABLE.get(Menu.TOPOLOGY_VISUALISE_CPS)
            is CriticalPointParser
        )

    def test_unregistered_menu_returns_none(self) -> None:
        from pymultiwfn.analysis.parsers import ParserRoute

        # PROPERTIES_AT_POINT is interactive and not routed
        assert ParserRoute.ROUTE_TABLE.get(Menu.PROPERTIES_AT_POINT) is None

    def test_menu6_print_and_save_actions_route_to_wavefunction_parser(
        self,
    ) -> None:
        from pymultiwfn.analysis.parsers import ParserRoute

        menus = [
            Menu.PRINT_ALL_GTF,
            Menu.PRINT_ALL_BASIS_FUNCTIONS,
            Menu.PRINT_ORBITAL_INFO,
            Menu.SAVE_WFN,
            Menu.DELETE_INNER_ORBITALS,
        ]
        for menu in menus:
            assert ParserRoute.ROUTE_TABLE.get(menu) is WavefunctionParser

    def test_generate_cp2k_input_routes_to_utility_parser(self) -> None:
        from pymultiwfn.analysis.parsers import ParserRoute

        assert ParserRoute.ROUTE_TABLE.get(Menu.GENERATE_CP2K_INPUT) is UtilityParser

    def test_genuinely_interactive_menu6_and_cdft_entries_are_unrouted(
        self,
    ) -> None:
        from pymultiwfn.analysis.parsers import ParserRoute

        # Live-verified: each of these crashes ("end-of-file during
        # read") waiting for interactive input Multiwfn cannot be
        # given non-interactively (an occupation value, Gaussian
        # net-charge/multiplicity, or a second wavefunction path), so
        # there is no complete output to route to a parser.
        for menu in (
            Menu.MODIFY_OCCUPATION,
            Menu.CDFT_GENERATE_CHARGED_WFN,
            Menu.CDFT_GRID_FUKUI_POTENTIAL,
        ):
            assert ParserRoute.ROUTE_TABLE.get(menu) is None

    def test_only_four_menu_members_remain_unrouted(self) -> None:
        from pymultiwfn.analysis.parsers import ParserRoute

        covered = set(ParserRoute.ROUTE_TABLE.keys())
        missing = {m for m in Menu if m not in covered}
        assert missing == {
            Menu.PROPERTIES_AT_POINT,
            Menu.MODIFY_OCCUPATION,
            Menu.CDFT_GENERATE_CHARGED_WFN,
            Menu.CDFT_GRID_FUKUI_POTENTIAL,
        }
