"""Physical & thermodynamic regression tests for aska-physical/funi-frame.

Covers:
  1. Exact parabolic arch-tie virtual-work crown deflection (replacing the legacy 0.18 beam factor):
     includes dead load, arch shortening, and tie elongation, proves virtual-work invariance, and
     certifies both active Fe-SMA pass (34.96 mm, L/686.4) and un-prestressed / post-activation
     incremental live-load L/360 failures (152.15 mm / 128.18 mm).
  2. Film-convection-limited mPCM wall ride-through ODE (q = h*(T - T_pc)) and isothermal finite-NTU
     PCM skid heat-transfer rate (Q = eps * mdot*c_p*dT): proves passive natural convection at
     h = 8 W/(m^2*K) hits 35 C in 2.08 min (1.14% sensible+latent utilization; needs h >= 20.55 W/(m^2*K))
     and that finite NTU = 1.5 (eps = 0.7769) caps heat rejection at 7.769 MW < 10.0 MW, discharging
     the full 2,500 kWh_th inventory in 83.33 m^3 over 19.31 min (requires 1.2872x flow capacity or higher UA).
  3. State-dependent sorption COP_ads(T_a) = eta_II * COP_rev(T_h, T_a, T_c): proves COP collapses
     from 0.620 at 30 C to 0.2174 at 40 C (alpha* = 0.8214, Q_drive = 6.900 MW) and 0.1279 at 45 C
     (fails self-cooling at alpha = 0.85; ceiling T_a <= 42.04 C at fixed eta_II = 0.3022).

Runs stdlib-only across ``-I -B``, ``-I -O -B``, and ``-I -B -W error`` (zero bare ``assert``s).
"""

from __future__ import annotations

import ast
import importlib
import math
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _import_verify():
    if "verify" in sys.modules:
        return sys.modules["verify"]
    try:
        return importlib.import_module("verify")
    except ModuleNotFoundError as exc:
        if exc.name != "sympy":
            raise
        stub = types.ModuleType("sympy")
        with mock.patch.dict(sys.modules, {"sympy": stub}):
            return importlib.import_module("verify")


class TestFuniFramePhysicalHardening(unittest.TestCase):
    def test_roof_virtual_work_deflection_replaces_0p18_and_exposes_tie_stretch(self) -> None:
        verify = _import_verify()
        rf = verify.solve_funicular_roof()
        # Exact virtual-work invariance between cut-tie and tied virtual force systems
        self.assertLess(rf["virtual_work_invariance_err_mm"], 1e-6)
        # Active Fe-SMA net service crown deflection (dead + plant + hung) is 34.96 mm, NOT 11.74 mm
        self.assertAlmostEqual(rf["delta_vault_mm"], 34.96, places=2)
        self.assertGreaterEqual(rf["span_over_delta_vault"], 360.0)
        # Honest tie-stretch crown deflections: un-prestressed (152.15 mm) and post-activation live (128.18 mm) fail L/360
        self.assertAlmostEqual(rf["delta_vault_unprestressed_mm"], 152.15, places=2)
        self.assertLess(rf["span_over_delta_unprestressed"], 360.0)
        self.assertAlmostEqual(rf["delta_vault_incremental_live_mm"], 128.18, places=2)
        self.assertLess(rf["span_over_delta_incremental_live"], 360.0)
        # Turning off SMA activation reproduces unprestressed deflection identically
        cold = verify.solve_funicular_roof(T_sma_C=25.0)
        self.assertAlmostEqual(cold["delta_vault_mm"], rf["delta_vault_unprestressed_mm"], places=2)
        for bad in (0.0, -5.0, float("nan"), float("inf"), True):
            with self.assertRaises(ValueError):
                verify.solve_funicular_roof(q_plant=bad)
            with self.assertRaises(ValueError):
                verify.solve_funicular_roof(q_hung=bad)

    def test_pcm_wall_convection_limit_and_forced_airflow_threshold(self) -> None:
        verify = _import_verify()
        wl = verify.solve_vibro_thermal_wall()
        self.assertEqual(wl["h_conv_passive_W_m2K"], 8.0)
        self.assertAlmostEqual(wl["wall_heat_flux_W_m2"], 750000.0 / 5215.0, places=2)
        self.assertAlmostEqual(wl["h_conv_min_for_35C_hold_W_m2K"], 750000.0 / (5215.0 * 7.0), places=2)
        # At h = 8 W/(m^2*K), equilibrium is ~45.98 C >> 35 C and aisle hits 35 C in ~2.08 min (~1.14% latent used)
        self.assertGreater(wl["T_aisle_equilibrium_passive_C"], 45.0)
        self.assertLess(wl["passive_ride_through_convection_limited_min"], 3.0)
        self.assertGreater(wl["passive_ride_through_convection_limited_min"], wl["passive_ride_through_air_only_min"])
        self.assertLess(wl["mpcm_latent_utilization_before_35C_pct"], 2.0)
        self.assertAlmostEqual(wl["passive_ride_through_bucket_ideal_min"], 28.30, places=2)
        # Forced convection above h_min unlocks the full bucket ceiling
        forced = verify.solve_vibro_thermal_wall(h_conv_passive_W_m2K=25.0)
        self.assertLessEqual(forced["T_aisle_equilibrium_passive_C"], 35.0)
        self.assertEqual(forced["passive_ride_through_convection_limited_min"], forced["passive_ride_through_bucket_ideal_min"])
        self.assertEqual(forced["mpcm_latent_utilization_before_35C_pct"], 100.0)
        for bad in (0.0, -8.0, float("nan"), float("inf"), True):
            with self.assertRaises(ValueError):
                verify.solve_vibro_thermal_wall(h_conv_passive_W_m2K=bad)

    def test_sorption_cop_recooling_dependence_and_pcm_skid_ntu(self) -> None:
        verify = _import_verify()
        base = verify.solve_three_tier_tes_and_cooling()
        self.assertAlmostEqual(base["cop_ads_operating"], 0.62, places=4)
        self.assertAlmostEqual(base["cop_ads_second_law_eta_II"], 0.3022, places=4)
        self.assertAlmostEqual(base["cop_ads_at_40C_recooling"], 0.2174, places=4)
        self.assertAlmostEqual(base["alpha_liq_critical_at_40C"], 0.8214, places=4)
        self.assertAlmostEqual(base["Q_drive_at_40C_recooling_MW"], 6.900, places=3)
        self.assertFalse(base["adsorption_self_cools_at_45C_recooling"])
        self.assertAlmostEqual(base["T_recooling_max_for_self_cooling_at_eta_II_C"], 42.04, places=2)

        # Calling solver directly at T_a = 40 C and 45 C dynamically recomputes COP, Q_drive, and self-cooling
        at40 = verify.solve_three_tier_tes_and_cooling(T_recooling_design_C=40.0)
        self.assertAlmostEqual(at40["cop_ads_operating"], base["cop_ads_at_40C_recooling"], places=4)
        self.assertAlmostEqual(at40["Q_drive_for_adsorption_MW"], base["Q_drive_at_40C_recooling_MW"], places=3)
        self.assertTrue(at40["adsorption_self_cools_air_load"])

        at45 = verify.solve_three_tier_tes_and_cooling(T_recooling_design_C=45.0)
        self.assertAlmostEqual(at45["cop_ads_operating"], base["cop_ads_at_45C_recooling"], places=4)
        self.assertFalse(at45["adsorption_self_cools_air_load"])

        # Isothermal HX NTU rate-vs-inventory physics: eps scales heat-transfer rate (7.769 MW < 10 MW),
        # discharging the full 2,500 kWh_th inventory over 19.31 min rather than destroying stored energy
        expected_eps = 1.0 - math.exp(-1.5)
        self.assertAlmostEqual(base["pcm_skid_effectiveness"], round(expected_eps, 4), places=4)
        self.assertAlmostEqual(base["pcm_skid_peak_cooling_rate_MW"], round(10.0 * expected_eps, 3), places=3)
        self.assertAlmostEqual(base["pcm_skid_rate_shortfall_MW"], round(10.0 * (1.0 - expected_eps), 3), places=3)
        self.assertAlmostEqual(base["pcm_skid_discharge_duration_at_hx_rate_min"], round((2.5 / (10.0 * expected_eps)) * 60.0, 2), places=2)
        self.assertAlmostEqual(base["pcm_skid_flow_capacity_multiplier_for_full_rate"], round(1.0 / expected_eps, 4), places=4)
        self.assertFalse(base["pcm_skid_meets_full_it_rate_at_nominal_flow"])
        self.assertAlmostEqual(
            base["pcm_skid_peak_cooling_rate_MW"] * (base["pcm_skid_discharge_duration_at_hx_rate_min"] / 60.0),
            2.5,
            delta=0.01,
        )

        # Out-of-domain temperatures & NTU fail closed
        for bad_t in (18.0, 10.0, 58.0, 70.0, float("nan"), float("inf"), True):
            with self.assertRaises(ValueError):
                verify.solve_three_tier_tes_and_cooling(T_recooling_design_C=bad_t)
        for bad_ntu in (0.0, -1.5, float("nan"), float("inf"), True):
            with self.assertRaises(ValueError):
                verify.solve_three_tier_tes_and_cooling(pcm_skid_ntu=bad_ntu)

    def test_zero_bare_asserts(self) -> None:
        for rel in ("verify.py", "tests/test_physical_hardening.py"):
            tree = ast.parse((REPO_ROOT / rel).read_text(encoding="utf-8"))
            self.assertEqual([n.lineno for n in ast.walk(tree) if isinstance(n, ast.Assert)], [])


if __name__ == "__main__":
    unittest.main()
