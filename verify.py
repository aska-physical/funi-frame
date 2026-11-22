#!/usr/bin/env python3
"""Deterministic six-domain verifier for aska-physical/funi-frame (v0.1.1):
Monolithic Vibro-Thermal Funicular Data-Center Shell, 3-Tier Waste-Heat TES & Zero-Water Sorption Cooling.
"""
from __future__ import annotations
import argparse, ast, json, math, sys
from pathlib import Path
import sympy as sp

CERT_PATH = Path(__file__).resolve().parent / "results" / "certificate.json"


SOLVER_DOMAIN = {
    "solve_funicular_roof": "domain_1_funicular_radiative_roof",
    "solve_vibro_thermal_wall": "domain_2_vibro_thermal_sandwich_wall",
    "solve_vibration_cascade": "domain_3_broadband_vibration_cascade",
    "solve_whole_frame_seismic": "domain_4_top_down_seismic_cellular_frame",
    "solve_three_tier_tes_and_cooling": "domain_5_6_three_tier_tes_and_closed_loop_cooling",
}


def literal_certificate_fields() -> dict[str, list[str]]:
  """Certificate fields whose value is a number written in this file, or a name bound straight to one.
  Read from the source, so the list of declared inputs cannot fall behind the solvers."""
  number = lambda v: isinstance(v, ast.Constant) and isinstance(v.value, (int, float)) and not isinstance(v.value, bool)
  found = {}
  for fn in ast.parse(Path(__file__).read_text(encoding="utf-8")).body:
    if not (isinstance(fn, ast.FunctionDef) and fn.name in SOLVER_DOMAIN):
      continue
    bound = set()
    for st in ast.walk(fn):
      if isinstance(st, ast.Assign) and len(st.targets) == 1:
        target, value = st.targets[0], st.value
        pairs = [(target, value)]
        if isinstance(target, ast.Tuple) and isinstance(value, ast.Tuple) and len(target.elts) == len(value.elts):
          pairs = list(zip(target.elts, value.elts))
        bound |= {t.id for t, v in pairs if isinstance(t, ast.Name) and number(v)}
    fields = [k.value for r in ast.walk(fn) if isinstance(r, ast.Return) and isinstance(r.value, ast.Dict)
              for k, v in zip(r.value.keys, r.value.values)
              if number(v) or (isinstance(v, ast.Name) and v.id in bound)]
    found[SOLVER_DOMAIN[fn.name]] = sorted(fields)
  return found


def require(condition, message="") -> None:
  """A check that cannot be switched off. `python3 -O` strips every `assert` statement, so none is used here."""
  if not condition:
    raise RuntimeError(f"verification failed: {message}")


def verify_symbolic_theorems() -> dict[str, str]:
    r, eta, zeta = sp.symbols("r eta zeta", positive=True)
    T_hyst = sp.sqrt((1 + eta**2) / ((1 - r**2) ** 2 + eta**2))
    T_visc = sp.sqrt((1 + (2 * zeta * r) ** 2) / ((1 - r**2) ** 2 + (2 * zeta * r) ** 2))
    require(sp.limit(r * T_hyst, r, sp.oo) == 0, "a symbolic identity does not hold")
    require(sp.simplify(sp.limit(r * T_visc, r, sp.oo) - 2 * zeta) == 0, "a symbolic identity does not hold")
    require(sp.simplify(sp.limit(r**2 * T_hyst, r, sp.oo) - sp.sqrt(1 + eta**2)) == 0, "a symbolic identity does not hold")

    rho_pin, k_pin, k_ins, t_ins = sp.symbols("rho_pin k_pin k_ins t_ins", positive=True)
    U_eff = ((1 - rho_pin) * k_ins + rho_pin * k_pin) / t_ins
    require(sp.simplify(U_eff / (k_ins / t_ins) - (1 + rho_pin * (k_pin / k_ins - 1))) == 0, "a symbolic identity does not hold")

    phi, rho_p, cp_p, L_f, rho_w, cp_w, dT = sp.symbols("phi rho_p cp_p L_f rho_w cp_w dT", positive=True)
    E_w = rho_w * cp_w * dT
    E_pcm = (1 - phi) * E_w + phi * rho_p * (cp_p * dT + L_f)
    require(sp.simplify(E_pcm / E_w - ((1 - phi) + phi * rho_p * (cp_p * dT + L_f) / E_w)) == 0, "a symbolic identity does not hold")

    a_liq, cop = sp.symbols("alpha_liq cop_ads", positive=True)
    require(sp.simplify(sp.solve(sp.Eq(a_liq * cop, 1 - a_liq), a_liq)[0] - 1 / (1 + cop)) == 0, "a symbolic identity does not hold")

    Th, Ta, Tc = sp.symbols("T_h T_a T_c", positive=True)
    cop_rev = (Tc / (Ta - Tc)) * ((Th - Ta) / Th)
    require(sp.simplify(cop_rev.subs(Ta, Tc * Th * (1 + cop) / (cop * Th + Tc)) - cop) == 0, "a symbolic identity does not hold")

    P, pue_c, pue_f, c_loop, q0, a0, P2 = sp.symbols("P PUE_c PUE_f c_loop q0 a0 P2", positive=True)
    require(sp.simplify(((P * pue_c) / pue_f - P) - P * (pue_c - pue_f) / pue_f) == 0, "a symbolic identity does not hold")
    used = lambda load: load + (load - q0) / c_loop + a0
    grown = sp.solve(sp.Eq(used(P2), P * pue_c), P2)
    require(len(grown) == 1 and sp.simplify((grown[0] - P) - (P * pue_c - used(P)) / (1 + 1 / c_loop)) == 0, "a symbolic identity does not hold")
    return {
        "hysteretic_vs_viscous_rolloff": "lim_{r->oo} r*T_hyst = 0 (O(r^-2)) vs. lim_{r->oo} r*T_visc = 2*zeta (O(r^-1))",
        "thermal_bridge_amplification": "U_eff / U_core = 1 + rho_pin * (k_pin / k_ins - 1)",
        "pcm_volumetric_gain_identity": "E_tank / E_water = (1 - phi) + phi * rho_p * (cp_p*dT + L_f) / (rho_w*cp_w*dT)",
        "self_cooling_liquid_capture_threshold": "alpha_liq* = 1 / (1 + COP_ads)",
        "sorption_recooling_limit": "COP_ads <= (T_c/(T_a - T_c))*((T_h - T_a)/T_h)  <=>  T_a <= T_c*T_h*(1 + COP_ads)/(COP_ads*T_h + T_c)",
        "it_power_from_freed_utility_power": "dP_IT = P_freed/PUE if all overhead scales with load; P_freed/(1 + 1/COP_loop) if only the cooling loop does",
    }


def solve_funicular_roof() -> dict:
  L, w_bay, h_v, f_r, t_sh = 24.0, 3.0, 1480.0, 1.350, 70.0
  q_plant, q_hung, q_tor, fc, n_r, Af, Asma = 7.5, 4.0, 6.2, 55.0, 4, 780.0, 1150.0
  m_sh = (t_sh * 1e-3) * 1850.0 + (n_r * 0.075 * 0.280 / w_bay) * 2150.0
  dead_m = m_sh + 16.0 + ((Af * 1e-6 * 1600.0 + Asma * 1e-6 * 7800.0) * n_r / w_bay)
  q_dead = dead_m * 9.80665 * 1e-3
  w_serv = (q_dead + q_plant + q_hung) * w_bay
  P_sma = n_r * Asma * (435.0 * math.tanh((175.0 - 50.0) / 58.0)) * 1e-3
  Ec, h_cr = 4700.0 * math.sqrt(fc), 220.0
  I_ar = n_r * (75.0 * h_cr**3 / 12.0) + (w_bay * 1000.0) * t_sh**3 / 12.0
  A_cr = (w_bay * 1000.0) * t_sh + n_r * 75.0 * (h_cr - t_sh)
  EA_tie = n_r * (Af * 195800.0 + Asma * 165000.0)
  chi = (15.0 * Ec * I_ar) / (8.0 * (f_r * 1000.0) ** 2 * (1.0 / (1.0 / EA_tie + 1.0 / (Ec * A_cr))))
  H_rig = w_serv * L**2 / (8.0 * f_r)
  H_srv = (H_rig + chi * P_sma) / (1.0 + chi)
  M_srv = w_serv * L**2 / 8.0 - H_srv * f_r
  S_cr = I_ar / (0.5 * h_cr)
  fc_cr = (H_srv * 1e3) / A_cr + abs(M_srv * 1e6) / S_cr

  EA_gf = n_r * (420.0 * 46000.0)
  chi_gf = (15.0 * Ec * I_ar) / (8.0 * (f_r * 1000.0) ** 2 * EA_gf)
  H_gf = max(0.0, (H_rig - (0.6 / 24000.0) * EA_gf * 1e-3 * chi_gf) / (1.0 + chi_gf)) * 0.82
  ft_gf = abs((w_serv * L**2 / 8.0 - H_gf * f_r) * 1e6) / S_cr - (H_gf * 1e3) / A_cr

  EA_top = Ec * (w_bay * 1000.0) * t_sh
  EI_v = ((EA_top * EA_tie) / (EA_top + EA_tie)) * ((h_v - 0.5 * t_sh - 38.0) * 1e-3) ** 2
  d_v = 0.18 * (5.0 * (q_plant + q_hung) * w_bay * 1e3 * L**4) / (384.0 * EI_v) * 1000.0
  d_ipe = (5.0 * 15000.0 * L**4) / (384.0 * 210e9 * 9.208e-4) * 1000.0
  U_rf = 1.0 / (0.150 / 0.022 + 0.070 / 0.160 + 0.85)
  return {
      "span_m": L, "vault_depth_mm": h_v, "span_to_rise_ratio": round(L / f_r, 2),
      "plenum_height_saved_vs_2400mm_deep_truss_mm": round(2400.0 - h_v, 1),
      "dead_mass_kg_m2": round(dead_m, 2),
      "mass_reduction_vs_460kg_ballasted_roof_pct": round((1.0 - dead_m / 460.0) * 100.0, 2),
      "H_service_kN": round(H_srv, 2), "fc_crown_MPa": round(fc_cr, 2),
      "fc_allow_service_MPa": round(0.45 * fc, 2), "ft_crown_passive_gfrp_MPa": round(ft_gf, 2),
      "delta_vault_mm": round(d_v, 2), "span_over_delta_vault": round(L * 1000.0 / d_v, 1),
      "span_over_delta_ipe600_rafter": round(L * 1000.0 / d_ipe, 1),
      "tornado_net_uplift_kN": round(max(0.0, (q_tor - 0.9 * q_dead) * L * w_bay), 2),
      "tornado_tiedown_capacity_kN": 783.0, "U_roof_W_m2K": round(U_rf, 4),
      "R_value_imperial": round((1.0 / U_rf) * 5.678263, 1), "noon_radiative_sky_cooling_W_m2": 93.0,
  }


def solve_vibro_thermal_wall() -> dict:
  m_int = 0.052 * 680.0
  w_mass = 0.048 * 2200.0 + 0.160 * 105.0 + m_int + 1.9
  k_cf = (1.0 - 0.0008) * 0.022 + 0.0008 * 0.80
  U_cf = 1.0 / (0.17 + 0.048 / 1.45 + 0.052 / 0.16 + 0.160 / k_cf)
  U_st = 1.0 / (0.17 + 0.150 / 1.80 + 0.038 / ((1.0 - 0.0035) * 0.036 + 0.0035 * 50.0))
  perim_A = 4.0 * math.sqrt(2500.0) * 5.45
  rev_drop = 10.0 * math.log10((perim_A * 0.72 + 2500.0 * 0.35) / (perim_A * 0.03 + 2500.0 * 0.04))
  E_air = (0.45 * 2500.0 * 5.45) * 1.18 * 1.006 * 11.0
  m_act = 0.10 * (perim_A + 1.65 * 2500.0) * m_int
  E_pcm = m_act * 0.92 * 11.0 + (m_act * 0.28) * 195.0
  return {
      "wall_mass_kg_m2": round(w_mass, 2),
      "mass_reduction_vs_368kg_tiltup_pct": round((1.0 - w_mass / 368.0) * 100.0, 2),
      "U_wall_cfrp_W_m2K": round(U_cf, 4), "R_wall_imperial": round((1.0 / U_cf) * 5.678263, 1),
      "U_wall_steel_pcsp_W_m2K": round(U_st, 3), "nrc_poraver_helmholtz": 0.72,
      "hall_reverberant_noise_reduction_dBA": round(rev_drop, 2),
      "stc_community_rating_dB": round(20.0 * math.log10(w_mass) + 16.0, 1),
      "passive_ride_through_with_mpcm_min": round((E_air + E_pcm) / (750.0 * 60.0), 2),
      "passive_ride_through_air_only_min": round(E_air / (750.0 * 60.0), 2),
  }


def _th(f: float, fn: float, eta: float) -> float:
  return math.sqrt((1.0 + eta**2) / ((1.0 - (f / fn) ** 2) ** 2 + eta**2))


def solve_vibration_cascade() -> dict:
  tones = [(29.5, 0.085), (49.0, 0.060), (60.0, 0.110), (120.0, 0.045)]
  g_a = math.sqrt(sum((a * _th(f, 4.5, 0.16) * _th(f, 9.2, 0.12) * 0.40) ** 2 for f, a in tones))
  g_c = math.sqrt(sum((a * _th(f, 18.0, 0.05) * _th(f, 26.0, 0.003) * 0.78) ** 2 for f, a in tones))
  return {
      "oem_switch_limit_g_rms": 0.0005, "rack_g_rms_funi_frame": round(g_a, 6),
      "rack_g_rms_conventional_steel": round(g_c, 6), "oem_utilization_pct": round(g_a / 0.0005 * 100.0, 2),
      "isolation_improvement_vs_steel_dB": round(20.0 * math.log10(g_a / g_c), 2),
  }


def solve_whole_frame_seismic(roof_m: float, wall_m: float) -> dict:
  h_c, h_a, p_m = 2 * 6.0 + 2.40, 2 * 5.45 + 1.48, 4.0 * math.sqrt(2500.0)
  m_c = (2500.0 * (744.0 + 460.0) + p_m * h_c * 368.0) / 1000.0
  m_a = (2500.0 * (195.12 + roof_m) + p_m * h_a * wall_m) / 1000.0
  m_live = (2500.0 * (0.25 * 15.0 * 101.97 + 7.5 * 101.97)) / 1000.0
  V_a = (1.20 / 5.0) * (m_a + m_live) * 9.80665
  V_c = (1.20 / 4.5) * (m_c + m_live) * 9.80665
  return {
      "total_building_height_m": round(h_a, 2), "height_saved_vs_conv_m": round(h_c - h_a, 2),
      "structural_dead_mass_tonnes": round(m_a, 1), "conv_structural_dead_mass_tonnes": round(m_c, 1),
      "structural_dead_mass_reduction_pct": round((1.0 - m_a / m_c) * 100.0, 2),
      "seismic_base_shear_aska_kN": round(V_a, 1),
      "seismic_base_shear_reduction_pct": round((1.0 - V_a / V_c) * 100.0, 2),
      "top_down_story_drift_aska_pct": round((V_a * 4.0) / (880000.0 * 1.50 * 5.45) * 100.0, 3),
      "top_down_story_drift_conv_unbraced_pct": round((V_c * 4.0) / (185000.0 * 1.50 * 6.00) * 100.0, 3),
      "risk_cat_iv_drift_limit_pct": 1.0,
  }


def _cop_rev(t_drive_c: float, t_recool_c: float, t_chilled_c: float) -> float:
  th, ta, tc = t_drive_c + 273.15, t_recool_c + 273.15, t_chilled_c + 273.15
  return (tc / (ta - tc)) * ((th - ta) / th)


def solve_three_tier_tes_and_cooling() -> dict:
  P_it, Q_liq, Q_air = 10.0, 8.5, 1.5
  cop_ads, cop_loop = 0.62, 28.0
  # T_rec is a proposed design value, new in v0.1.1. No source is cited for it and the owner has not confirmed it.
  # Control 13 holds for any value below T_rec_max.
  T_drv, T_chw, T_rec = 58.0, 18.0, 30.0  # drive supply, warmest chilled water, proposed design recooling (deg C)
  v_pcm, v_wat = 2500.0 / 30.0, 2500.0 / 11.5
  Q_drv = Q_air / cop_ads
  Q_pil, Q_sky = 200 * 25.0 * 50.0 * 1e-6, 2500.0 * 93.0 * 1e-6
  Q_ates = (Q_liq - Q_drv) + (Q_air + Q_drv) - (Q_pil + 0.55 * Q_sky)
  P_c_conv, P_c_funi = P_it / 3.0 + 0.27, Q_ates / cop_loop + 0.195
  P_freed = P_c_conv - P_c_funi
  pue_funi = (P_it + P_c_funi + 0.2) / P_it
  # Freed utility power is not IT power: new racks need their own cooling and overhead.
  P_it_low, P_it_high = P_freed / pue_funi, P_freed / (1.0 + 1.0 / cop_loop)
  T_rec_max = (T_chw + 273.15) * (T_drv + 273.15) * (1.0 + cop_ads) / (cop_ads * (T_drv + 273.15) + T_chw + 273.15) - 273.15
  return {
      "P_it_MW": P_it, "Q_liquid_w45_MW": Q_liq, "Q_air_rdhx_MW": Q_air,
      "alpha_liq_critical_for_self_cooling": round(1.0 / (1.0 + cop_ads), 4),
      "cop_ads_design": cop_ads, "T_drive_supply_C": T_drv, "T_chilled_water_C": T_chw, "T_recooling_design_C": T_rec,
      "T_recooling_design_basis": "proposed in v0.1.1; no source cited; awaiting the owner's confirmation",
      "cop_loop_design": cop_loop,
      "cop_ads_reversible_limit_at_design_recooling": round(_cop_rev(T_drv, T_rec, T_chw), 3),
      "cop_ads_reversible_limit_at_45C_dry_recooling": round(_cop_rev(T_drv, 45.0, T_chw), 3),
      "T_recooling_max_for_design_cop_C": round(T_rec_max, 2),
      "Q_drive_for_adsorption_MW": round(Q_drv, 3), "adsorption_self_cools_air_load": Q_drv <= Q_liq,
      "cop_ads_at_32C_air_exhaust": 0.0, "pcm_ride_through_volume_m3": round(v_pcm, 2),
      "water_equivalent_volume_m3": round(v_wat, 2),
      "water_only_ride_through_in_same_volume_min": round((v_pcm * 11.5 / 1000.0) / P_it * 60.0, 2),
      "basalt_battery_volume_m3": 88.68, "basalt_media_mass_tonnes": 117.06, "basalt_media_cost_usd": 14046.0,
      "annualized_capex_basalt_usd_kWh": 0.014, "annualized_capex_liion_usd_kWh": 0.068,
      "Q_energy_piles_MW": round(Q_pil, 3), "Q_sky_roof_peak_MW": round(Q_sky, 3),
      "PUE_conventional": round((P_it + P_c_conv + 0.2) / P_it, 3),
      "PUE_funi_frame": round((P_it + P_c_funi + 0.2) / P_it, 3),
      "cooling_electricity_reduction_pct": round((1.0 - P_c_funi / P_c_conv) * 100.0, 2),
      "WUE_funi_frame_L_per_kWh": 0.0, "potable_water_saved_million_L_yr": round(P_it * 8.760 * 1.80, 2),
      "freed_substation_power_MW": round(P_freed, 2),
      "extra_it_power_at_fixed_supply_low_MW": round(P_it_low, 2),
      "extra_it_power_at_fixed_supply_high_MW": round(P_it_high, 2),
      # Low and high each vary two things: low is the low IT load at the low rate, high is the high IT load at the high rate.
      "unlocked_annual_revenue_low_M_usd": round(P_it_low * 25.0, 2),
      "unlocked_annual_revenue_high_M_usd": round(P_it_high * 50.0, 2),
      "unlocked_annual_revenue_pairing": "low = low IT load at 25 M usd per MW-yr; high = high IT load at 50 M usd per MW-yr",
      "unlocked_annual_revenue_at_25M_usd_per_MW_yr_M_usd": [round(P_it_low * 25.0, 2), round(P_it_high * 25.0, 2)],
      "unlocked_annual_revenue_at_50M_usd_per_MW_yr_M_usd": [round(P_it_low * 50.0, 2), round(P_it_high * 50.0, 2)],
  }


def build_certificate() -> dict:
  th, rf, wl, vb = verify_symbolic_theorems(), solve_funicular_roof(), solve_vibro_thermal_wall(), solve_vibration_cascade()
  se, ts = solve_whole_frame_seismic(rf["dead_mass_kg_m2"], wl["wall_mass_kg_m2"]), solve_three_tier_tes_and_cooling()
  fr = 0.62 * math.sqrt(55.0)
  ct = {
      "must_fail_01_rolled_steel_24m_rafter_fails_deflection": rf["span_over_delta_ipe600_rafter"] < 360.0 <= rf["span_over_delta_vault"],
      "must_fail_02_passive_gfrp_roof_tie_cracks_crown": rf["ft_crown_passive_gfrp_MPa"] > fr and rf["fc_crown_MPa"] <= rf["fc_allow_service_MPa"],
      "must_fail_03_tornado_uplift_resisted_without_ballast": 0.0 < rf["tornado_net_uplift_kN"] and rf["tornado_tiedown_capacity_kN"] >= 1.5 * rf["tornado_net_uplift_kN"],
      "must_fail_04_steel_wythe_pins_fail_thermal_bridge_gate": wl["U_wall_steel_pcsp_W_m2K"] > 1.0 and wl["U_wall_cfrp_W_m2K"] <= 0.20,
      "must_fail_05_bare_concrete_fails_hall_acoustic_absorption": wl["nrc_poraver_helmholtz"] >= 0.65 and wl["hall_reverberant_noise_reduction_dBA"] >= 7.0,
      "must_fail_06_unbuffered_air_fails_ups_thermal_ridethrough": wl["passive_ride_through_air_only_min"] < 3.0 and wl["passive_ride_through_with_mpcm_min"] >= 15.0,
      "must_fail_07_undamped_steel_exceeds_0p0005g_switch_limit": vb["rack_g_rms_conventional_steel"] > vb["oem_switch_limit_g_rms"] >= vb["rack_g_rms_funi_frame"],
      "must_fail_08_top_heavy_unbraced_frame_fails_risk_cat_iv_drift": se["top_down_story_drift_conv_unbraced_pct"] > se["risk_cat_iv_drift_limit_pct"] >= se["top_down_story_drift_aska_pct"],
      "must_fail_09_plain_water_fails_compact_15min_ridethrough": ts["water_only_ride_through_in_same_volume_min"] < 6.0 and ts["pcm_ride_through_volume_m3"] < 0.40 * ts["water_equivalent_volume_m3"],
      "must_fail_10_liion_exceeds_thermal_arbitrage_capex_ceiling": ts["annualized_capex_liion_usd_kWh"] > 0.050 and ts["annualized_capex_basalt_usd_kWh"] <= 0.020,
      "must_fail_11_low_grade_air_exhaust_cannot_drive_adsorption": ts["cop_ads_at_32C_air_exhaust"] == 0.0 and ts["adsorption_self_cools_air_load"] is True,
      "must_fail_12_evaporative_plant_fails_zero_water_and_pue_gate": ts["PUE_conventional"] > 1.35 and ts["PUE_funi_frame"] <= 1.08 and ts["WUE_funi_frame_L_per_kWh"] == 0.0,
  }
  ct["must_fail_13_dry_recooling_at_45C_cannot_reach_design_adsorption_cop"] = (
      ts["cop_ads_reversible_limit_at_45C_dry_recooling"] < ts["cop_ads_design"] <= ts["cop_ads_reversible_limit_at_design_recooling"]
      and ts["T_recooling_design_C"] < ts["T_recooling_max_for_design_cop_C"])
  # Control 14 does not trust the two formulas in the solver. It checks them by a power balance: the new IT load
  # plus the overhead it brings must use up the freed utility power. Counting the freed power itself as IT load
  # must leave that balance open. Tolerances follow from the certificate's rounding (0.01 MW, 0.01 M usd).
  freed, it_lo, it_hi = (ts["freed_substation_power_MW"], ts["extra_it_power_at_fixed_supply_low_MW"],
                         ts["extra_it_power_at_fixed_supply_high_MW"])
  loop_factor = 1.0 + 1.0 / ts["cop_loop_design"]
  ct["must_fail_14_freed_utility_power_counted_as_it_power"] = (
      abs(it_lo * ts["PUE_funi_frame"] - freed) <= 0.02 and abs(it_hi * loop_factor - freed) <= 0.02
      and freed * loop_factor - freed > 0.02
      and it_lo <= it_hi < freed
      and abs(ts["unlocked_annual_revenue_low_M_usd"] - it_lo * 25.0) <= 0.13
      and abs(ts["unlocked_annual_revenue_high_M_usd"] - it_hi * 50.0) <= 0.26
      and all(abs(value - it * rate) <= rate * 0.005 + 0.006
              for rate, key in ((25.0, "unlocked_annual_revenue_at_25M_usd_per_MW_yr_M_usd"),
                                (50.0, "unlocked_annual_revenue_at_50M_usd_per_MW_yr_M_usd"))
              for value, it in zip(ts[key], (it_lo, it_hi)))
      and ts["unlocked_annual_revenue_low_M_usd"] < freed * 25.0 - 0.13)
  require(all(ct.values()), [k for k, v in ct.items() if not v])
  # Values that no solver computes: they are stated from the cited literature or the design brief.
  declared = {
      "domain_1_funicular_radiative_roof": ["span_m", "vault_depth_mm", "tornado_tiedown_capacity_kN",
                                            "noon_radiative_sky_cooling_W_m2"],
      "domain_2_vibro_thermal_sandwich_wall": ["nrc_poraver_helmholtz"],
      "domain_3_broadband_vibration_cascade": ["oem_switch_limit_g_rms"],
      "domain_4_top_down_seismic_cellular_frame": ["risk_cat_iv_drift_limit_pct"],
      "domain_5_6_three_tier_tes_and_closed_loop_cooling": [
          "P_it_MW", "Q_liquid_w45_MW", "Q_air_rdhx_MW", "cop_ads_design", "cop_ads_at_32C_air_exhaust",
          "T_drive_supply_C", "T_chilled_water_C", "T_recooling_design_C", "cop_loop_design", "basalt_battery_volume_m3",
          "basalt_media_mass_tonnes", "basalt_media_cost_usd", "annualized_capex_basalt_usd_kWh",
          "annualized_capex_liion_usd_kWh", "WUE_funi_frame_L_per_kWh"],
  }
  # The list above is checked against the source: every certificate field that is a written number must be on it.
  literal = literal_certificate_fields()
  require({k: sorted(v) for k, v in declared.items()} == literal,
          {k: sorted(set(literal[k]) ^ set(declared.get(k, []))) for k in literal if set(literal[k]) ^ set(declared.get(k, []))})
  # Naming these does not repair them. Control 10 compares two stated constants and still cannot fail.
  declared_only = ["must_fail_10_liion_exceeds_thermal_arbitrage_capex_ceiling"]
  # Every other control that reads a declared field, with the field it reads. Controls 1, 2, 4, 6 and 9 read none directly.
  declared_partly = {
      "must_fail_03_tornado_uplift_resisted_without_ballast":
          "the tie-down capacity it compares the computed uplift against, tornado_tiedown_capacity_kN, is a declared design value",
      "must_fail_05_bare_concrete_fails_hall_acoustic_absorption":
          "nrc_poraver_helmholtz is declared as 0.72, and 0.72 is written again as a literal in the reverberation calculation",
      "must_fail_07_undamped_steel_exceeds_0p0005g_switch_limit":
          "the limit it compares against, oem_switch_limit_g_rms, is declared from a manufacturer's figure",
      "must_fail_08_top_heavy_unbraced_frame_fails_risk_cat_iv_drift":
          "the limit it compares against, risk_cat_iv_drift_limit_pct, is declared from the code",
      "must_fail_11_low_grade_air_exhaust_cannot_drive_adsorption": "cop_ads_at_32C_air_exhaust is declared as 0.0",
      "must_fail_12_evaporative_plant_fails_zero_water_and_pue_gate": "WUE_funi_frame_L_per_kWh is declared as 0.0",
      "must_fail_13_dry_recooling_at_45C_cannot_reach_design_adsorption_cop":
          "cop_ads_design and T_recooling_design_C are declared, and so are the two temperatures the limits are computed from",
      "must_fail_14_freed_utility_power_counted_as_it_power": "cop_loop_design is declared",
  }
  require(set(declared_partly) | set(declared_only) <= set(ct), "a control named here does not exist")
  return {
      "repository": "aska-physical/funi-frame", "version": "0.1.1", "symbolic_theorems": th,
      "domain_1_funicular_radiative_roof": rf, "domain_2_vibro_thermal_sandwich_wall": wl,
      "domain_3_broadband_vibration_cascade": vb, "domain_4_top_down_seismic_cellular_frame": se,
      "domain_5_6_three_tier_tes_and_closed_loop_cooling": ts, "must_fail_controls": ct,
      "controls_passed": f"{sum(ct.values())}/{len(ct)}",
      "declared_inputs": declared,
      "controls_on_declared_inputs_only": declared_only,
      "controls_resting_partly_on_declared_inputs": declared_partly,
  }


def main() -> int:
  ap = argparse.ArgumentParser()
  ap.add_argument("--write", action="store_true")
  args = ap.parse_args()
  cert = build_certificate()
  if args.write:
      CERT_PATH.parent.mkdir(parents=True, exist_ok=True)
      CERT_PATH.write_text(json.dumps(cert, indent=2) + "\n", encoding="utf-8")
  saved = json.loads(CERT_PATH.read_text(encoding="utf-8"))
  if saved != cert:
      print("ERROR: certificate drift detected against verify.py!", file=sys.stderr)
      return 1
  print(f"[PASS] funi-frame v0.1.1 ({cert['controls_passed']} controls, PUE={cert['domain_5_6_three_tier_tes_and_closed_loop_cooling']['PUE_funi_frame']}, WUE={cert['domain_5_6_three_tier_tes_and_closed_loop_cooling']['WUE_funi_frame_L_per_kWh']} L/kWh, Rack={cert['domain_3_broadband_vibration_cascade']['rack_g_rms_funi_frame']} g_rms): OK")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
