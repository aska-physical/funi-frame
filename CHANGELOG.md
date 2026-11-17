# Changelog

All notable machine-verified engineering milestones for `aska-physical/funi-frame`.

## [0.1.1] — 2026-11-17

Corrections from the cross-repository review in [`aska-physical/optifab`](/aska-physical/optifab) (`HANDOFF.md` §3). Prepared by an AI agent; awaiting human sign-off.

### Changed
- **Unlocked compute is reported in IT power, not utility power.** The `3.06 MW` of freed substation power carries **`2.85–2.96 MW`** of new IT load (`P_freed / PUE` if every overhead scales with load; `P_freed / (1 + 1/COP_loop)` if only the `ATES` loop does). `unlocked_annual_revenue_low_M_usd` falls from `76.62` to **`71.35`**, and `unlocked_annual_revenue_high_M_usd` (**`147.95`**) is now in the certificate. `README.md` §6 retitled from "Freed Compute Power" to "Freed Utility Power". This assumes the `$25–$50M/MW-yr` value is per MW of IT load.

### Added
- **Re-cooling temperature and its second-law limit (`Theorem 5`, `Must-fail #13`).** The certificate now carries `T_drive_supply_C = 58`, `T_chilled_water_C = 18`, `T_recooling_design_C = 30`, the reversible COP limit at design re-cooling (`2.051`) and at `45 °C` dry re-cooling (`0.423`), and `T_recooling_max_for_design_cop_C = 41.46`. The `30 °C` design value is a choice made in this version and needs the owner's confirmation.
- **`Must-fail #14`:** freed utility power may not be counted as IT power.
- **`declared_inputs` and `controls_on_declared_inputs_only`** in the certificate: every value that is stated from the literature or the design brief rather than computed, and the one control (`#10`) that compares only such values.

### Not changed
- No structural, acoustic, vibration, seismic or storage result. `PUE`, `WUE` and all `12` original controls are as in `v0.1.0`.

## [0.1.0] — 2026-11-16

### Added
- **`verify.py` (`v0.1.0`) & `results/certificate.json`:** Six-domain thermo-structural & vibro-acoustic verification suite guarded by **4 SymPy symbolic theorems** and **12/12 must-fail controls**:
  - **Domain 1 (`24.0 m` Tied Funicular Radiative Roof):** Castigliano `Fe-SMA + T1100G` arch-tie solver (`L/f = 17.78`), `L/2,044` service stiffness vs. `IPE 600` rafter failure (`L/71.6`), tornado uplift tie-down check (`783.0 >= 307.0 kN`), and `93 W/m²` radiative sky-cooling + `U = 0.1234 W/(m²·K)` (`R-46.0`) envelope model.
  - **Domain 2 (`Carbon-TRC + Poraver + mPCM` Sandwich Shear Wall):** Parallel-path CFRP vs. steel wythe thermal-bridge `U`-value (`0.1316` vs. `2.307 W/(m²·K)`), Sabine reverberant fan-noise drop (`NRC = 0.72`, `-10.97 dBA`), `STC 60.1 dB` generator containment, and transient Biot-limited `28 °C` `mPCM` UPS thermal ride-through solver.
  - **Domain 3 (Three-Stage Broadband Vibration Isolation):** Coupled `Fe-Mn-Si` hysteretic plinth (`4.5 Hz`, `eta = 0.16`) + `±45°` CFRP `CLD` rib (`9.2 Hz`, `eta = 0.12`) + nodal collar transmissibility (`0.000090 <= 0.000500 g_rms`, `-63.7 dB` improvement over steel).
  - **Domain 4 (Top-Down 2-Story Seismic Cellular Frame):** `ASCE 7` Risk Category IV inverted-pendulum base-shear (`-44.2%`) and story-drift (`0.563% <= 1.00%`) verification under top-floor-first `313 psf` rack occupancy (`-64.8%` dead mass).
  - **Domains 5 & 6 (Three-Tier Waste-Heat TES & Closed-Loop Sorption Cooling):** Integrated sizing and `PUE`/`WUE` model for `83.3 m³` `HS24` `PCM` skids, `88.7 m³` (`20 MWh_th`) crushed-basalt batteries (`$0.014/kWh_th`), `200` thermo-active foundation piles (`250 kW_th`), `ATES`/`RTES` wellfields, and direct `58 °C` `ASHRAE W45` silica-gel adsorption chilling (`PUE = 1.074`, `WUE = 0.0 L/kWh`, `+3.06 MW` freed substation power).
- **`README.md` & `HANDOFF.md`:** Full architectural specification, negative findings capture, and `v0.2.0–v0.4.0` issue queue.
