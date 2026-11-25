# Changelog

All notable machine-verified engineering milestones for `aska-physical/funi-frame`.

## [0.1.2] — 2026-11-25

Physical hardening of roof virtual-work stiffness, wall/skid PCM heat-transfer kinetics, and state-dependent adsorption COP.

### Fixed
- **Exact parabolic arch-tie virtual-work deflection (`solve_funicular_roof`, `SymPy Theorem 7`, `Must-fails #1, #15`):** Replaced the legacy `0.18` beam factor (`11.74 mm`, `L/2,043.6`) with exact curved parabolic compatibility and two-system Principle of Virtual Forces invariance (`< 1e-6 mm`), including `2.15 kPa` dead load, arch axial shortening, and tie elongation (`H_service = 2,178.90 kN`, `f_c = 10.75 <= 24.75 MPa`). Net total-load crown sag is **`34.96 mm` (`L/686.4`) only after `-93.21 mm` upward `dead + Fe-SMA prestress` camber** (`-117.19 mm` prestress-only, `+23.97 mm` dead-load take-up, `4.05 mm` net service tie elongation). Explicitly flags live-load `L/360` failure: post-activation incremental `11.5 kPa` live/plant load stretches the tie **`31.90 mm`** for **`128.18 mm` (`L/187.2 > L/360` and `> L/240`, `passes_live_load_L360 = False`)** downward excursion (`152.15 mm` / `37.86 mm` tie elongation un-prestressed).
- **Film-convection-limited `mPCM` wall ride-through ODE (`solve_vibro_thermal_wall`, `SymPy Theorem 8`, `Must-fail #6`):** Added exact lumped air-to-liner convection (`q = 143.82 W/m²` across `5,215 m²`). Passive natural convection (`h = 8.0 W/(m²·K)`) yields equilibrium `45.98 °C >> 35 °C`, hitting the `35 °C` ceiling in **`2.08 min`** (`1.14%` sensible+latent utilization vs. `28.30 min` idealized bucket); unlocking the bucket requires forced aisle airflow **`h >= 20.55 W/(m²·K)`**.
- **Finite-NTU isothermal `HS24` PCM skid sizing (`solve_three_tier_tes_and_cooling`, `Must-fails #9, #16`):** Modeled `eps = 1 - exp(-NTU)` (`0.7769` at `NTU = 1.5`) as a heat-transfer rate limit rather than energy loss. Counterfactual nominal `1.0x` loop flow is capped at **`7.769 MW < 10.0 MW`** (`2.231 MW` deficit, `19.31 min` partial discharge, fails full-rate gate); sizing secondary loop flow capacity to **`pcm_skid_flow_capacity_factor = 1.30`** (`> 1.2872x` minimum) delivers **`10.099 MW >= 10.0 MW`** for a verified **`15.00 min`** full-load ride-through with monotonically positive gates (raising `NTU` or flow factor passes every control).
- **State-dependent sorption $\text{COP}_{\text{ads}}(T_a)$ (`Must-fail #16`):** Anchoring $\eta_{\text{II}} = 0.3022$ at `30 °C` makes `COP_ads` drop dynamically to **`0.2174`** at `40 °C` (`alpha* = 82.14%`, `Q_drive = 6.900 MW`) and **`0.1279`** at `45 °C` (`Q_drive = 11.725 MW > 8.5 MW`, self-cooling fails above **`42.04 °C`**).
- Added stdlib zero-bare-assert regression suite [`tests/test_physical_hardening.py`](tests/test_physical_hardening.py), updated `declared_inputs` / `declared_partly` (16/16 controls), and regenerated `results/certificate.json` (`v0.1.2`).

## [0.1.1] — 2026-11-17

Corrections from the cross-repository review in [`aska-physical/optifab`](/aska-physical/optifab) (`HANDOFF.md` §3). Prepared by an AI agent; awaiting human sign-off.

### Changed
- **Unlocked compute is reported in IT power, not utility power.** The `3.06 MW` of freed substation power carries **`2.85–2.96 MW`** of new IT load (`P_freed / PUE` if every overhead scales with load; `P_freed / (1 + 1/COP_loop)` if only the `ATES` loop does). `unlocked_annual_revenue_low_M_usd` falls from `76.62` to **`71.35`**, and `unlocked_annual_revenue_high_M_usd` (**`147.95`**) is now in the certificate. The low figure pairs the low IT load with `$25M/MW-yr` and the high figure pairs the high IT load with `$50M/MW-yr`; `unlocked_annual_revenue_pairing` says so, and two new fields give the range at each rate alone (`71.35–73.98` and `142.69–147.95`). `README.md` §6 retitled from "Freed Compute Power" to "Freed Utility Power". This assumes the `$25–$50M/MW-yr` value is per MW of IT load, which the owner has still to confirm.

### Added
- **Re-cooling temperature and its second-law limit (`Theorem 5`, `Must-fail #13`).** The certificate now carries `T_drive_supply_C = 58`, `T_chilled_water_C = 18`, `T_recooling_design_C = 30`, the reversible COP limit at design re-cooling (`2.051`) and at `45 °C` dry re-cooling (`0.423`), and `T_recooling_max_for_design_cop_C = 41.46`. The `30 °C` design value is a proposal made in this version with no source cited. `T_recooling_design_basis` records that it awaits the owner's confirmation.
- **`Must-fail #14`:** freed utility power may not be counted as IT power. It checks both IT figures by a power balance (new IT load plus its overhead equals the freed utility power, to within rounding) and both revenue figures against them, so a wrong formula that kept the ordering fails. It checks the two single-rate revenue ranges the same way.
- **`declared_inputs`, `controls_on_declared_inputs_only` and `controls_resting_partly_on_declared_inputs`** in the certificate: every value that is stated from the literature or the design brief rather than computed, the one control (`#10`) that compares only such values, and the eight others (`#3`, `#5`, `#7`, `#8`, `#11`, `#12`, `#13`, `#14`) that read one, each with the value it reads. `declared_inputs` is checked against the source of `verify.py`, which found two it had missed (`span_m`, `vault_depth_mm`). This is disclosure, not a fix: `#10` still cannot fail, and a cost model for it is open work.

### Fixed
- **`python3 -O verify.py` skipped every check.** The verifier's checks were `assert` statements, which Python strips under `-O`: with a failing control it built a certificate reading `13/14` and raised nothing. All `10` are now calls to `require()`, which raises `RuntimeError` whatever the flags. Seven of them were in `v0.1.0`. The same hazard is tracked in `aska-physical/cfrp-masts#3` and was found again in the review of `aska-physical/optifab#13`.

- `tests/test_verify_o_mode.py` (three tests, standard library only) pins the fix: no `assert` statement in `verify.py`, a clean tree passes with and without `-O`, and a failing control stops the run with and without `-O`.
- The `[PASS]` line reads `WUE` from the certificate. It was a written `0.0`.

### Not changed
- No structural, acoustic, vibration, seismic or storage result. `PUE`, `WUE` and all `12` original controls are as in `v0.1.0`.

## [0.1.0] — 2026-11-16

### Added
- **`verify.py` (`v0.1.0`) & `results/certificate.json`:** Six-domain thermo-structural & vibro-acoustic verification suite guarded by **4 SymPy symbolic theorems** and **12/12 must-fail controls**:
  - **Domain 1 (`24.0 m` Tied Funicular Radiative Roof):** Castigliano `Fe-SMA + T1100G` arch-tie solver (`L/f = 17.78`), legacy `0.18` beam-factor `L/2,044` estimate (superseded in `v0.1.2` by exact virtual-work `34.96 mm` net of camber / `128.18 mm` `L/187.2` live-load failure) vs. `IPE 600` rafter failure (`L/71.6`), tornado uplift tie-down check (`783.0 >= 307.0 kN`), and `93 W/m²` radiative sky-cooling + `U = 0.1234 W/(m²·K)` (`R-46.0`) envelope model.
  - **Domain 2 (`Carbon-TRC + Poraver + mPCM` Sandwich Shear Wall):** Parallel-path CFRP vs. steel wythe thermal-bridge `U`-value (`0.1316` vs. `2.307 W/(m²·K)`), Sabine reverberant fan-noise drop (`NRC = 0.72`, `-10.97 dBA`), `STC 60.1 dB` generator containment, and transient Biot-limited `28 °C` `mPCM` UPS thermal ride-through solver.
  - **Domain 3 (Three-Stage Broadband Vibration Isolation):** Coupled `Fe-Mn-Si` hysteretic plinth (`4.5 Hz`, `eta = 0.16`) + `±45°` CFRP `CLD` rib (`9.2 Hz`, `eta = 0.12`) + nodal collar transmissibility (`0.000090 <= 0.000500 g_rms`, `-63.7 dB` improvement over steel).
  - **Domain 4 (Top-Down 2-Story Seismic Cellular Frame):** `ASCE 7` Risk Category IV inverted-pendulum base-shear (`-44.2%`) and story-drift (`0.563% <= 1.00%`) verification under top-floor-first `313 psf` rack occupancy (`-64.8%` dead mass).
  - **Domains 5 & 6 (Three-Tier Waste-Heat TES & Closed-Loop Sorption Cooling):** Integrated sizing and `PUE`/`WUE` model for `83.3 m³` `HS24` `PCM` skids, `88.7 m³` (`20 MWh_th`) crushed-basalt batteries (`$0.014/kWh_th`), `200` thermo-active foundation piles (`250 kW_th`), `ATES`/`RTES` wellfields, and direct `58 °C` `ASHRAE W45` silica-gel adsorption chilling (`PUE = 1.074`, `WUE = 0.0 L/kWh`, `+3.06 MW` freed substation power).
- **`README.md` & `HANDOFF.md`:** Full architectural specification, negative findings capture, and `v0.2.0–v0.4.0` issue queue.
