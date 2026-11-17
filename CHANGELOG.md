# Changelog

All notable machine-verified engineering milestones for `aska-physical/funi-frame`.

## [0.1.0] — 2026-11-16

### Added
- **`verify.py` (`v0.1.0`) & `results/certificate.json`:** Six-domain thermo-structural & vibro-acoustic verification suite guarded by **4 SymPy symbolic theorems** and **12/12 must-fail controls**:
  - **Domain 1 (`24.0 m` Tied Funicular Radiative Roof):** Castigliano `Fe-SMA + T1100G` arch-tie solver (`L/f = 17.78`), `L/2,044` service stiffness vs. `IPE 600` rafter failure (`L/71.6`), tornado uplift tie-down check (`783.0 >= 307.0 kN`), and `93 W/m²` radiative sky-cooling + `U = 0.1234 W/(m²·K)` (`R-46.0`) envelope model.
  - **Domain 2 (`Carbon-TRC + Poraver + mPCM` Sandwich Shear Wall):** Parallel-path CFRP vs. steel wythe thermal-bridge `U`-value (`0.1316` vs. `2.307 W/(m²·K)`), Sabine reverberant fan-noise drop (`NRC = 0.72`, `-10.97 dBA`), `STC 60.1 dB` generator containment, and transient Biot-limited `28 °C` `mPCM` UPS thermal ride-through solver.
  - **Domain 3 (Three-Stage Broadband Vibration Isolation):** Coupled `Fe-Mn-Si` hysteretic plinth (`4.5 Hz`, `eta = 0.16`) + `±45°` CFRP `CLD` rib (`9.2 Hz`, `eta = 0.12`) + nodal collar transmissibility (`0.000090 <= 0.000500 g_rms`, `-63.7 dB` improvement over steel).
  - **Domain 4 (Top-Down 2-Story Seismic Cellular Frame):** `ASCE 7` Risk Category IV inverted-pendulum base-shear (`-44.2%`) and story-drift (`0.563% <= 1.00%`) verification under top-floor-first `313 psf` rack occupancy (`-64.8%` dead mass).
  - **Domains 5 & 6 (Three-Tier Waste-Heat TES & Closed-Loop Sorption Cooling):** Integrated sizing and `PUE`/`WUE` model for `83.3 m³` `HS24` `PCM` skids, `88.7 m³` (`20 MWh_th`) crushed-basalt batteries (`$0.014/kWh_th`), `200` thermo-active foundation piles (`250 kW_th`), `ATES`/`RTES` wellfields, and direct `58 °C` `ASHRAE W45` silica-gel adsorption chilling (`PUE = 1.074`, `WUE = 0.0 L/kWh`, `+3.06 MW` freed substation power).
- **`README.md` & `HANDOFF.md`:** Full architectural specification, negative findings capture, and `v0.2.0–v0.4.0` issue queue.
