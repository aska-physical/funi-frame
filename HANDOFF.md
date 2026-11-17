# Inter-Agent Handoff — `aska-physical/funi-frame` (`v0.1.1`)

Prepared per the Aska group ethos (*maximal autonomy and agency; capture every computation, decision, and negative finding so downstream agents build forward without re-deriving*).

## 1. What Was Learned, Proved, and Verified (`v0.1.0`–`v0.1.1`)

- **Theorem 1 — Hysteretic (`Fe-Mn-Si` / `CLD`) vs. Viscous High-Frequency Isolation Roll-Off (`SymPy + Must-fail #7`):**
  - Classic Kelvin-Voigt viscous dampers suffer high-frequency damper lock-up: $\lim_{r\to\infty} r\,T_{\text{visc}}(r,\zeta) = 2\zeta > 0$ (`-20 dB/decade`), leaking `29.5–120+ Hz` chiller, pump, and genset harmonics into server cabinets (`0.137 g_rms` on conventional steel frames). Rate-independent martensitic (`Fe-Mn-Si`, `4.5 Hz`, $\eta_1 = 0.16$) and viscoelastic `±45°` CFRP Constrained-Layer Damping (`CLD`, `9.2 Hz`, $\eta_2 = 0.12$) follow $\lim_{r\to\infty} r^2 T_{\text{hyst}}(r,\eta) = \sqrt{1+\eta^2}$ (`-40 dB/decade`), holding white-space rack vibration to **`0.000090 g_rms <= 0.000500 g_rms`** (`17.91%` of `STRUCTURE Magazine` optical/director-switch ceiling, `-63.69 dB` improvement).
- **Theorem 2 — Wythe Thermal-Bridge Amplification (`SymPy + Must-fail #4`):**
  - $U_{\text{eff}}/U_{\text{core}} = 1 + \rho_{\text{pin}}(k_{\text{pin}}/k_{\text{ins}} - 1)$. Even `0.35%` steel wythe pins (`k = 50 W/m·K`) across `38 mm` foam degrade precast sandwich panels to `U = 2.307 W/(m²·K)`; `45°` CFRP grid connectors (`0.08%` area, transverse `k = 0.80 W/m·K`) preserve **`U = 0.1316 W/(m²·K)` (`R-43.2`)**.
- **Theorem 3 — Latent `PCM` Volumetric Gain & Transient Biot Wall Buffering (`SymPy + Must-fails #6, #9`):**
  - Plain chilled water over a `9 K` delta-T stores only `11.5 kWh_th/m³`, exhausting an `83.33 m³` mechanical buffer in `5.75 min`. Encapsulated `savE HS24` salt hydrate (`30.0 kWh_th/m³`, Osterman et al. 2023) delivers a full **`15.0 min` (`2.5 MWh_th`)** ride-through in `-61.67%` less tank volume, while the outer `~5.2 mm` active Biot penetration skin (`eta_biot = 0.10`) of `28 °C` `mPCM`-doped Poraver wall/ceiling liners passively buffers cold-aisle air `<= 35 °C` (`ASHRAE A2`) for **`28.30 min`** vs. `1.78 min` for unbuffered room air.
- **Theorem 4 — Waste-Heat Sorption Self-Cooling Threshold (`SymPy + Must-fail #11`):**
  - Critical liquid-cooling capture fraction $\alpha_{\text{liq}}^* = 1/(1 + \text{COP}_{\text{ads}}) = 61.73\%$ at silica-gel $\text{COP}_{\text{ads}} = 0.62$. Capturing `85%` of IT heat at `58 °C` (`ASHRAE W45`) self-generates `100%` of the `15–18 °C` chilled water needed for the `1.50 MW_th` air/RDHx load using `2.419 MW_th` of server reject heat with **zero mechanical chiller compressors**, dropping campus `PUE` from `1.38` to **`1.074`** at **`WUE = 0.000 L/kWh`**.
- **Theorem 5 — Sorption Re-Cooling Limit (`SymPy + Must-fail #13`, added in `v0.1.1`):**
  - A heat-driven chiller between a drive at $T_h$, re-cooling at $T_a$ and chilled water at $T_c$ obeys $\text{COP}_{\text{ads}} \le \frac{T_c}{T_a - T_c}\cdot\frac{T_h - T_a}{T_h}$, so $T_a \le \frac{T_c T_h (1 + \text{COP}_{\text{ads}})}{\text{COP}_{\text{ads}} T_h + T_c}$. At `58 °C` drive, `18 °C` chilled water and `COP_ads = 0.62` the re-cooling temperature must stay below **`41.46 °C`**. The design value is `30 °C` (reversible limit `2.051`, so the machine runs at `30%` of the reversible limit). The certificate now carries the three temperatures.
- **Identity 6 — IT Power From Freed Utility Power (`SymPy + Must-fail #14`, added in `v0.1.1`):**
  - Freed utility power $P_{\text{freed}}$ carries $P_{\text{freed}}/\text{PUE}$ of new IT load if every overhead scales with load and $P_{\text{freed}}/(1 + 1/\text{COP}_{\text{loop}})$ if only the `ATES` loop does: **`2.85–2.96 MW`** for the `3.06 MW` freed here. `v0.1.0` multiplied the full `3.06 MW` by the compute value per MW; `v0.1.1` uses the IT figure (`$71.35M–$147.95M/yr`). The low end pairs the low IT load with `$25M/MW-yr` and the high end pairs the high IT load with `$50M/MW-yr`; at one rate the range is `$71.35M–$73.98M` (`$25M`) or `$142.69M–$147.95M` (`$50M`). This assumes the `$25–$50M/MW-yr` value is per MW of IT load, which the owner has still to confirm.
  - `Must-fail #14` checks both figures by a power balance: the new IT load plus the overhead it brings must use up the freed utility power to within the certificate's rounding, and counting the freed power itself as IT load must leave the balance open. The two formulas are also proved as SymPy identities.
- **Negative Findings Captured (Do Not Re-Derive):**
  1. **Shallow `L/f > 30` (`f = 0.68 m`) on `24 m` heavy-plant roof vaults over-stresses the crown:** Keeping `L/f <= 18` (`f = 1.35 m`, `L/f = 17.78`, `h_vault = 1,480 mm`) cuts horizontal arch thrust in half (`2,181.65 kN`, `f_c = 9.60 <= 24.75 MPa`) while still saving `920 mm` against `2,400 mm` deep heavy roof trusses.
  2. **Stage-2 rib isolation above `15 Hz` amplifies `29.5 Hz` screw chillers past `0.0005 g_rms`:** Both isolation stages must sit below `f_min / 3` (`4.5 Hz` plinth + `9.2 Hz` `CLD` cassette seat) to prevent `29.5 Hz` chiller bleed-through.
  3. **Air-cooled server exhaust (`32 °C`) and `ASHRAE W17` (`25.3 °C`) cannot drive sorption chillers without heat-pump lift (`Must-fail #11`):** Direct heat-pump-free regeneration requires `>= 55 °C` (`W45`/`W+`) for silica-gel or `>= 70–75 °C` (`HTHP`-boosted, ORNL OSTI 3740567) for single-effect `LiBr-H2O`.
  4. **Li-ion battery arbitrage (`$0.068/kWh` annualized pack capex) is `4.86x` costlier than crushed basalt (`$0.014/kWh_th`, `Must-fail #10`).** Both figures are declared inputs: `Must-fail #10` compares two stated constants and cannot fail (listed under `controls_on_declared_inputs_only`). Listing it is disclosure, not a fix. A cost model that derives the two figures is open work, and until it exists `#10` guards nothing. `Must-fail #5`, `#11` and `#12` rest partly on declared inputs (`0.72` absorption, written twice in `verify.py`; `COP = 0.0` at `32 °C` exhaust; `WUE = 0.0`) and are listed under `controls_resting_partly_on_declared_inputs`.
  5. **`45 °C` dry re-cooling cannot reach the `0.62` adsorption COP (`Must-fail #13`):** the reversible limit there is `0.423`. Zero-water operation at this COP needs the `ATES` / energy-pile / sky sinks.

## 2. Repository Artifacts (`main`)

- [`verify.py`](verify.py) (`v0.1.1`) — 6 SymPy identities + 6 coupled thermo-structural solvers + **14/14 must-fail controls**.
- [`results/certificate.json`](results/certificate.json) (`v0.1.1`) — Canonical machine-checked outputs, with `declared_inputs` naming every value that is stated rather than computed.
- [`docs/ENGINEERING_STUDIES.md`](docs/ENGINEERING_STUDIES.md) — Complete archival monograph of the Whole-Building Vibro-Thermal Shell Study and the Three-Tier Waste-Heat TES & Closed-Loop Sorption Cooling Study.
- [`README.md`](README.md) & [`CHANGELOG.md`](CHANGELOG.md) — System specification and version history.

## 3. Next-Agent Build Queue (`v0.2.0` – `v0.4.0`)

1. **Issue #1 — Wavenumber-dependent Ross-Kerwin-Ungar (`RKU`) 4th-order `CLD` solver & acoustic transfer-matrix `STC`/`NRC` verification.**
2. **Issue #2 — Hourly `8,760-h` hydro-thermal aquifer plume balance (`HT-ATES` / brackish `RTES`) & ERCOT negative-price basalt-battery dispatch.**
3. **Issue #3 — Modular wall-to-vault corner collar joint, `Fe-SMA` ring-spring plinth FEA parameterization & EF-3 tornado projectile/blast shear gate.**
