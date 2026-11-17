# Engineering Studies Monograph — `aska-physical/funi-frame` (`v0.1.0`)

Complete archival capture of both session engineering studies for [`aska-physical/funi-frame`](/aska-physical/funi-frame), machine-verified by [`verify.py`](../verify.py) against [`results/certificate.json`](../results/certificate.json) (`12/12` must-fail controls):
1. **Study I:** *Whole-Building Monolithic Vibro-Thermal Data-Center Exoskeleton (Funicular Radiative Roof, Zero-Bridge Poraver-mPCM Sandwich Shear Walls, Three-Stage Hysteretic Vibration Cascade & Top-Down Seismic Cellular Box)*
2. **Study II:** *Capturing AI Data-Center Waste Heat, Three-Tier High-Density Thermal Energy Storage (PCM, Mineral Packed Beds & Subsurface ATES/RTES/Energy Piles), and Closed-Loop Sorption Cooling*

## Study I — Whole-Building Monolithic Vibro-Thermal Data-Center Exoskeleton

### I.1 Four Failure Modes of Conventional Precast + Steel Truss Data-Center Shells
Modern `80–250 kW/rack` AI halls (`12–15+ kN/m²` live loads, `$200,000+` per rack) break conventional tilt-up concrete + steel Pratt truss shells (*STRUCTURE Magazine* Linwood 2023; *CalcSteel* 2026):
1. **Span vs. Plenum Depth & Tornado Uplift:** Across `24.0 m` column-free bays carrying `7.5 kPa` rooftop plant + `4.0 kPa` hung busways/manifolds (`~15 kN/m`), rolled `IPE 600` rafters fail deflection at `L/71.6` (`< L/360`), forcing `2.4 m`-deep trusses (`460 kg/m²` ballasted against `6.2 kPa` EF-2 tornado suction).
2. **Wythe Thermal Bridging & Reverberant Noise:** `280 mm` precast sandwich panels (`368 kg/m²`) with `0.35%` steel wythe pins (`k = 50 W/m·K`) across `38 mm` foam degrade to `U = 2.307 W/(m²·K)` (`R-2.5`), while bare concrete (`NRC = 0.03`) reflects `85–95+ dBA` fan noise.
3. **High-Frequency Vibration Lock-Up:** Steel frames (`eta ~ 0.003`) + `18 Hz` pads leak `29.5–120 Hz` chiller/pump/genset tones into racks at **`0.137 g_rms`** (`274x` over the **`0.0005 g_rms`** `10–200 Hz` optical/director-switch limit).
4. **Top-Down Seismic Inverted Pendulum & Dust:** Upper-floor-first tenant phasing creates a top-heavy inverted pendulum (`4.359% > 1.00%` `ASCE 7` Risk Cat IV drift on unbraced heavy frames), while exposed concrete soffits shed alkaline dust.

### I.2 Domains 1–4 Architectural Layup & Verified Performance (`10 MW`, `2,500 m²` 2-Story Campus)
- **Domain 1 — `24.0 m` Tied Double-Skin Funicular Radiative Roof (`h = 1,480 mm`, `f = 1,350 mm`, `L/f = 17.78`):**
  - *Layup:* Photonic/fluoropolymer hydronic sky emitter (`96%` solar reflectance, `94.6%` `8–13 µm` emittance, **`93.0 W/m²`** noon cooling, `4.9 °C` sub-ambient) + `70 mm` `55 MPa` Poraver-TRC compression vault (`4` ribs `75x280 mm`) + `150 mm` PIR/aerogel (`k = 0.022 W/m·K`) + `4` [`cfrp-truss`](/aska-physical/cfrp-truss) (`Af = 780 mm²`) & `175 °C` Joule-activated `Fe-SMA` (`Asma = 1,150 mm²`, `423.5 MPa` recovery) bottom ties.
  - *Verified (`Must-fails #1–#3`):* Mass **`219.32 kg/m²`** (`-52.32%`), saves **`920 mm`** plenum height, `H = 2,181.65 kN/bay`, crown **`fc = 9.60 <= 24.75 MPa`** (passive `GFRP` tie cracks at `ft = 266.03 MPa`), sag **`11.74 mm` (`L/2,043.6`)**, tornado tie-down **`783.0 >= 307.03 kN`**, **`U = 0.1234 W/(m²·K)` (`R-46.0`)**.
  - **Domain 2 — Zero-Bridge `Carbon-TRC + Poraver + mPCM` Sandwich Shear Wall (`260 mm`, `159.66 kg/m²`, `-56.61%`):**
    - *SymPy Theorem 2:* $U_{\text{eff}}/U_{\text{core}} = 1 + \rho_{\text{pin}}(k_{\text{pin}}/k_{\text{ins}} - 1)$. Using `45°` CFRP grids (`rho = 0.08%`, `k = 0.80 W/m·K`) across `160 mm` foam between `48 mm` exterior Carbon-TRC (`STC 60.1 dB`) and `52 mm` interior Helmholtz-slotted Poraver expanded glass (`lambda = 0.16 W/m·K`, **`NRC = 0.72`**, **`-10.97 dBA`** hall noise drop) yields **`U = 0.1316 W/(m²·K)` (`R-43.2`)** (`Must-fails #4–#5`).
    - *Transient Biot Ride-Through (`Must-fail #6`):* Doping the interior Poraver liner with `28 wt%` `28 °C` micro-encapsulated `PCM` (`~5.2 mm` active Fourier skin, `eta_biot = 0.10`) extends cold-aisle passive ride-through (`<= 35 °C` under `750 kW` heat leak) from **`1.78 min` to `28.30 min`**.
    - **Domain 3 — Three-Stage Hysteretic Vibration Cascade (`SymPy Theorem 1` & `Must-fail #7`):**
      - *SymPy Proof:* $\lim_{r\to\infty} r\,T_{\text{hyst}} = 0$ and $\lim_{r\to\infty} r^2 T_{\text{hyst}} = \sqrt{1+\eta^2}$ (`-40 dB/dec`) vs. viscous $\lim_{r\to\infty} r\,T_{\text{visc}} = 2\zeta$ (`-20 dB/dec` lock-up).
      - *Cascade (`29.5/49/60/120 Hz` plant tones):* `Fe-Mn-Si` ring-spring plinth (`4.5 Hz`, `eta=0.16`) + `±45°` CFRP `CLD` cassette seat (`9.2 Hz`, `eta=0.12`) + `0.40` nodal collar break -> **`0.000090 <= 0.000500 g_rms`** (`17.91%` utilization, **`-63.69 dB`** vs. steel).
      - **Domain 4 — Top-Down Seismic Cellular Box (`Must-fail #8`):**
        - Total height **`12.38 m`** (`-2.02 m` vs. `14.40 m`), shell dead mass **`1,431.4 t`** vs. `4,069.8 t` (**`-64.83%`**), `ASCE 7` Risk Cat IV base shear **`10,118.9 kN`** (**`-44.23%`**), top-down story drift **`0.563% <= 1.00%`** (unbraced heavy frame `4.359%`), and **zero concrete spalling dust** (`SIP` composite soffits).

        ## Study II — AI Data-Center Waste-Heat Capture, Three-Tier TES & Closed-Loop Sorption Cooling

        ### II.1 Server Heat Capture Regimes (`180 TWh` 2024 -> `325–580 TWh` 2028; `$25–$50M/MW-yr` Compute Value)

        | Capture Architecture | Rack Power | Liquid Share | Return Temp. | Thermodynamic Role (`Bytes2Heat` 2024; `Sandia/LBNL` 2023; `nVent CDU`) |
        | :--- | :---: | :---: | :---: | :--- |
        | `CRAC / CRAH` Air (`ASHRAE T1`) | `< 12 kW_th` | `0%` | `20–40 °C` | Low exergy; `COP_ads = 0.0` at `32 °C` (`Must-fail #11`) |
        | Rear-Door HX (`RDHx`) | `15–55 kW_th` | `50–70%` | `20–30 °C` | Ambient loop / `15%` (`1.5 MW_th`) air trim in `funi-frame` |
        | Direct Cold Plates (`W17–W45+`) | `80–200 kW_th` | `60–90%` (`85%`) | `30–60 °C` (`58 °C`) | Drives silica-gel adsorption (`>=55 °C`), `HT-ATES` & `HTHP` |
        | Single / Two-Phase Immersion | `<= 250 kW_th` | `Up to 100%` | `30 to >60 °C` | Fan-free isothermal capture; highest Carnot availability |

        ### II.2 Three-Tier Thermal Energy Storage (`TES`) Sizing & Benchmarks

        | TES Tier & Modality | Operating Band | Vessel Density | Cost & Sizing on `10 MW` (`87.6 GWh_th/yr`) AI Hall |
        | :--- | :---: | :---: | :--- |
        | **Tier 1: `savE HS24` Salt-Hydrate `PCM`** (`40 vol%`) | `23–32 °C` (`9 K`) | **`30.0 kWh_th/m³`** | **`83.33 m³`** (`15 min` / `2.5 MWh_th`) vs. `217.39 m³` water (`11.5 kWh_th/m³`, `5.75 min` in `83.33 m³`, `Must-fail #9`; `4–5 yr` payback) |
        | **Tier 2a: Crushed Basalt Bed** (`2 x 10 MWh_th`) | `100–600 °C` (`500 K`) | **`225.5 kWh_th/m³`** | **`88.68 m³`** (`117.06 t`, `$14,046` media), **`$0.014/kWh_th`** annualized (`4.86x < $0.068/kWh` Li-ion, `Must-fail #10`) |
        | **Tier 2b: Quartzite Silica Sand** (`10 MWh_th` alt.) | `100–600 °C` (`500 K`) | `144.1 kWh_th/m³` | `69.39 m³` (`86,737 kg`), **`$3,470` media/bed** (`$0.35/kWh_th`), `$0.017/kWh_th` annualized |
        | **Tier 2c: Magnetite Ore Bed** (`10 MWh_th` alt.) | `100–600 °C` (`500 K`) | **`300.2 kWh_th/m³`** | **`33.31 m³`** (`84,597 kg`, most compact vessel), `$42,298` media, `$0.018/kWh_th` annualized |
        | **Tier 3: `ATES` / Brackish `RTES` + `200` Energy Piles** | `10–65 °C` | Formation scale | **`250 kW_th`** piles (`25 m x 50 W/m`) + `$2–$9/MWh_th` `ATES` (`COP >10–25`) / deep brackish `RTES` (`COP 38–54`, `-78%` power) |

        ### II.3 Closed-Loop Sorption Cooling, Structural Conduit Coupling & `10 MW` Campus Balance
        - **SymPy Theorem 4 (`Must-fail #11`):** Self-cooling liquid-capture threshold $\alpha_{\text{liq}}^* = 1/(1 + \text{COP}_{\text{ads}}) = 61.73\%$ at silica-gel `COP_ads = 0.62` (`55–100 °C` drive). Capturing `85%` (`8.50 MW_th`) at `58 °C` (`W45`) feeds **`2.419 MW_th`** directly into silica-gel adsorption chillers to produce `1.50 MW_th` of `15–18 °C` `RDHx` chilled water with **zero compressors and zero heat-pump lift** (optional `R744/R717` `HTHP` at `COP 3.5–6.0` lifts `58 °C` -> `75–95 °C` for `20.5 MW` `LiBr-H2O` chillers or `5–8 MW_th` `4GDH` district export).
        - **Structural Manifold & `50 kHz` EIS Leak Twin:** Insulated `TCS/FWS` headers route inside `200 mm` `45°` Torayca T1100G rib apertures of [`funi-floors`](/aska-physical/funi-floors) with [`impedence-analyzers`](/aska-physical/impedence-analyzers) `50 kHz` Cole-Cole moisture sensing (`+46.16%` reactance shift, `0.0000%` load crosstalk).
        - **Verified `10 MW` Campus Result (`Must-fail #12`):** Cuts `PUE` from `1.380` to **`1.074`** (**`-85.05%` cooling electricity**), achieves **`WUE = 0.000 L/kWh`** (**`157.68M L/yr` potable water saved**), and frees **`3.06 MW` of substation power, which carries `2.85–2.96 MW` of new IT load (`$71.35M–$147.95M/yr` unlocked AI compute revenue; corrected in `v0.1.1`, see [`CHANGELOG.md`](../CHANGELOG.md))**.

        ## Primary Engineering References
        1. **Linwood, E.** (2023), *STRUCTURE Magazine* (*Design Parameters for Data Center Facilities*); **CalcSteel** (2026) (*Steel for AI Data Centers*); **S. A. Miro** (2026) & **ISTA Engineers**.
        2. **O'Hegarty & Kinnane** (2022) (*PCSP Review*, `ResearchGate 364675724`); *Buildings* `12(9):1429` (`DOI: 10.3390/buildings12091429`); *CarbonCast*; *Expanded Glass* (`ResearchGate 368676656`, `358829471`).
        3. **Udovenko et al.** (*High-Damping Fe-Mn-Si Alloys*, `25–30%` decrement); **Uni Augsburg** (`OPUS 122830`, `±45°` CFRP `CLD`); **Radiative Sky Cooling** (*ACS AMI* `16(25):32679`; *EnerarXiv* `5084`; `PMC12986632`; *CEC* `CEC-500-2025-002` *SkyCool*).
        4. **Pandey, Stumpf & Lin** (2026), *Groundwater* `64(4):394–404` (`PMC13353728`, `ATES`/`RTES`); **Bhattacharyya** (2023), *IJREC* `9(2):14–21` (Packed-Bed Thermal Batteries); **Osterman et al.** (2023), *Energies* `16(17):6127` (`savE HS24` `PCM`).
        5. **Ott & Radgen** (2024), Univ. Stuttgart *Bytes2Heat*; **Wang et al.** (2026), ORNL `OSTI 3740567` (*HTHP + LiBr-H2O*); **Sandia/LBNL** (2023); **nVent** `SB-H85398`; **Calentix**; **Johnson Controls** (2026); **Brandl** (`ResearchGate 257726802`, Energy Piles); **Gupta & Puri** (`ResearchGate 353463282`); **Composite Moisture EIS** (`ResearchGate 321206389`).
