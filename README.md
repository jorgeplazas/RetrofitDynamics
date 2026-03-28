# RetrofitDynamics
Code and data for the article *Integrated Dynamic Modeling of $\mathrm{CO}_{2}$ Emission in Aviation Retrofit Strategies* by Francisco A. Buendía Hernández and Jorge Plazas.

This repository packages the MATLAB and R code used to (i) run the retrofit **scenario simulations**
and (ii) perform the **time‑resolved global sensitivity analysis** reported in the manuscript.

> **Reproducibility first.** This README provides environment requirements, a file inventory,
> step‑by‑step execution, expected inputs/outputs, and a checklist of items to provide if any
> file is missing.

---

## 1. Repository layout (proposed)

```text
code/
  matlab/
    run_scenarios.m                  # Entry point for all scenario runs  [REQUIRED]
    src/
      model_odes.m                   # Coupled ODE system                [REQUIRED]
      params_default.m               # Baseline parameter set            [REQUIRED]
      demand_module.m                # Demand dynamics                   [REQUIRED]
      retrofit_adoption.m            # Retrofit adoption dynamics        [REQUIRED]
      emission_intensity.m           # Emission‑intensity dynamics       [REQUIRED]
      utils/
        configure_paths.m            # Path handling                     [RECOMMENDED]
        save_results.m               # I/O helpers                       [REQUIRED]
        plot_trajectories.m          # Figures for main text             [RECOMMENDED]

  r/
    run_sensitivity.R                # Entry point for sensitivity runs  [REQUIRED]
    src/
      sensitivity_setup.R            # Factor space & sampling           [REQUIRED]
      sobol_indices.R                # Sobol/GSA routines                [REQUIRED]
      time_resolved_sensitivity.R    # Phase‑aware index computation     [REQUIRED]
      plots_sensitivity.R            # Sensitivity figures               [RECOMMENDED]
      utils_io.R                     # I/O helpers                       [RECOMMENDED]

data/
  input/
    fleet_baseline.csv               # Baseline fleet composition        [REQUIRED]
    retrofit_options.csv             # Retrofit levers/options           [REQUIRED]
    demand_baseline.csv              # Traffic/demand baseline           [REQUIRED]
    emission_factors.csv             # Fuel/engine factors               [REQUIRED]
  output/                            # Auto‑created; results & figures

docs/
  figures_map.csv                    # Map: figure/table -> script       [RECOMMENDED]
  reproduction_log.md                # Notes from a full rerun           [AUTO]

env/
  matlab_versions.txt                # `ver -support` dump               [AUTO]
  R_sessionInfo.txt                  # `sessionInfo()` dump              [AUTO]
```

*If your local structure differs, keep the functional equivalence and update the paths in the
launcher scripts (`run_scenarios.m`, `run_sensitivity.R`).*

---

## 2. Software requirements

- **MATLAB**: R2023b or newer, with Signal Processing and Optimization Toolboxes if used by your
  scripts. Exact version/Toolboxes should be recorded (see §9).  
- **R**: 4.2+ (tested ≥ 4.3). Required packages will be installed on first run (see §5.2).

> **Citing software**: please cite MATLAB (The MathWorks, Inc., Natick, MA) and R (R Core Team)
> in the manuscript references. Example citations are provided at the end of this README (§10).

---

## 3. Quick start (default paths)

### 3.1 MATLAB — scenario simulations
```matlab
% From the repository root
cd code/matlab
addpath(genpath('src'));
run('run_scenarios.m');
```
Expected outputs (created if absent):
- `data/output/scenarios/*.mat` – trajectories and summary metrics per scenario
- `data/output/figures/*.png` – figures used in the main text and SI

### 3.2 R — time‑resolved global sensitivity
```r
# From the repository root
setwd("code/r")
source("run_sensitivity.R")
```
Expected outputs:
- `data/output/sensitivity/indices_*.rds` – time‑resolved indices
- `data/output/figures_sensitivity/*.png` – sensitivity plots

---

## 4. Configuration

Edit the configuration blocks at the top of `run_scenarios.m` and `run_sensitivity.R`:

- **I/O**: `DATA_IN = ../../data/input`, `DATA_OUT = ../../data/output`
- **Random seeds**: set fixed seeds for reproducibility (`rng(1234,'twister')` in MATLAB;
  `set.seed(1234)` in R) where applicable.
- **Parallelism**: enable/disable parallel pools (MATLAB `parpool`; R `future`, `parallel`).

---

## 5. Reproducing figures & tables

### 5.1 Main text figures (MATLAB)
- **Fig. 1–3** (system trajectories under baseline and policy levers):
  `run_scenarios.m` → calls `model_odes.m`, `params_default.m`, plotting in `plot_trajectories.m`.
- **Fig. 4–5** (scenario comparison & cumulative CO₂):
  `run_scenarios.m` with scenario list defined in the script; aggregation via `save_results.m`.

### 5.2 Sensitivity figures (R)
- **Fig. 6–8** (phase‑aware sensitivity indices): `run_sensitivity.R` → `time_resolved_sensitivity.R`.
- **Fig. 9–10** (parameter dominance maps): `sobol_indices.R` + `plots_sensitivity.R`.

*On first run, `run_sensitivity.R` will attempt to install missing CRAN packages automatically
(e.g., `sobolro`/`sensitivity`, `data.table`, `ggplot2`). If your environment forbids automatic
install, pre‑install manually and re‑run.*

---

## 6. Data inputs

Place the following CSV files in `data/input/` (UTF‑8, comma‑separated, header in first row):

| File | Required columns (minimal) |
|---|---|
| `fleet_baseline.csv` | `aircraft_type, units, avg_age, utilization_hrs` |
| `retrofit_options.csv` | `option_id, eff_reduction, adoption_cap, lead_time_yrs, cost_proxy` |
| `demand_baseline.csv` | `year, rpk, gdp_proxy` |
| `emission_factors.csv` | `aircraft_type, fuel_burn_kg_per_hr, co2_factor_kg_per_kg` |

*If your manuscript uses different schemas, adjust the loaders in `save_results.m` / `utils_io.R`.*

---

## 7. Outputs

- **Numerical results**: MATLAB `.mat` files and R `.rds` objects under `data/output/`.
- **Graphics**: `.png`/`.pdf` figures saved under `data/output/figures*/` with filenames prefixed by
  figure numbers to ease cross‑checking during review.

---

## 8. Execution time & hardware

Typical runtime (reference: 8 cores @ ≥3.0 GHz, 32 GB RAM):
- MATLAB scenarios: 5–20 min depending on the number of levers/scenarios.
- R sensitivity: 10–40 min depending on the sampling size and time grid.

---

## 9. Environment capture (strongly recommended)

After a successful run, please record exact environments:

- **MATLAB**
  ```matlab
  cd ../../env
  ver -support > matlab_versions.txt
  ```
- **R**
  ```r
  writeLines(capture.output(sessionInfo()), "../../env/R_sessionInfo.txt")
  ```

These files accompany the submission to support reproducibility.

---

## 10. How to cite

Please cite the article and the software used. Examples for software:

- **MATLAB** — *The MathWorks, Inc. (2023). MATLAB (R2023b), Natick, Massachusetts: The MathWorks, Inc.*
- **R** — *R Core Team (2024). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria.*

If you deposit this code in a repository (e.g., Zenodo/Mendeley Data), include the DOI here and add a
`CITATION.cff` file at the repository root.

---

## 11. License

**Code license**: *to be specified by the author* (e.g., MIT/BSD/Apache‑2.0).  
**Documentation** (this README): CC BY 4.0.

---

## 12. Contact

Corresponding author: *Francisco A. Buendía Hernández*  
Email: *franciscoantonio.bue@uah.es*

---

## 13. Missing files checklist (to be completed by the author)

Mark ✅ when provided; otherwise attach the file(s) and keep the relative paths used above.

- [ ] `code/matlab/run_scenarios.m`
- [ ] `code/matlab/src/model_odes.m`
- [ ] `code/matlab/src/params_default.m`
- [ ] `code/matlab/src/demand_module.m`
- [ ] `code/matlab/src/retrofit_adoption.m`
- [ ] `code/matlab/src/emission_intensity.m`
- [ ] `code/matlab/src/utils/configure_paths.m`
- [ ] `code/matlab/src/utils/save_results.m`
- [ ] `code/matlab/src/utils/plot_trajectories.m`
- [ ] `code/r/run_sensitivity.R`
- [ ] `code/r/src/sensitivity_setup.R`
- [ ] `code/r/src/sobol_indices.R`
- [ ] `code/r/src/time_resolved_sensitivity.R`
- [ ] `code/r/src/plots_sensitivity.R`
- [ ] `code/r/src/utils_io.R`
- [ ] `data/input/fleet_baseline.csv`
- [ ] `data/input/retrofit_options.csv`
- [ ] `data/input/demand_baseline.csv`
- [ ] `data/input/emission_factors.csv`

```
Tip: keep file names short (ASCII) and avoid spaces/symbols to prevent upload issues in Editorial
Manager. Use one “.” (the extension) and avoid duplicates.
```

