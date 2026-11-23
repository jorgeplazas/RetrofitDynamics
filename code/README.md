# RetrofitDynamics code



## Model and interpretation

```text
x' = a11*x + a12*q
q' = a21*x + a22*q - beta*q*z
z' = eta*(1-z) - gamma*z + (delta/Nx)*q
```


| Symbol | Meaning | Units | Baseline |
|---|---|---|---|
| x | Annual passenger activity | passenger-km/year | 5e11 initially |
| q / R y | Annual operational CO2 emissions rate | kg/year | 1e9 initially |
| z | Deployment index; share interpretation within [0,1] | dimensionless | 0.1 initially |
| a11 | Activity growth coefficient | 1/year | 0.05 |
| a12 | Emissions-to-activity feedback | passenger-km/(kg year) | -0.02 |
| a21 | Activity-to-emissions forcing | kg/(passenger-km year) | 0.03 |
| a22 | Background emissions decay | 1/year | -0.04 |
| beta / R Beta | Deployment-dependent emissions decay | 1/year | 0.8 |
| eta | Deployment entry rate, not aircraft efficiency | 1/year | 0.6; scenarios 0.3, 0.6, 0.9 |
| gamma | Deployment exit rate | 1/year | 0.1 |
| chi = delta/Nx | Emissions-responsive deployment coefficient | 1/kg | 2e-14 |

These are dimensionally consistent illustrative assumptions, not established
empirical calibration. `delta=0.02` and `Nx=1e12` remain in MATLAB to preserve the
original parameterisation.



## MATLAB usage

Use MATLAB R2019b or newer. Open one .mlx and select Run, or run the corresponding
.m file. The files were rebuilt from the original MLX XML containers and supplied
with identical plain-text exports. Native Live Editor opening/execution has not
been tested in this environment. If a particular release cannot open a rebuilt
MLX, open the corresponding .m and save it as a Live Script in MATLAB.

Each script uses `ode45` with relative tolerance 1e-9, absolute tolerances
`[1e-2, 1e-3, 1e-12]`, and maximum internal step 0.1 year. All use initial state
`[5e11;1e9;0.1]` and 500 requested output times from 0 to 50 years. The output
mesh is not a fixed integration step. `deval` evaluates the exact 20-year endpoint
for the secant slope, replacing the old nearest-grid-time approximation.

`outputs_matlab/` contains trajectory CSVs and PNGs. Mod3/Mod4 also export
summary CSVs with 20-year slopes, final annual rates and cumulative mass. CO2
plots use Mt/year, with 1 Mt = 1e9 kg. Cumulative mass is `cumtrapz` of the annual
rate over the output mesh; it is approximate numerical post-processing and is
not a fourth dynamic state. Running the same MATLAB script again overwrites its
own prefixed output files; the four script prefixes are distinct.

## R usage

Install the sole non-base dependency once in R:

```r
install.packages("deSolve")
```

Run one entry point from a terminal:

```sh
Rscript paper_sobol_reproducible.R
# Equivalent alternatives:
Rscript modelo_SA_retrofit_dinamico.R
Rscript R_Analisis_Sobol.txt
```

Or in RStudio:

```r
source("paper_sobol_reproducible.R")
# Alternatively: source("R_Analisis_Sobol.txt")
```

The default study integrates 20,000 sample trajectories: N=1000 rows times
(8 parameters + 2 base designs), for each of two independent seeds. It uses
200 paired-row bootstrap resamples per seed and compares nested N=250, 500 and
1000 estimates. Runtime depends on the machine. Baseline scenario integrations
are additional. The default output directory is `outputs_sobol`; it must be
new or empty, preventing accidental mixing with historical results.

For a quick pipeline check, configure a smaller run (this is not publication
analysis):

```r
Sys.setenv(RETROFIT_N="32", RETROFIT_NBOOT="20",
           RETROFIT_SEEDS="20250916,20250917", RETROFIT_OUT="outputs_smoke")
source("paper_sobol_reproducible.R")
```

Use a fresh R session or unset these variables before the larger run. For
custom distributions, tolerances or sample sizes, use the core API:

```r
source("modelo_SA_retrofit_dinamico.R")
cfg <- retrofit_config()
# Edit cfg$lower and cfg$upper only with a justified uncertainty specification.
# Names/order must remain identical to cfg$base. All inputs are independent uniforms.
# Example solver-tolerance check with an otherwise identical design:
cfg$rtol <- 1e-10
cfg$atol <- cfg$atol / 10
run_study(cfg, "outputs_tighter_tolerance")
```


## Saved R outputs

The output root records configuration in both RDS and readable R form, input
bounds/provenance, session information, scenario curves/summary, and a combined
sample-size/seed table. Each seed directory contains the sampled designs, raw
model outputs, first/total-order CSVs, pointwise intervals and time-series plots.
Time is an explicit CSV column. Wide matrices are time-by-parameter; the long
index table also identifies state, parameter, N and seed. 


