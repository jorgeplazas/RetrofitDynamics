# Corrected RetrofitDynamics code

The mathematical right-hand sides from **Mod3 and Mod4_c1 are preserved**. Mod1,
Mod2 and all R entry points now use that same model. MATLAB scenario labels and
post-processing have been corrected. Original saved Live Script figures were
cleared so they cannot be mistaken for outputs from the revised code.

## Model and interpretation

```text
x' = a11*x + a12*q
q' = a21*x + a22*q - beta*q*z
z' = eta*(1-z) - gamma*z + (delta/Nx)*q
```

The R state named `y` is the manuscript's `q`, and R `Beta` is MATLAB `beta`.
R uses `chi = delta/Nx` directly because these two historical constants enter
only through their ratio. Sampling them separately would obscure that dependence.

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

These are dimensionally consistent **illustrative assumptions**, not established
empirical calibration. `delta=0.02` and `Nx=1e12` remain in MATLAB to preserve the
original parameterisation, but delta is not an obsolescence rate in this model.
No claim of 30%, 60% or 90% aircraft efficiency is made. Correcting that label
does not independently establish a realistic mitigation magnitude.

## Files and changes

- `Mod1.mlx`, `Mod2.mlx`: original alternative equations replaced with the retained
  Mod3/Mod4 equations. Both show the three state trajectories at eta=0.6.
- `Mod3.mlx`, `Mod4_c1.mlx`: original equations, numerical parameter values and
  initial states retained. Updated rate labels, explicit line-handle legends,
  exact 20-year slope calculation, units, solver options and domain checks.
  Unsupported stabilization shading and "Growth resumes" annotations removed.
- `Mod1.m`, `Mod2.m`, `Mod3.m`, `Mod4_c1.m`: matching plain-text exports. Each is
  self-contained with local functions; use either .m or .mlx, not both at once.
- `modelo_SA_retrofit_dinamico.R`: canonical R model, configuration, integration,
  Martinez estimator, bootstrap and exports. It defines functions when sourced
  and runs the study when invoked directly with Rscript.
- `paper_sobol_reproducible.R`: paper workflow entry point, sourcing the core.
- `R_Analisis_Sobol.txt`: R-sourceable entry point for the same workflow. It no
  longer contains a separate, potentially divergent model implementation.
- `validation/`: executable R tests, independent SciPy checks and numerical
  reference values. `VALIDATION.md` records what was actually tested.

Keep the three R/text files together. The entry points locate the core relative
to their own file, so they can be run from another working directory. Outputs
are written beneath the current working directory unless configured otherwise.

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
not a fourth dynamic state. Refine the output mesh before relying on precise
cumulative differences. Running the same MATLAB script again overwrites its
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

The shipped scenarios deliberately retain eta=0.3, 0.6, 0.9 even when sensitivity
bounds are edited. The delivered study enforces a 50-year horizon. To generalise it, also adapt
`check_config`, the fixed 20-year slope, and final-year column labels in
`write_scenarios`; do not reuse those labels under a different horizon.

## What changed in the sensitivity calculation

1. **Same model, same initial values and same 50-year horizon as MATLAB.**
   deSolve `lsodar` includes the initial time zero for integration, while indices
   are evaluated at years 1 through 50. Fixed initial states have zero variance,
   so Sobol indices at time zero would be undefined.
2. **Eight identifiable inputs.** The shipped independent uniform ranges are
   newly declared exploratory ranges: +/-10% around the retained baseline for
   all inputs except eta, which spans 0.3 to 0.9 /year. They are not empirical
   confidence intervals and are not the previous manuscript's nine-input ranges.
   Results change when these distributional assumptions change.
3. **Explicit Martinez estimator.** Independent matrices A and B are sampled;
   C_i is A with column i replaced by B's column i. First order is
   `cor(f(B), f(C_i))`; total order is `1-cor(f(A), f(C_i))`. This is the
   correlation-based estimator, not Jansen's squared-difference estimator.
   Its formulas match the [sensitivity source](https://github.com/cran/sensitivity/blob/master/R/sobolmartinez.R).
4. **Direct solver control.** The fragile ODEsensitivity wrapper is removed.
   Solver tolerances and maximum step are passed directly to
   [deSolve](https://cran.r-project.org/web/packages/deSolve/deSolve.pdf).
5. **Raw estimates and uncertainty.** No clipping to [0,1], no forced first-order
   sum of one, and no forced interpretation of sensitivity as fleet share.
   Undefined variance returns NA. The pointwise percentile bootstrap resamples
   the same row indices across A, B and all hybrids. Intervals describe sampling
   uncertainty, not physical-model uncertainty or simultaneous time bands.
6. **Evidence instead of preset phases.** No hard-coded rankings, efficiency
   interpretations or autonomous gamma/(gamma+delta) equilibrium map. Plots
   show all eight inputs, with ranges allowing out-of-[0,1] estimates. Compare
   sample-size/seed variation and interval widths before claiming convergence.

## Saved R outputs

The output root records configuration in both RDS and readable R form, input
bounds/provenance, session information, scenario curves/summary, and a combined
sample-size/seed table. Each seed directory contains the sampled designs, raw
model outputs, first/total-order CSVs, pointwise intervals and time-series plots.
Time is an explicit CSV column. Wide matrices are time-by-parameter; the long
index table also identifies state, parameter, N and seed. Sensitivity outputs
are intentionally not labelled as reproducing the submitted paper.

## Domain and interpretation limits

The unchanged model does not preserve every physically plausible state for
arbitrary parameters. At x=0 with q>0, x' is negative; at z=1, z'=-gamma+chi*q
can be positive. Both languages stop at detected domain boundaries and check
returned states. R saves a failed design row and aborts instead of deleting or
resampling it, which would change the independent-input probability model.
These checks are numerical safeguards, not a proof of global invariance or a
validation of arbitrary parameter ranges.

Do not use old sensitivity rankings to explain the newly unified model. Run a
justified final design, inspect convergence, and update manuscript figures and
claims from those results. The empirical calibration and economic interpretation
issues raised by the reviewers still require evidence beyond code consistency.
