> Historical record from 11 September 2026. The full R studies have since
> completed; see the parent archive README and results dated 12 September.

# Validation record

Validated against the supplied seven files on 11 September 2026. No remote
GitHub repository was changed. The submitted manuscript's numerical sensitivity
results were not used as validation targets, because they use another model.

## Completed checks

- Four MATLAB .m files contain the same three derivative expressions as the
  original Mod3/Mod4_c1. Baseline coefficients and scenario entry rates are retained.
- Each rebuilt .mlx passes ZIP integrity and XML parsing checks; its extracted
  code equals the corresponding .m export. Stale embedded figure outputs are removed.
- Native R execution: **R 4.2.2 Patched, deSolve 1.34**, Linux. The runtime was
  installed into a temporary local directory for testing, not included in this package.
- Native R baseline integrations at eta=0.3, 0.6, 0.9 agree with independent SciPy
  DOP853 integrations. SciPy Radau independently agrees with DOP853. Checks use
  all three states; reference points include exactly years 0, 20 and 50.
- Native R solver-tolerance refinement agrees at the test threshold; direct
  sign checks confirm that increased deployment lowers the local emissions tendency.
- An intentionally out-of-domain trajectory is rejected.
- The Martinez implementation passes additive and interaction benchmarks with
  known analytical first-order and total-order indices, using 100,000 samples.
- Zero-variance correlations produce NA, not fabricated zero indices.
- Full workflow smoke test: two seeds, N=32, sample-size levels 8/16/32,
  20 bootstrap replicates, all 50 annual evaluation times and all three states.
  This evaluates 640 sampled ODE trajectories plus baseline scenarios. The test
  checks dimensions, finite indices/intervals, explicit time columns, and refusal
  to mix outputs with an existing run.
- Paper and .txt entry points were additionally exercised from another working
  directory, including sourcing the .txt file, with N=8 and four bootstrap replicates.

The exact executable checks are `validation/run_tests.R` and
`validation/check_independent.py`. Run from this package directory:

```sh
Rscript validation/run_tests.R
python3 validation/check_independent.py
```

To compare R output trajectories on the full 500-point mesh, pass the R output
folder as the Python script's argument. The Python check needs NumPy and SciPy;
the R production workflow needs only deSolve beyond base R.

## Numerical reference (independent integration)

| eta (/year) | q(50) (kg/year) | Exact 0-20 year secant slope (kg/year^2) |
|---|---:|---:|
| 0.3 | 2.50852373567e11 | 2.850708214e9 |
| 0.6 | 2.25640002042e11 | 2.537957688e9 |
| 0.9 | 2.16884913127e11 | 2.431345629e9 |

The 20-year slope now uses exactly t=20; the older MATLAB grid endpoint was
19.93987976 years. This intentional correction changes the last digits of the
slope without changing the ODE.

## Remaining validation limits

MATLAB and Live Editor were not available. MATLAB source/container checks and
independent numerical agreement do not establish native .mlx execution or UI
compatibility. Matching .m scripts are included for review and native execution.

The default N=1000, two-seed, 200-bootstrap study was not run here. Smoke-test
indices are not publication results and establish neither convergence nor the
empirical validity of the illustrative uncertainty ranges. No empirical fit,
held-out validation, or reproduction of every submitted figure is claimed.
