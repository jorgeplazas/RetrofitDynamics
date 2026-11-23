# Independent numerical checks for Reviewer 6

These supplementary calculations were reconstructed from Equations (1)-(3),
Tables 1 and 4, and the initial conditions printed in the supplied revised
manuscript. They are new checks prepared on 1 November 2026. They do not replace
the authors' original R/MATLAB implementation, time-dependent sensitivity
figures, or the full reproducibility archive mentioned in the manuscript.

## Run

Python 3 with NumPy and SciPy is sufficient:

```bash
python3 analysis.py
```

The supplied outputs were generated with NumPy 2.3.5 and SciPy 1.17.0. The
`results` directory contains all evaluated input designs, cumulative outputs,
bootstrap estimates, summary indices and factorial results. No empirical
calibration is performed or claimed.

## Cumulative-output sensitivity

The parameter order is a11, a12, a21, a22, beta, eta, gamma, delta. The main
independent uniform distributions are exactly the printed Table 4 ranges.
The equal-relative-range variant replaces eta in [0.30, 0.90] with eta in
[0.54, 0.66], retaining the other ranges.

Each distribution set is analysed with (1) fixed initial states, (2) uncertain
initial activity only, and (3) all three initial states uncertain. The assumed
initial-state ranges are x0 in [4.5e11, 5.5e11] pkm/year, y0 in [0.9e9, 1.1e9]
kg/year, and z0 in [0.09, 0.11]. These are exploratory +/-10% perturbations,
not empirically estimated distributions. Dependence between states and
parameters is not inferred or tested.

Two independent NumPy default_rng seeds, 20261101 and 20261102, are used with
N=1000 rows in each of A and B. These are new Python designs, not the original
R designs. AB_i replaces column i of A with column i of B. Each run requires
N(k+2) trajectories. The Martinez estimates are correlation(B, AB_i) for S1
and 1-correlation(A, AB_i) for ST, as specified in the manuscript. Bootstrap
intervals are pointwise 95% percentile intervals from 200 paired-row resamples,
using a separate RNG with seed increased by 100. Raw estimates are retained.
Nonzero S1 estimates for nearly irrelevant inputs can be dominated by Monte
Carlo error. They are not interpreted as evidence of effects exceeding ST.

E(50) is computed by augmenting the numerical integration with E'=y, E(0)=0.
The original three governing equations are unchanged. Scaled state variables
x/Nx, y/Nx, E/Nx are used numerically, with Nx=1e12 fixed; multiplying the
scaled cumulative output by 1000 converts it to Mt (1 Mt = 1e9 kg).
Classical fourth-order Runge-Kutta uses a 0.05-year step. Four trajectories
from each design are checked against DOP853 with rtol=1e-11, atol=1e-13 and
maximum step 0.1 year. A full-design step-refinement check at 0.025 year is
included for the main-range, all-initial-states design with seed 20261101.
The observed maximum relative difference in this full-design refinement was
2.77e-10. This is numerical verification, not empirical validation.

Every integration step is checked for the physical domain. Domain extrema
and solver comparisons are recorded in analysis_results.json. Do not
interpret stepwise checks alone as a theorem of domain preservation; see
the sufficient conditions below.

Each design_outputs_*.npz contains:

- A and B: the independent scaled input matrices.
- lower and upper: bounds for the input columns.
- E50_scaled: cumulative output, with rows for A, B, then AB_i in input order.
- S1_bootstrap and ST_bootstrap: all 200 bootstrap index vectors.

## Complete factorial comparison

G denotes a11 changed from 0.05 to 0.03; T denotes a22 changed from -0.04 to
-0.06; D denotes eta changed from 0.6 to 0.8. All other values remain at their
printed baselines. All eight combinations are solved with DOP853 and independently
checked with Radau, integrating E as an additional output state. The results
are in factorial_scenarios.csv and analysis_results.json.

These independently integrated cumulative values differ slightly in their last
reported decimal places from some values printed in Table 3. For example,
the present integral gives 2237.69510 Mt for G alone, while the manuscript prints
2237.68 Mt. Reconcile numerical quadrature and rounding in the original
pipeline before adopting a replacement table. Percentage conclusions agree
at the manuscript's stated precision. Do not mix implementations silently.

On the absolute cumulative-emissions scale, define

I = E111 - E100 - E010 - E001 + 2 E000.

I is approximately +105.86213 Mt. Thus the combined reduction is smaller
than the sum of the three reductions applied separately to baseline. This
comparison does not demonstrate positive additive synergy. Other interaction
scales require their own explicit definitions.

## Sufficient conditions for physical admissibility over a finite horizon

Assume x0>0, a11>=0, a12<=0, a21>=0, beta>=0, eta>=0, gamma>0, delta>0,
and Nx>0. Suppose a positive q (in kg/pkm) satisfies

1. 0 <= y0 <= q*x0 and 0 <= z0 <= 1;
2. a21 + (a22-a11)*q - a12*q*q <= 0;
3. q*x0*exp(a11*T) <= gamma*Nx/delta.

Then x>0, y>=0 and 0<=z<=1 on [0,T].

Proof: On y=q*x,

  (y-q*x)' = x*(a21 + (a22-a11-beta*z)*q - a12*q*q) <= 0

while z>=0, so the cone 0<=y<=q*x is preserved. In that cone,
x' >= (a11+a12*q)*x, so positive x cannot reach zero in finite time,
and x' <= a11*x, so x<=x0*exp(a11*t). At y=0, y'=a21*x>=0.
At z=0, z'=eta+delta*y/Nx>=0. At z=1,
z'=-gamma+delta*y/Nx<=0 by condition 3. These estimates also bound every
state over the finite interval, allowing the local solution to extend to T.

For all main Table 4 combinations, fixed initial states and T=50, q=0.5
works uniformly. The worst value of condition 2 is
0.033 + (-0.036-0.045)*0.5 + 0.022*0.5^2 = -0.002.
The upper bound on y is 3.91066e12 kg/year, below
min(gamma*Nx/delta)=4.09091e12 kg/year.

For the +/-10% initial-state ranges, q=0.47 works uniformly over Table 4:
the condition-2 bound is -0.0002102, and the condition-3 y bound is
4.04362e12 kg/year, again below 4.09091e12 kg/year. These statements
also cover the narrowed eta range.

For the eight fixed factorial combinations, q=0.5 works: the largest
condition-2 value is zero, and the y bound is 3.04562e12 kg/year, below
gamma*Nx/delta=5e12 kg/year. This establishes admissibility for these
domains and this horizon, not for unrestricted coefficients or all time.

## Interpretation of deployment parameters

Ignoring the emissions-responsive term, the deployment equation has
z_bar=eta/(eta+gamma), relaxation time tau=1/(eta+gamma), and time to
complete 90% of the initial-to-equilibrium adjustment t90=log(10)*tau.
For gamma=0.1 and eta=0.3,0.6,0.9, these give tau=2.5,1.42857,1 years,
t90=5.75646,3.28941,2.30259 years, and z_bar=0.75,0.85714,0.90.
These are not times to retrofit 90% of the fleet. The full model also
contains the emissions-dependent driving term and has no fixed z_bar
while emissions vary.

## Files not available for these checks

The separate high-resolution figure files and the authors' original
submission ZIP were not supplied with this request. No claim is made that
they have been repaired, regenerated, audited or replaced. The proposed
figure consistency and availability statements require those actual files
to be checked and the full submission archive to be supplied.
