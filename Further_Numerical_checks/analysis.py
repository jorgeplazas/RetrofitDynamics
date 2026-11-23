"""Independent checks reconstructed from the equations printed in the manuscript.

This supplementary calculation does not replace the authors' original R pipeline.
"""
import json
import time
import csv
from pathlib import Path
from itertools import product

import numpy as np
from scipy.integrate import solve_ivp

BASE = np.array([.05, -.02, .03, -.04, .8, .6, .1, .02])
LO = np.array([.045, -.022, .027, -.044, .72, .30, .09, .018])
HI = np.array([.055, -.018, .033, -.036, .88, .90, .11, .022])
SEEDS = [20261101, 20261102]
N = 1000
B = 200
OUTPUT = Path(__file__).resolve().parent / 'results'
OUTPUT.mkdir(exist_ok=True)


def rhs(w, p):
    a, b, c, d, beta, eta, gamma, delta = p[:, :8].T
    X, Y, z, C = w
    return np.array([a*X+b*Y, c*X+(d-beta*z)*Y,
                     eta*(1-z)-gamma*z+delta*Y, Y])


def integrate_many(p, dt=.05):
    # X=x/Nx, Y=y/Nx, C=E/Nx, with fixed Nx=10^12.
    w = np.zeros((4, len(p)))
    w[0] = p[:, 8] if p.shape[1] > 8 else .5
    w[1] = p[:, 9] if p.shape[1] > 9 else .001
    w[2] = p[:, 10] if p.shape[1] > 10 else .1
    domain = [float(w[0].min()), float(w[1].min()), float(w[2].min()), float(w[2].max())]
    for _ in range(round(50/dt)):
        k1 = rhs(w, p)
        k2 = rhs(w+dt*k1/2, p)
        k3 = rhs(w+dt*k2/2, p)
        k4 = rhs(w+dt*k3, p)
        w += (dt/6)*(k1+2*k2+2*k3+k4)
        domain[0] = min(domain[0], float(w[0].min()))
        domain[1] = min(domain[1], float(w[1].min()))
        domain[2] = min(domain[2], float(w[2].min()))
        domain[3] = max(domain[3], float(w[2].max()))
    if domain[0] < 0 or domain[1] < 0 or domain[2] < 0 or domain[3] > 1:
        raise ValueError(domain)
    return w[3], domain


def corr(a, b):
    ac = a-a.mean(axis=-1, keepdims=True)
    bc = b-b.mean(axis=-1, keepdims=True)
    return (ac*bc).sum(axis=-1)/np.sqrt((ac*ac).sum(axis=-1)*(bc*bc).sum(axis=-1))


def gsa(narrow, initial, seed):
    lo, hi = LO.copy(), HI.copy()
    if narrow:
        lo[5], hi[5] = .54, .66
    if initial:
        lo, hi = np.r_[lo, .45], np.r_[hi, .55]
    if initial == 3:
        lo, hi = np.r_[lo, .0009, .09], np.r_[hi, .0011, .11]
    k = len(lo)
    rng = np.random.default_rng(seed)
    A = rng.uniform(lo, hi, (N, k))
    D = rng.uniform(lo, hi, (N, k))
    matrices = [A, D]
    for i in range(k):
        h = A.copy()
        h[:, i] = D[:, i]
        matrices.append(h)
    P = np.vstack(matrices)
    values, domain = integrate_many(P)
    values = values.reshape(k+2, N)
    si = corr(values[1], values[2:])
    st = 1-corr(values[0], values[2:])
    rb = np.random.default_rng(seed+100)
    bs, bt = [], []
    for _ in range(B):
        j = rb.integers(N, size=N)
        v = values[:, j]
        bs.append(corr(v[1], v[2:]))
        bt.append(1-corr(v[0], v[2:]))
    cis = np.percentile(bs, [2.5, 97.5], axis=0).T
    cit = np.percentile(bt, [2.5, 97.5], axis=0).T
    # Check a random subset against an adaptive high-accuracy solver.
    verify = []
    for i in [0, N-1, 2*N+3, len(P)-1]:
        p = P[i:i+1]
        w0 = [p[0, 8] if initial else .5,
              p[0, 9] if initial == 3 else .001,
              p[0, 10] if initial == 3 else .1, 0]
        sol = solve_ivp(lambda t, w: rhs(w[:, None], p)[:, 0],
                        [0, 50], w0, method='DOP853',
                        rtol=1e-11, atol=1e-13, max_step=.1)
        verify.append(abs(values.ravel()[i]/sol.y[3, -1]-1))
    mode = 'all_initial' if initial == 3 else ('uncertain_x0' if initial else 'fixed_x0')
    case = f'{"equal_ranges" if narrow else "main_ranges"}_{mode}'
    np.savez_compressed(OUTPUT / f'design_outputs_{case}_{seed}.npz',
                        A=A, B=D, E50_scaled=values, S1_bootstrap=bs,
                        ST_bootstrap=bt, lower=lo, upper=hi)
    refinement_error = None
    if initial == 3 and not narrow and seed == SEEDS[0]:
        refined, _ = integrate_many(P, dt=.025)
        refinement_error = float(np.max(np.abs(refined/values.ravel()-1)))
    return {
        'seed': seed, 'N': N, 'bootstrap': B, 'domain': domain,
        'max_relative_quadrature_check': max(verify),
        'max_relative_step_refinement_error': refinement_error,
        'indices': {name: {'S1':float(s), 'ST':float(t), 'S1_CI':cs.tolist(), 'ST_CI':ct.tolist()}
                    for name, s, t, cs, ct in zip(
                        ['a11','a12','a21','a22','beta','eta','gamma','delta','x0','y0','z0'], si, st, cis, cit)},
        'E50_Mt_quantiles': (np.percentile(values[:2].ravel(), [2.5,50,97.5])*1000).tolist(),
    }


def factorial():
    rows = []
    verification = []
    for g, t, d in product([0, 1], repeat=3):
        p = BASE.copy()
        p[0], p[3], p[5] = [.05,.03][g], [-.04,-.06][t], [.6,.8][d]
        def f(time, w):
            return rhs(w[:,None], p[None,:])[:,0]
        primary = solve_ivp(f,[0,50],[.5,.001,.1,0],method='DOP853',
                            rtol=1e-11,atol=1e-13,max_step=.1)
        check = solve_ivp(f,[0,50],[.5,.001,.1,0],method='Radau',
                          rtol=1e-10,atol=1e-12,max_step=.1)
        relative = float(np.max(np.abs(primary.y[:,-1]/check.y[:,-1]-1)))
        verification.append(relative)
        # Nx=10^12, 1 Mt=10^9 kg, so the conversion factor is 1000.
        rows.append({'G':g,'T':t,'D':d,'a11':p[0],'a22':p[3],'eta':p[5],
                     'y50_Mt_per_year':primary.y[1,-1]*1000,
                     'E50_Mt':primary.y[3,-1]*1000})
    base = rows[0]['E50_Mt']
    for row in rows:
        row['reduction_percent'] = 100*(1-row['E50_Mt']/base)
    with open(OUTPUT / 'factorial_scenarios.csv','w') as f:
        writer=csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    e = {(r['G'],r['T'],r['D']):r['E50_Mt'] for r in rows}
    interaction = e[1,1,1]-e[1,0,0]-e[0,1,0]-e[0,0,1]+2*e[0,0,0]
    return {'rows':rows, 'max_solver_relative_difference':max(verification),
            'combined_additivity_residual_Mt':interaction,
            'additional_reduction_from_growth_only_percent':100*(1-e[1,1,1]/e[1,0,0])}


if __name__ == '__main__':
    out = {'description':'Independent reconstruction; Python RNG designs differ from the authors original R designs.',
           'parameter_order':['a11','a12','a21','a22','beta','eta','gamma','delta','x0','y0','z0'],
           'normalization':'x/Nx, y/Nx and E/Nx; Nx=1e12; scaled E multiplied by 1000 gives Mt.',
           'factorial':factorial(), 'runs':{}}
    started = time.time()
    for narrow, initial in product([False, True], [0,1,3]):
        mode = 'all_initial' if initial == 3 else ('uncertain_x0' if initial else 'fixed_x0')
        key = f'{"equal_ranges" if narrow else "main_ranges"}_{mode}'
        out['runs'][key] = []
        for seed in SEEDS:
            result = gsa(narrow, initial, seed)
            out['runs'][key].append(result)
            print(key, seed, round(time.time()-started,2),
                  {k:round(v['ST'],4) for k,v in result['indices'].items()}, flush=True)
            with open(OUTPUT / 'analysis_results.json','w') as f:
                json.dump(out, f, indent=2)
    with open(OUTPUT / 'sensitivity_indices.csv','w') as f:
        names=['case','seed','parameter','S1','S1_low','S1_high','ST','ST_low','ST_high']
        writer=csv.DictWriter(f,fieldnames=names); writer.writeheader()
        for case, runs in out['runs'].items():
            for run in runs:
                for parameter, entry in run['indices'].items():
                    writer.writerow({'case':case,'seed':run['seed'],'parameter':parameter,
                                     'S1':entry['S1'],'S1_low':entry['S1_CI'][0],'S1_high':entry['S1_CI'][1],
                                     'ST':entry['ST'],'ST_low':entry['ST_CI'][0],'ST_high':entry['ST_CI'][1]})
