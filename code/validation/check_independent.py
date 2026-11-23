"""Independent SciPy check of the retained MATLAB equations and R baseline CSVs.
Run from the package root: python3 validation/check_independent.py <R-output-directory>
Requires numpy and scipy. Does not execute MATLAB or the MLX user interface.
"""
from pathlib import Path
import sys, csv, re
import numpy as np
from scipy.integrate import solve_ivp
ROOT = Path(__file__).resolve().parents[1]

def rhs(t, s, eta):
    x, q, z = s
    return [.05*x-.02*q, .03*x-.04*q-.8*q*z, eta*(1-z)-.1*z+2e-14*q]

def reference():
    rows=[]
    for eta in [.3,.6,.9]:
        a=solve_ivp(lambda t,s:rhs(t,s,eta),(0,50),[5e11,1e9,.1],
                    method='DOP853',rtol=1e-11,atol=[1e-3,1e-4,1e-13],t_eval=[0,20,50])
        b=solve_ivp(lambda t,s:rhs(t,s,eta),(0,50),[5e11,1e9,.1],
                    method='Radau',rtol=1e-11,atol=[1e-3,1e-4,1e-13],t_eval=[0,20,50])
        assert a.success and b.success
        np.testing.assert_allclose(a.y,b.y,rtol=2e-9,atol=1e-11)
        rows.extend(dict(eta=eta,time=t,x=x,q=q,z=z) for t,(x,q,z) in zip(a.t,a.y.T))
    return rows

def main():
    rows=reference()
    # Every delivered MATLAB file contains identical retained derivative expressions.
    expected=['dxdt=a11*x+a12*co2;',
              'dco2dt=a21*x+a22*co2-beta*co2*z;',
              'dzdt=eta*(1-z)-gamma*z+delta*co2/Nx;']
    for name in ['Mod1','Mod2','Mod3','Mod4_c1']:
        text=re.sub(r'\s+','',(ROOT/(name+'.m')).read_text())
        assert all(line in text for line in expected),name
    if len(sys.argv)>1:
        folder=Path(sys.argv[1])
        for eta in [.3,.6,.9]:
            a=np.genfromtxt(folder/f'scenario_eta_{eta:.1f}.csv',delimiter=',',names=True)
            sol=solve_ivp(lambda t,s:rhs(t,s,eta),(0,50),[5e11,1e9,.1],
                          method='DOP853',rtol=1e-11,atol=[1e-3,1e-4,1e-13],
                          t_eval=a['time_year'])
            for field,v in zip(['x_pkm_per_year','q_kg_per_year','z'],sol.y):
                np.testing.assert_allclose(a[field],v,rtol=2e-7,atol=1e-10)
        print('PASS: native R scenario outputs agree with independent SciPy integrations.')
    print('PASS: two SciPy solvers agree and all four MATLAB RHS expressions match.')

if __name__=='__main__':main()
