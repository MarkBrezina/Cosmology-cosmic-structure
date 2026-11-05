import warnings
exec(open('/home/claude/lab/lab_analysis.py').read().split("t0 = time.time()\nBASE, _")[0])
BASE, _ = evaluate(UNIFORM)
warnings.simplefilter("error")
for xh, n in [(31.6, 2.0), (17.78, 2.0), (10.0, 2.0), (5623, 3.0), (100, 4.0)]:
    law = DarkEnergyLaw(x_half=xh, n=n)
    try:
        r, lm = evaluate(law)
    except Exception as e:
        print("ERROR", xh, n, repr(e)); 
        warnings.simplefilter("ignore"); r, lm = evaluate(law); warnings.simplefilter("error")
    print(f"x_half={xh} n={n}: Om={r['Om']:.4f} d_exp={r['chi2_exp']-BASE['chi2_exp']:.2f} d_vel={r['chi2_vel']-BASE['chi2_vel']:.2f} k={r['k_pred']:.3f} kappa0={r['kappa0']:.4f} meanS0={lm.mean_S0:.4f}")
    print("   kappa at z=0..3:", np.round(lm.kappa(np.array([0, 0.5, 1, 2, 3, 5, 10, 100])), 3))
    print("   f at z:", np.round(lm.f(np.array([0, 0.5, 1, 2, 3, 5, 10, 100])), 4))
