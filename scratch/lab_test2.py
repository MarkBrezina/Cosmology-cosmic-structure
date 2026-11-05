import time
exec(open('/home/claude/lab/lab_analysis.py').read().split("t0 = time.time()\nBASE, _")[0])
BASE, _ = evaluate(UNIFORM)

def total_delta(r):
    return (r["chi2_exp"] - BASE["chi2_exp"]) + (r["chi2_vel"] - BASE["chi2_vel"])

# grid scan timing
t0 = time.time()
rows = []
for xh in np.geomspace(1, 1e4, 17):
    for n in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0):
        r, _ = evaluate(DarkEnergyLaw(x_half=xh, n=n))
        rows.append((xh, n, total_delta(r), r["chi2_exp"] - BASE["chi2_exp"], r["chi2_vel"] - BASE["chi2_vel"]))
print(f"scan {len(rows)} laws in {time.time()-t0:.1f}s")
df = pd.DataFrame(rows, columns=["x_half", "n", "d_total", "d_exp", "d_vel"])
print(df.pivot(index="n", columns="x_half", values="d_total").round(1).to_string())

# p scan with uniform shape
for p in (0.85, 0.9, 0.95, 0.97, 1.0, 1.03, 1.05):
    r, lm = evaluate(replace(UNIFORM, p=p, name=f"p={p}"))
    print(f"p={p:.2f} d_total={total_delta(r):7.2f} d_exp={r['chi2_exp']-BASE['chi2_exp']:7.2f} w0={r['w0']:.3f} Om={r['Om']:.3f}")

# optimiser over (log x_half, n, p), cut fixed at 200
def objective(th):
    lxh, n, p = th
    if not (0 <= lxh <= 5 and 0.3 <= n <= 6 and 0.7 <= p <= 1.2): return 1e6
    r, _ = evaluate(DarkEnergyLaw(x_half=10**lxh, n=n, p=p))
    return total_delta(r)
t0 = time.time()
best = None
for start in ([1.0, 1.0, 0.97], [1.5, 2.0, 0.95], [2.5, 1.0, 1.0], [0.7, 0.8, 0.9]):
    res = optimize.minimize(objective, start, method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-3, maxiter=400))
    print("start", start, "->", np.round(res.x, 3), round(res.fun, 2), res.nfev)
    if best is None or res.fun < best.fun: best = res
print(f"optimiser {time.time()-t0:.1f}s; best x_half={10**best.x[0]:.2f} n={best.x[1]:.2f} p={best.x[2]:.3f} delta={best.fun:.2f}")
r, lm = evaluate(DarkEnergyLaw(x_half=10**best.x[0], n=best.x[1], p=best.x[2]))
print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
zz = np.array([0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5])
print("DE(z)/DE0:", np.round(lm.f(zz), 3), " w(z):", np.round(lm.w_eff(zz), 3))
print("volume fractions today: x<0.3:", round(WEB.volume_fraction(0, 0.3), 3), " 0.3-3:", round(WEB.volume_fraction(0.3, 3), 3), " 3-30:", round(WEB.volume_fraction(3, 30), 3), " >30:", round(WEB.volume_fraction(30, 1e12), 4))
