import sys, math
src = open('/home/claude/lab/lab2.py').read()
exec(src.split("FAMILIES = {")[0])
_ZD2 = ZH[ZH < 2e4]
def void_fraction_history(x_void, beta, p):
    s = np.sqrt(np.log1p((sigma_tophat(2.0) * D_LCDM) ** 2))
    fv = special.ndtr((math.log(x_void) + s ** 2 / 2) / s)
    f = (fv / fv[0]) ** beta * (1 + _ZD2) ** (6 * (1 - p))
    z = np.concatenate([_ZD2, np.geomspace(_ZD2[-1] * 1.01, 1e8, 30)])
    return History(z, np.concatenate([f, np.full(30, f[-1])])), fv
for cname, use in {"SN+BAO+flows": ("sn", "bao", "vel"), "SN+BAO+CMB+flows": ("sn", "bao", "cmb", "vel")}.items():
    base = fit_family(lambda: UNI, [[0.31]], [], use)["chi2"]
    r = fit_family(lambda lx, beta, p: void_fraction_history(10 ** lx, beta, p)[0], [[0.31, -0.5, 0.5, 0.95], [0.31, -1.0, 1.0, 0.9]],
                   [(-3.0, -0.05), (0.0, 6.0), (0.6, 1.3)], use)
    lx, beta, p = r["params"]; h, fv = void_fraction_history(10 ** lx, beta, p)
    i1, i2 = int(np.argmin(abs(_ZD2 - 1))), int(np.argmin(abs(_ZD2 - 2)))
    print(f"{cname}: void-share law Δχ²={r['chi2'] - base:+.2f}  x_void={10**lx:.3g}, beta={beta:.2f}, p={p:.3f}; void share today {fv[0]:.2f}, z=1 {fv[i1]:.2f}, z=2 {fv[i2]:.2f}; "
          f"w0={float(h.w_eff(0)):+.3f} w(1)={float(h.w_eff(1)):+.3f}; ρDE(0.5)={float(h.f(0.5)):.3f} ρDE(2)={float(h.f(2)):.3f}")
