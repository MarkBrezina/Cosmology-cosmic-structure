import sys, math, time
src = open('/home/claude/lab/lab2.py').read()
exec(src.split("FAMILIES = {")[0])            # data, likelihoods, CMB, flows, History, cpl, void_surface, fit_family

def average_density(lxh, n, p, x_cut=200.0):
    """Dark energy set by the universe's MEAN matter density only: rho_DE ∝ cap(z) * S(xbar), xbar = (1+z)^3."""
    law = DarkEnergyLaw(x_half=10 ** lxh, n=n, x_cut=x_cut, p=p)
    xbar = (1 + ZH) ** 3
    S = law.shape(xbar)
    S0 = float(law.shape(1.0))
    if not S0 > 0:
        raise ValueError
    return History(ZH, (1 + ZH) ** (6 * (1 - p)) * S / S0)

for cname, use in {"SN+BAO+flows": ("sn", "bao", "vel"), "SN+BAO+CMB+flows": ("sn", "bao", "cmb", "vel")}.items():
    base = fit_family(lambda: UNI, [[0.31]], [], use)["chi2"]
    r = fit_family(average_density, [[0.31, 1.0, 1.0, 0.93], [0.31, 0.5, 2.0, 0.95], [0.31, 1.5, 1.0, 0.9], [0.31, 2.0, 3.0, 0.97]],
                   [(-1.0, 4.0), (0.3, 6.0), (0.6, 1.3)], use)
    h = r["hist"]; lxh, n, p = r["params"]
    print(f"{cname}: average-density switch-on Δχ²={r['chi2'] - base:+.2f}  x_half={10**lxh:.3g}× today's mean, n={n:.2f}, p={p:.3f}, Ωm={r['Om']:.3f}, flow slope={r['parts']['k_pred']:.2f}")
    print(f"     w(0)={float(h.w_eff(0)):+.3f} w(0.5)={float(h.w_eff(0.5)):+.3f} w(1)={float(h.w_eff(1)):+.3f} w(2)={float(h.w_eff(2)):+.3f}  "
          f"ρDE(0.5)={float(h.f(0.5)):.3f} ρDE(1)={float(h.f(1.0)):.3f} ρDE(2)={float(h.f(2.0)):.3f}; S at today's mean = {DarkEnergyLaw(x_half=10**lxh, n=n).shape(1.0):.3f} of the cap")
    r2 = fit_family(lambda lx, lR: void_surface(10 ** lx, 10 ** lR), [[0.31, -1.3, -0.5], [0.31, -1.0, -0.9], [0.31, -2.0, -0.7]],
                    [(-3.0, -0.05), (-1.0, 1.5)], use)
    print(f"{cname}: void-surface law with cells down to 0.1 Mpc: Δχ²={r2['chi2'] - base:+.2f}, x_void={10**r2['params'][0]:.3g}, R={10**r2['params'][1]:.3g} Mpc")

# uncertainty on the effective growth dimension D with the CMB
use = ("sn", "bao", "cmb")
best = fit_family(lambda D: cpl(1 - 2 * D / 3, 0.0), [[0.31, 2.98]], [(2.0, 3.6)], use)
prof = []
for D in np.arange(2.86, 3.101, 0.02):
    r = fit_family(lambda: cpl(1 - 2 * D / 3, 0.0), [[best["Om"]]], [], use)
    prof.append((D, r["chi2"] - best["chi2"]))
inside = [D for D, c in prof if c <= 1.0]
print("D profile:", " ".join(f"{D:.2f}:{c:.1f}" for D, c in prof))
print(f"effective growth dimension with the CMB: D = {best['params'][0]:.3f}, 1σ range {min(inside):.2f}–{max(inside):.2f}")
