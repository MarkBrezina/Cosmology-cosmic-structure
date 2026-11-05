import math, sys
src = open('/home/claude/lab/head_to_head.py').read()
exec(src.split("models = [")[0])            # data, likelihoods, CPL, velocity slope, fit helpers
exec("def timescape_dl" + src.split("def timescape_dl")[1].split("# checks:")[0])   # timescape distances + SN chi2
c_lcdm = {k: fit(make_lcdm, [[0.31]], [(0.05, 0.8)], use).fun for k, use in COMBOS.items()}

# 1) timescape with uncorrected CMB-frame redshifts (no 2M++ peculiar-velocity correction), both models alike
def chi2_zframe(dl_of_z, zcol):
    mu100 = 5 * np.log10((1 + SN.zhel) / (1 + zcol) * (C_KMS / 100) * dl_of_z(zcol)) + 25
    y = SN.mb - np.where(SN.cal, SN.ceph, mu100); Wy = SN.W @ y; b = SN.A.T @ Wy; th = SN.Finv @ b
    return float(y @ Wy - b @ th)
zcmb = pp.zCMB.to_numpy()[SN.idx]
best_l = min(chi2_zframe(lambda z, om=om: (1 + z) * np.interp(z, Z_GRID, background(UNI, om)[1]), zcmb) for om in np.linspace(0.2, 0.45, 101))
best_t = min(chi2_zframe(lambda z, fv=fv: timescape_dl(z, fv), zcmb) for fv in np.linspace(0.4, 0.95, 111))
print(f"uncorrected CMB-frame redshifts: Δχ² timescape − ΛCDM = {best_t - best_l:+.2f}")

# 2) can the physical-density switch-on + slow dilution reproduce DESI's evolving dark energy?
W0, WA = -0.858, -0.58          # DESI DR2 + CMB + Pantheon+ (arXiv:2508.10514, Table IV)
desi = CPL(W0, WA)
zz = np.linspace(0, 2.5, 120)
target = np.log(desi.f(zz))
def mismatch(th):
    lxh, n, p = th
    if not (0 <= lxh <= 4 and 0.3 <= n <= 6 and 0.6 <= p <= 1.2): return 1e9
    lm = LawModel(DarkEnergyLaw(x_half=10 ** lxh, n=n, x_cut=200.0, p=p), WEB)
    return float(np.mean((np.log(np.maximum(lm.f(zz), 1e-300)) - target) ** 2))
best = None
for s in ([0.7, 1.0, 0.85], [1.0, 2.0, 0.9], [0.5, 1.5, 0.8], [1.3, 3.0, 0.92]):
    r = optimize.minimize(mismatch, s, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-8, maxiter=3000))
    if best is None or r.fun < best.fun: best = r
lxh, n, p = best.x
mimic = DarkEnergyLaw(name="mimic", x_half=10 ** lxh, n=n, x_cut=200.0, p=p)
lm_m = LawModel(mimic, WEB)
print(f"\nbest mimic of DESI's history: x_half = {10**lxh:.2f}x mean, n = {n:.2f}, p = {p:.3f}; rms mismatch in ln ρ_DE over 0<z<2.5 = {math.sqrt(best.fun):.4f}")
for z in (0, 0.5, 1.0, 2.0):
    print(f"  z={z}: ρ_DE/ρ_DE0 law {float(lm_m.f(z)):.3f} vs DESI {float(desi.f(z)):.3f};  w law {float(lm_m.w_eff(z)):+.3f} vs DESI {W0 + WA * z / (1 + z):+.3f}")
print(f"  void excess over average: {100 * (mimic.shape(0.2) / mimic.shape(1.0) - 1):.1f}%,  κ today {lm_m.kappa_web[0]:.3f}")
for name, lm in (("DESI published w0wa (fixed)", desi), ("your law tuned to it", lm_m)):
    for cname, use in COMBOS.items():
        r = optimize.minimize_scalar(lambda Om: sum(v for k, v in chi2_parts(lm, Om).items() if k in use), bounds=(0.1, 0.6), method="bounded")
        if cname == "SN+BAO+flows":
            f, D = growth(lm, float(r.x))
            print(f"  {name:28s} {cname}: Δχ² vs ΛCDM {r.fun - c_lcdm[cname]:+.2f}  (Ωm {r.x:.3f}, predicted flow slope {f * D / REF_FD:.2f} vs measured {VEL['k']:.2f} ± {VEL['k_err']:.2f})")
        else:
            print(f"  {name:28s} {cname}: Δχ² vs ΛCDM {r.fun - c_lcdm[cname]:+.2f}")
