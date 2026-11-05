"""Head-to-head: Mark's laws against community models on the same data and likelihoods."""
import math
import sys
import time

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from lab_core import *  # noqa: F401,F403

t_start = time.time()
pp = pd.read_csv(fetch("pantheon"), sep=r"\s+")
_c = pd.read_csv(fetch("pantheon_cov"), header=None).to_numpy().ravel()
PP_COV = _c[1:].reshape(int(_c[0]), int(_c[0]))
_rows = [ln.split() for ln in open(fetch("desi_mean")) if ln.strip() and not ln.startswith("#")]
BAO_TABLE = pd.DataFrame(dict(z=[float(r[0]) for r in _rows], value=[float(r[1]) for r in _rows], kind=[r[2] for r in _rows]))
BAO_COV = np.loadtxt(fetch("desi_cov"))
SN, BAO, WEB = Supernovae(pp, PP_COV), BAOData(BAO_TABLE, BAO_COV), CosmicWeb(R_mpc=2.0)
UNI = LawModel(UNIFORM, WEB)


class CPL:
    """Smooth dark energy with w(z) = w0 + wa z/(1+z) (DESI's evolving dark energy)."""
    def __init__(self, w0, wa):
        self.w0, self.wa = w0, wa
        self.f_grid = self.f(Z_GRID)

    def f(self, z):
        z = np.asarray(z, float)
        return (1 + z) ** (3 * (1 + self.w0 + self.wa)) * np.exp(-3 * self.wa * z / (1 + z))

    def g(self, z):
        return np.zeros_like(np.asarray(z, float))


# ---------------------------------------------------------------- measured velocity slope (as in the notebook)
def measure_velocity_slope(M, H0, Om, n_boot=1000, sig_int=0.103):
    d = pp.copy()
    d["mu"] = d.m_b_corr - M
    d["w"] = 1 / d.m_b_corr_err_RAW ** 2
    d["wmu"] = d.w * d.mu
    g = d.groupby("CID", sort=False)
    u = g.first()
    sums = g[["w", "wmu"]].sum()
    u["mu"] = sums.wmu / sums.w
    u["err"] = 1 / np.sqrt(sums.w)
    u = u[(u.zCMB > 0.01) & (u.zCMB < 0.06)]
    E, Dc = background(UNI, Om)
    chi = 10 ** ((u.mu - 25) / 5) / (1 + u.zHEL)
    zcos = np.interp(chi * H0 / C_KMS, Dc, Z_GRID)
    v = C_KMS * ((1 + u.zCMB) / (1 + zcos) - 1)
    sv = (1 + u.zCMB) / (1 + zcos) ** 2 * chi * LN10_5 * np.sqrt(u.err ** 2 + sig_int ** 2) * H0 * np.interp(zcos, Z_GRID, E)
    x, y, s = u.VPEC.to_numpy(), v.to_numpy(), sv.to_numpy()

    def nll(th, x, y, s):
        var = s ** 2 + math.exp(2 * th[2])
        return 0.5 * np.sum((y - th[0] - th[1] * x) ** 2 / var + np.log(var))
    best = optimize.minimize(nll, [0.0, 1.0, math.log(280.0)], args=(x, y, s), method="Nelder-Mead").x
    boots = [optimize.minimize(nll, best, args=(x[i], y[i], s[i]), method="Nelder-Mead").x[1]
             for i in (rng.integers(0, len(x), len(x)) for _ in range(n_boot))]
    return dict(k=float(best[1]), k_err=float(np.std(boots)))


def sn_only_uniform():
    res = optimize.minimize_scalar(lambda Om: SN.chi2(UNI, Om), bounds=(0.1, 0.6), method="bounded")
    c, p = SN.chi2(UNI, float(res.x), details=True)
    return float(res.x), p


Om_sn, p_sn = sn_only_uniform()
VEL = measure_velocity_slope(p_sn["M"], p_sn["H0"], Om_sn)
_f, _D = growth(UNI, 0.315)
REF_FD = _f * _D
print(f"velocity slope {VEL['k']:.3f} ± {VEL['k_err']:.3f}")


def chi2_parts(lm, Om):
    f, D = growth(lm, Om)
    return dict(sn=SN.chi2(lm, Om), bao=BAO.chi2(lm, Om), vel=((f * D / REF_FD - VEL["k"]) / VEL["k_err"]) ** 2)


def fit(make, starts, bounds, use):
    def obj(th):
        if any(not (lo <= v <= hi) for v, (lo, hi) in zip(th, bounds)):
            return 1e9
        lm, Om = make(th)
        if lm is None:
            return 1e9
        parts = chi2_parts(lm, Om) if "vel" in use else dict(sn=SN.chi2(lm, Om) if "sn" in use else 0.0,
                                                              bao=BAO.chi2(lm, Om) if "bao" in use else 0.0)
        return sum(parts[k] for k in use)
    best = None
    for s in starts:
        r = optimize.minimize(obj, s, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-6, maxiter=4000, maxfev=4000))
        if best is None or r.fun < best.fun:
            best = r
    return best


COMBOS = {"SN": ("sn",), "BAO": ("bao",), "SN+BAO": ("sn", "bao"), "SN+BAO+flows": ("sn", "bao", "vel")}


def make_lcdm(th):
    return UNI, th[0]


def make_wcdm(th):
    return CPL(th[1], 0.0), th[0]


def make_cpl(th):
    if th[1] + th[2] >= 0:
        return None, None
    return CPL(th[1], th[2]), th[0]


MARK = DarkEnergyLaw(name="Mark's capped law", x_half=30.0, n=1.0, x_cut=200.0, cut_width_dex=0.15, p=1.0)
MARK_LM = LawModel(MARK, WEB)


def make_mark(th):
    return MARK_LM, th[0]


models = [
    ("ΛCDM (uniform dark energy)", make_lcdm, [[0.31]], [(0.05, 0.8)], 0),
    ("wCDM = your compounding law (w = 1 − 2p)", make_wcdm, [[0.31, -1.0], [0.3, -0.9]], [(0.05, 0.8), (-2.0, 0.0)], 1),
    ("w0waCDM (DESI's evolving dark energy)", make_cpl, [[0.31, -0.9, -0.3], [0.3, -0.75, -0.8], [0.32, -1.0, 0.1]],
     [(0.05, 0.8), (-2.0, 0.0), (-4.0, 2.0)], 2),
    ("Your capped density law (30×, n=1, stop 200×)", make_mark, [[0.31]], [(0.05, 0.8)], 0),
]
results = {}
for name, make, starts, bounds, k in models:
    results[name] = {}
    for cname, use in COMBOS.items():
        r = fit(make, starts, bounds, use)
        results[name][cname] = (r.fun, r.x, k)
base = results["ΛCDM (uniform dark energy)"]
print("\nΔχ² relative to ΛCDM (extra parameters in brackets); best-fit parameters on SN+BAO+flows")
for name, d in results.items():
    line = f"{name:50s} [{d['SN'][2]}]"
    for cname in COMBOS:
        line += f"  {cname}: {d[cname][0] - base[cname][0]:+7.2f}"
    print(line, "  params:", np.round(d["SN+BAO+flows"][1], 3))


# ---------------------------------------------------------------- timescape (Wiltshire tracker solution), supernovae only
def timescape_dl(zq, fv0):
    """Dressed luminosity distance in units of c/H̄0 (bare Hubble constant)."""
    t0 = (2 + fv0) / 3
    b = 2 * (1 - fv0) * (2 + fv0) / (9 * fv0)
    t = np.geomspace(t0 * 1e-5, t0, 40000)
    fv = t / (t + 1.5 * b)
    onepz = (2 + fv) * fv ** (1 / 3) / (3 * fv0 ** (1 / 3) * t)
    cb = b ** (1 / 3)

    def F(tt):
        u = tt ** (1 / 3)
        return 2 * u + (cb / 6) * np.log((u + cb) ** 2 / (u * u - cb * u + cb * cb)) + (cb / math.sqrt(3)) * np.arctan((2 * u - cb) / (math.sqrt(3) * cb))
    dA = t ** (2 / 3) * (F(t0) - F(t))
    dL = onepz ** 2 * dA
    z = onepz - 1
    o = np.argsort(z)
    return np.interp(zq, z[o], dL[o])


def sn_chi2_from_dl(dl_hd):
    mu100 = 5 * np.log10((1 + SN.zhel) / (1 + SN.zhd) * (C_KMS / 100) * dl_hd) + 25
    y = SN.mb - np.where(SN.cal, SN.ceph, mu100)
    Wy = SN.W @ y
    b = SN.A.T @ Wy
    th = SN.Finv @ b
    return float(y @ Wy - b @ th), float(100 * 10 ** (th[1] / 5))


# checks: Einstein-de Sitter limit and low-z dressed Hubble constant
zt = np.array([0.1, 0.5, 1.0, 2.0])
eds = 2 * (1 + zt) * (1 - (1 + zt) ** -0.5)
print("\ntimescape check, f_v0 -> 0 vs Einstein-de Sitter:", np.round(timescape_dl(zt, 1e-6) / eds, 5))
for fv0 in (0.5, 0.7):
    slope = 1e-4 / timescape_dl(np.array([1e-4]), fv0)[0]
    print(f"  f_v0={fv0}: low-z dressed H0/H̄0 numeric {slope:.4f} vs formula {(4 * fv0 ** 2 + fv0 + 4) / (2 * (2 + fv0)):.4f}")

# ΛCDM SN-only check through the same generic function
Om_l = base["SN"][1][0]
E, Dc = background(UNI, Om_l)
c_l, _ = sn_chi2_from_dl((1 + SN.zhd) * np.interp(SN.zhd, Z_GRID, Dc))
print(f"  generic SN χ² for ΛCDM {c_l:.3f} vs class {SN.chi2(UNI, Om_l):.3f}")

fvs = np.linspace(0.30, 0.95, 66)
prof = [sn_chi2_from_dl(timescape_dl(SN.zhd, fv))[0] for fv in fvs]
i = int(np.argmin(prof))
r = optimize.minimize_scalar(lambda fv: sn_chi2_from_dl(timescape_dl(SN.zhd, fv))[0], bounds=(fvs[max(i - 1, 0)], fvs[min(i + 1, len(fvs) - 1)]), method="bounded")
fv_best = float(r.x)
c_ts, Hbare = sn_chi2_from_dl(timescape_dl(SN.zhd, fv_best))
inside = fvs[np.array(prof) <= c_ts + 1]
print(f"timescape, supernovae only: f_v0 = {fv_best:.3f} (1σ {inside.min():.2f}–{inside.max():.2f}), χ² = {c_ts:.2f}, "
      f"Δχ² vs ΛCDM = {c_ts - base['SN'][0]:+.2f}, bare H̄0 = {Hbare:.1f}, dressed H0 = {Hbare * (4 * fv_best ** 2 + fv_best + 4) / (2 * (2 + fv_best)):.1f}")
# same with a z > 0.075 Hubble-flow cut (calibrators kept)
keep = SN.cal | (SN.zhd > 0.075)


def sn_chi2_subset(dl_hd, mask):
    mu100 = 5 * np.log10((1 + SN.zhel) / (1 + SN.zhd) * (C_KMS / 100) * dl_hd) + 25
    y = (SN.mb - np.where(SN.cal, SN.ceph, mu100))
    idx = np.where(mask)[0]
    C = np.linalg.inv(SN.W)[np.ix_(idx, idx)]
    W = np.linalg.inv(C)
    A = SN.A[idx]
    yy = y[idx]
    Wy = W @ yy
    b = A.T @ Wy
    th = np.linalg.solve(A.T @ W @ A, b)
    return float(yy @ Wy - b @ th)


c_l75 = min(sn_chi2_subset((1 + SN.zhd) * np.interp(SN.zhd, Z_GRID, background(UNI, om)[1]), keep) for om in np.linspace(0.2, 0.45, 51))
c_t75 = min(sn_chi2_subset(timescape_dl(SN.zhd, fv), keep) for fv in np.linspace(0.4, 0.9, 51))
print(f"z > 0.075 Hubble flow ({int(keep.sum())} SNe incl. calibrators): Δχ² timescape − ΛCDM = {c_t75 - c_l75:+.2f}")


# ---------------------------------------------------------------- relative-density variant of the capped law
class RelWeb(CosmicWeb):
    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.x = np.exp(self.lnx_rel)          # density relative to the mean AT THAT TIME


RELWEB = RelWeb(R_mpc=2.0)


def total_chi2_law(law, web):
    lm = LawModel(law, web)
    r = optimize.minimize_scalar(lambda Om: sum(chi2_parts(lm, Om).values()), bounds=(0.1, 0.6), method="bounded", options=dict(xatol=1e-5))
    return float(r.fun), float(r.x), lm


c_base_all = base["SN+BAO+flows"][0]
rows = []
for xh in np.geomspace(0.5, 300, 22):
    for n in (0.5, 1.0, 2.0, 3.0):
        law = DarkEnergyLaw(x_half=float(xh), n=n, x_cut=200.0, p=1.0)
        for label, web in (("physical", WEB), ("relative", RELWEB)):
            c, Om, lm = total_chi2_law(law, web)
            rows.append(dict(variant=label, x_half=xh, n=n, d=c - c_base_all, void=float(law.shape(0.2) / law.shape(1.0) - 1),
                             w0=float(lm.w_eff(0.0)), kappa0=float(lm.kappa_web[0])))
tab = pd.DataFrame(rows)
for label in ("physical", "relative"):
    ok = tab[(tab.variant == label) & (tab.d < 4)]
    best = tab[tab.variant == label].sort_values("d").iloc[0]
    print(f"{label:8s} density: allowed laws {len(ok)}/{(tab.variant == label).sum()}, largest void excess allowed {100 * ok.void.max():.1f}%, "
          f"best Δχ² {best.d:+.2f} at x_half={best.x_half:.3g}, n={best.n}")
print(f"done in {time.time() - t_start:.0f}s")
