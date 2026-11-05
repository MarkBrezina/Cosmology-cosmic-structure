"""Round 2: CMB distance priors + Mark's surface-growth laws, scored with SN + BAO + flows (+ CMB)."""
import math
import sys
import time

sys.path.insert(0, "/home/claude/lab")
from lab_core import *  # noqa: F401,F403

T0 = time.time()
pp = pd.read_csv(fetch("pantheon"), sep=r"\s+")
_c = pd.read_csv(fetch("pantheon_cov"), header=None).to_numpy().ravel()
PP_COV = _c[1:].reshape(int(_c[0]), int(_c[0]))
_rows = [ln.split() for ln in open(fetch("desi_mean")) if ln.strip() and not ln.startswith("#")]
BAO_TABLE = pd.DataFrame(dict(z=[float(r[0]) for r in _rows], value=[float(r[1]) for r in _rows], kind=[r[2] for r in _rows]))
BAO_COV = np.loadtxt(fetch("desi_cov"))
SN, BAO, WEB = Supernovae(pp, PP_COV), BAOData(BAO_TABLE, BAO_COV), CosmicWeb(R_mpc=2.0)
UNI = LawModel(UNIFORM, WEB)


# ---------------------------------------------------------------- histories
class History:
    """Any dark-energy history: f(z) = rho_DE(z)/rho_DE(0) on a grid, g(z) = kappa*f (0 = smooth in space)."""
    def __init__(self, z, f, g=None, name=""):
        self._lz = np.log1p(z)
        self.f_web = np.maximum(np.asarray(f, float), 0.0)
        self.g_web = np.zeros_like(self.f_web) if g is None else np.asarray(g, float)
        self.f_grid = self.f(Z_GRID)
        self.name = name

    def f(self, z):
        return np.interp(np.log1p(z), self._lz, self.f_web)

    def g(self, z):
        return np.interp(np.log1p(z), self._lz, self.g_web)

    def w_eff(self, z):
        dlnf = np.gradient(np.log(np.maximum(self.f_web, 1e-300)), self._lz)
        return -1 + np.interp(np.log1p(z), self._lz, dlnf) / 3


ZH = np.concatenate([np.linspace(0, 3, 301), np.geomspace(3, 1e8, 400)[1:]])


def cpl(w0, wa):
    return History(ZH, (1 + ZH) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * ZH / (1 + ZH)), name=f"w0={w0:.3f}, wa={wa:.3f}")


def from_lawmodel(lm, smooth=False):
    z = lm.web_z
    zz = np.concatenate([z, np.geomspace(z[-1] * 1.01, 1e8, 50)])
    f = np.concatenate([lm.f_web, np.full(50, lm.f_web[-1])])
    g = None if smooth else np.concatenate([lm.g_web, np.full(50, lm.g_web[-1])])
    return History(zz, f, g)


D_LCDM = growth_lcdm(ZH[ZH < 2e4])                      # web growth history, computed once
_ZD = ZH[ZH < 2e4]


def void_surface(x_void, R_mpc):
    """Mark's surface law: void volume grows through void boundary area, dV_void/dt = alpha * S_void.
    The boundary of 'genuinely empty' regions (density below x_void times the mean) per unit volume
    follows the Rice formula for the lognormal web, S_V ∝ exp(-nu^2/2) per comoving volume, nu set by
    how far the threshold sits in the distribution. Physical area per volume = comoving × (1+z).
    H_DE ∝ S_V,phys, so rho_DE ∝ S_V,phys^2."""
    sig = sigma_tophat(R_mpc) * D_LCDM
    s = np.sqrt(np.log1p(sig ** 2))
    nu = (math.log(x_void) + s ** 2 / 2) / s
    lnSV = -nu ** 2 / 2 + np.log1p(_ZD)
    f = np.exp(2 * (lnSV - lnSV[0]))
    zz = np.concatenate([_ZD, np.geomspace(_ZD[-1] * 1.01, 1e8, 30)])
    return History(zz, np.concatenate([f, np.zeros(30)]), name=f"void surface, x_void={x_void:.3g}, R={R_mpc:.3g}")


# ---------------------------------------------------------------- CMB distance priors (Planck 2018, Zhai & Wang 2019)
OMEGA_GAMMA = 2.469e-5
OMEGA_R = OMEGA_GAMMA * (1 + 0.2271 * 3.046)
CMB_MEAN = np.array([1.74963, 301.80845, 0.02237])
CMB_COV = 1e-8 * np.array([[1598.9554, 17112.007, -36.311179],
                           [17112.007, 811208.45, -494.79813],
                           [-36.311179, -494.79813, 2.1242182]])
CMB_W = np.linalg.inv(CMB_COV)


def z_star(wb, wm):
    g1 = 0.0783 * wb ** -0.238 / (1 + 39.5 * wb ** 0.763)
    g2 = 0.560 / (1 + 21.1 * wb ** 1.81)
    return 1048 * (1 + 0.00124 * wb ** -0.738) * (1 + g1 * wm ** g2)


_X1 = np.linspace(0.0, 1.0, 3001)          # fractions of ln(1+z*) for the distance integral
_X2 = np.linspace(0.0, 1.0, 3001)


def _cmb_raw(hist, Om, h, wb):
    Or = OMEGA_R / h ** 2
    zs = z_star(wb, Om * h * h)
    lzs = math.log1p(zs)

    def E(z):
        return np.sqrt(Om * (1 + z) ** 3 + Or * (1 + z) ** 4 + (1 - Om - Or) * hist.f(z))
    x = _X1 * lzs
    z = np.expm1(x)
    DM = integrate.trapezoid((1 + z) / E(z), x)                       # c/H0 units
    x2 = lzs + _X2 * (math.log(1e8) - lzs)
    z2 = np.expm1(x2)
    Rb = 3 * wb / (4 * OMEGA_GAMMA * (1 + z2))
    rs = integrate.trapezoid((1 + z2) / (E(z2) * np.sqrt(3 * (1 + Rb))), x2)
    return np.array([math.sqrt(Om) * DM, math.pi * DM / rs, wb])


# The Hu & Sugiyama fitting formula puts z* at 1091.9 for Planck's best fit, where Planck's full
# recombination code gives 1089.92; that alone moves l_A by 0.47 (5σ). Calibrate once so Planck's
# best-fit ΛCDM reproduces the prior means exactly; model-to-model differences come from D_M(z*).
CMB_CAL = CMB_MEAN / _cmb_raw(UNI, 0.3153, 0.6736, 0.02237)


def cmb_vector(hist, Om, h, wb):
    return _cmb_raw(hist, Om, h, wb) * CMB_CAL


def cmb_chi2(hist, Om, h, wb):
    d = cmb_vector(hist, Om, h, wb) - CMB_MEAN
    return float(d @ CMB_W @ d)


print("CMB calibration factors (R, l_A, ωb):", np.round(CMB_CAL, 6))
_r = optimize.minimize(lambda th: cmb_chi2(UNI, *th) if 0.6 < th[1] < 0.8 else 1e9, [0.30, 0.69, 0.0222], method="Nelder-Mead",
                       options=dict(xatol=1e-9, fatol=1e-10, maxiter=4000))
print(f"ΛCDM from the CMB priors alone: Ωm={_r.x[0]:.4f} h={_r.x[1]:.4f} ωb={_r.x[2]:.5f} χ²={_r.fun:.4f}  (Planck 2018: 0.3153, 0.6736, 0.02237)")


# ---------------------------------------------------------------- flows
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


_r = optimize.minimize_scalar(lambda Om: SN.chi2(UNI, Om), bounds=(0.1, 0.6), method="bounded")
_c, _p = SN.chi2(UNI, float(_r.x), details=True)
VEL = measure_velocity_slope(_p["M"], _p["H0"], float(_r.x))
_f, _D = growth(UNI, 0.315)
REF_FD = _f * _D
print(f"flow slope {VEL['k']:.3f} ± {VEL['k_err']:.3f}")


def parts(hist, Om, h=None, wb=None, use=("sn", "bao", "vel", "cmb")):
    out = {}
    if "sn" in use:
        out["sn"] = SN.chi2(hist, Om)
    if "bao" in use:
        out["bao"] = BAO.chi2(hist, Om)
    if "vel" in use:
        f, D = growth(hist, Om)
        out["vel"] = ((f * D / REF_FD - VEL["k"]) / VEL["k_err"]) ** 2
        out["k_pred"] = f * D / REF_FD
    if "cmb" in use:
        out["cmb"] = cmb_chi2(hist, Om, h, wb)
    return out


def fit_family(make, starts, bounds, use):
    """Fit Ωm (+ h, ωb when the CMB is used) and the family's own parameters.
    Stage 1 fits without the flows term (fast); stage 2 polishes with everything included."""
    with_cmb = "cmb" in use
    use_all = use
    use = tuple(u for u in use_all if u != "vel")

    def obj(th):
        Om, rest = th[0], th[1:]
        if with_cmb:
            h, wb, rest = rest[0], rest[1], rest[2:]
            if not (0.5 < h < 0.9 and 0.019 < wb < 0.026):
                return 1e9
        else:
            h = wb = None
        if not (0.05 < Om < 0.7) or any(not (lo <= v <= hi) for v, (lo, hi) in zip(rest, bounds)):
            return 1e9
        try:
            hist = make(*rest)
        except Exception:
            return 1e9
        p = parts(hist, Om, h, wb, use)
        return sum(v for k, v in p.items() if k != "k_pred")
    best = None
    for s in starts:
        x0 = [s[0]] + ([0.674, 0.02237] if with_cmb else []) + list(s[1:])
        r = optimize.minimize(obj, x0, method="Nelder-Mead",
                              options=dict(xatol=1e-6, fatol=1e-6, maxiter=3000, maxfev=3000, adaptive=True))
        if best is None or r.fun < best.fun:
            best = r
    if "vel" in use_all:
        use = use_all
        best = optimize.minimize(obj, best.x, method="Nelder-Mead",
                                 options=dict(xatol=1e-6, fatol=1e-6, maxiter=1200, maxfev=1200, adaptive=True))
    th = best.x
    Om = th[0]
    h, wb, rest = (th[1], th[2], th[3:]) if with_cmb else (None, None, th[1:])
    hist = make(*rest)
    p = parts(hist, Om, h, wb, use)
    return dict(chi2=best.fun, Om=Om, h=h, wb=wb, params=list(rest), parts=p, hist=hist)


FAMILIES = {
    "ΛCDM = your volume compounding (p = 1)": (lambda: UNI, [[0.31]], []),
    "Smooth surface growth, whole universe (p = 2/3)": (lambda: cpl(-1 / 3, 0.0), [[0.3]], []),
    "Rough surface of dimension D (p = D/3)": (lambda D: cpl(1 - 2 * D / 3, 0.0), [[0.3, 2.9], [0.31, 3.0]], [(2.0, 3.6)]),
    "DESI's evolving dark energy (w0, wa)": (lambda w0, wa: cpl(w0, wa) if w0 + wa < 0 else (_ for _ in ()).throw(ValueError()),
                                             [[0.31, -0.9, -0.3], [0.31, -0.75, -0.8], [0.31, -1.0, 0.0]], [(-2.0, 0.0), (-4.0, 2.0)]),
    "Void-surface law (empty below x_void, web cell R)": (lambda lx, lR: void_surface(10 ** lx, 10 ** lR),
                                                          [[0.31, -1.0, 0.3], [0.31, -1.5, 0.7], [0.31, -0.7, 1.0], [0.31, -2.0, 0.0]],
                                                          [(-3.0, -0.05), (-0.3, 1.5)]),
}

MARK_LOCAL = from_lawmodel(LawModel(DarkEnergyLaw(x_half=30.0, n=1.0, x_cut=200.0, p=1.0), WEB))


def switch_on(lxh, n, p, smooth=True):
    return from_lawmodel(LawModel(DarkEnergyLaw(x_half=10 ** lxh, n=n, x_cut=200.0, p=p), WEB), smooth=smooth)


FAMILIES["Your capped law as specified (local density)"] = (lambda: MARK_LOCAL, [[0.31]], [])
FAMILIES["Switch-on by average density + dilution p"] = (lambda lxh, n, p: switch_on(lxh, n, p, True),
                                                         [[0.31, 0.0, 0.6, 0.83], [0.31, 0.5, 1.0, 0.9], [0.31, 1.5, 2.0, 0.95]],
                                                         [(-1.0, 4.0), (0.3, 6.0), (0.6, 1.3)])

COMBOS = {"SN+BAO+flows": ("sn", "bao", "vel"), "SN+BAO+CMB+flows": ("sn", "bao", "cmb", "vel")}
res = {}
for cname, use in COMBOS.items():
    print(f"\n=== {cname}")
    for name, (make, starts, bounds) in FAMILIES.items():
        t = time.time()
        r = fit_family(make, starts, bounds, use)
        res[(cname, name)] = r
        base = res[(cname, "ΛCDM = your volume compounding (p = 1)")]["chi2"]
        p = r["parts"]
        extra = f" h={r['h']:.4f} ωb={r['wb']:.5f}" if r["h"] else ""
        print(f"{name:52s} Δχ²={r['chi2'] - base:+8.2f}  Ωm={r['Om']:.4f}{extra}  params={np.round(r['params'], 4)}  "
              f"flow slope={p['k_pred']:.2f}  [" + ", ".join(f"{k}={v:.2f}" for k, v in p.items() if k != 'k_pred') + f"]  ({time.time() - t:.0f}s)")
        hist = r["hist"]
        if hasattr(hist, "w_eff"):
            print(f"{'':52s} w(0)={float(hist.w_eff(0.0)):+.3f}  w(0.5)={float(hist.w_eff(0.5)):+.3f}  w(1)={float(hist.w_eff(1.0)):+.3f}  "
                  f"ρDE(z=0.5)/ρDE0={float(hist.f(0.5)):.3f}  ρDE(1)/ρDE0={float(hist.f(1.0)):.3f}  ρDE(2)/ρDE0={float(hist.f(2.0)):.3f}")
print(f"\ntotal {time.time() - T0:.0f}s")
