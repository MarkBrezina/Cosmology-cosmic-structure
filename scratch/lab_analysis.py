# Prototype of the notebook's analysis cells.
import time
from lab_core import *  # noqa: F401,F403

# ---------------------------------------------------------------- load
t0 = time.time()
pp = pd.read_csv(fetch("pantheon"), sep=r"\s+")
cov_flat = pd.read_csv(fetch("pantheon_cov"), header=None).to_numpy().ravel()
n_sn = int(cov_flat[0])
PP_COV = cov_flat[1:].reshape(n_sn, n_sn)
rows = [ln.split() for ln in open(fetch("desi_mean")) if ln.strip() and not ln.startswith("#")]
BAO_TABLE = pd.DataFrame(dict(z=[float(r[0]) for r in rows], value=[float(r[1]) for r in rows], kind=[r[2] for r in rows]))
BAO_COV = np.loadtxt(fetch("desi_cov"))
SPARC_DIR = get_sparc()
print(f"loaded in {time.time()-t0:.1f}s: {len(pp)} supernova rows, {len(BAO_TABLE)} BAO measurements")

t0 = time.time()
SN = Supernovae(pp, PP_COV)
BAO = BAOData(BAO_TABLE, BAO_COV)
WEB = CosmicWeb(R_mpc=2.0)
print(f"likelihoods + web built in {time.time()-t0:.1f}s; sigma(R=2 Mpc, z=0) = {WEB.sigma[0]:.2f}")


def fit_expansion(lm, use_sn=True, use_bao=True):
    def total(Om):
        return (SN.chi2(lm, Om) if use_sn else 0.0) + (BAO.chi2(lm, Om) if use_bao else 0.0)
    res = optimize.minimize_scalar(total, bounds=(0.1, 0.6), method="bounded", options=dict(xatol=1e-5))
    Om = float(res.x)
    c_sn, p_sn = SN.chi2(lm, Om, details=True)
    c_bao, p_bao = BAO.chi2(lm, Om, details=True)
    return dict(Om=Om, chi2_exp=float(res.fun), chi2_sn=c_sn, chi2_bao=c_bao, H0=p_sn["H0"], M=p_sn["M"], H0rd=p_bao["H0rd"])


t0 = time.time()
UNI = LawModel(UNIFORM, WEB)
base_sn = fit_expansion(UNI, use_bao=False)
base_bao = fit_expansion(UNI, use_sn=False)
base = fit_expansion(UNI)
print(f"uniform, supernovae only: Om={base_sn['Om']:.3f} H0={base_sn['H0']:.2f} chi2={base_sn['chi2_sn']:.1f} (N={SN.n})")
print(f"uniform, BAO only: Om={base_bao['Om']:.3f} chi2={base_bao['chi2_bao']:.2f} (N=13) H0rd={base_bao['H0rd']:.0f}")
print(f"uniform, both: Om={base['Om']:.4f} chi2_sn={base['chi2_sn']:.1f} chi2_bao={base['chi2_bao']:.2f}; {time.time()-t0:.1f}s")


# ---------------------------------------------------------------- measured velocity slope
def measure_velocity_slope(M, H0, Om, n_boot=1000):
    d = pp.copy()
    d["mu"] = d.m_b_corr - M
    d["w"] = 1 / d.m_b_corr_err_RAW ** 2
    g = d.groupby("CID", sort=False)
    u = g.first()
    u["mu"] = g.apply(lambda x: np.sum(x.w * x.mu) / np.sum(x.w), include_groups=False)
    u["err_raw"] = g.apply(lambda x: 1 / np.sqrt(np.sum(x.w)), include_groups=False)
    u = u[(u.zCMB > 0.01) & (u.zCMB < 0.06)].copy()
    E, Dc = background(UNI, Om)
    chi = 10 ** ((u.mu - 25) / 5) / (1 + u.zHEL)
    zcos = np.interp(chi * H0 / C_KMS, Dc, Z_GRID)
    v = C_KMS * ((1 + u.zCMB) / (1 + zcos) - 1)
    sig_mu = np.sqrt(u.err_raw ** 2 + 0.103 ** 2)
    sv = C_KMS * (1 + u.zCMB) / (1 + zcos) ** 2 * chi * LN10_5 * sig_mu * H0 * np.interp(zcos, Z_GRID, E) / C_KMS
    x, y, s = u.VPEC.to_numpy(), v.to_numpy(), sv.to_numpy()

    def nll(th, x, y, s):
        b0, k, lns = th
        var = s ** 2 + np.exp(2 * lns)
        return 0.5 * np.sum((y - b0 - k * x) ** 2 / var + np.log(var))
    best = optimize.minimize(nll, [0.0, 1.0, math.log(280.0)], args=(x, y, s), method="Nelder-Mead").x
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, len(x), len(x))
        boots.append(optimize.minimize(nll, best, args=(x[i], y[i], s[i]), method="Nelder-Mead").x[1])
    return dict(k=float(best[1]), k_err=float(np.std(boots)), sigma_v=float(np.exp(best[2])), n=len(x))


t0 = time.time()
VEL = measure_velocity_slope(base_sn["M"], base_sn["H0"], base_sn["Om"])
print(f"velocity slope {VEL['k']:.3f} +/- {VEL['k_err']:.3f}, random motion {VEL['sigma_v']:.0f} km/s, {VEL['n']} supernovae; {time.time()-t0:.1f}s")

f_ref, D_ref = growth(UNI, 0.315)
REF_FD = f_ref * D_ref
f_b, D_b = growth(UNI, base["Om"])
print(f"reference f={f_ref:.3f}; uniform at fitted Om: predicted slope {f_b*D_b/REF_FD:.3f}")


def evaluate(law, fit_to=("sn", "bao")):
    lm = LawModel(law, WEB)
    fit = fit_expansion(lm, use_sn="sn" in fit_to, use_bao="bao" in fit_to)
    f, D = growth(lm, fit["Om"])
    k_pred = f * D / REF_FD
    chi2_vel = ((k_pred - VEL["k"]) / VEL["k_err"]) ** 2
    out = dict(law=law.name, x_half=law.x_half, n=law.n, x_cut=law.x_cut, p=law.p, **fit,
               k_pred=k_pred, chi2_vel=chi2_vel, kappa0=float(lm.kappa_web[0]),
               w0=float(lm.w_eff(0.0)), w05=float(lm.w_eff(0.5)), w1=float(lm.w_eff(1.0)),
               de_z1=float(lm.f(1.0)), de_z2=float(lm.f(2.0)))
    return out, lm


t0 = time.time()
BASE, _ = evaluate(UNIFORM)
print(f"evaluate(uniform) in {time.time()-t0:.2f}s:", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in BASE.items()})
for law in [DarkEnergyLaw(), DarkEnergyLaw(x_half=10), DarkEnergyLaw(x_half=100), DarkEnergyLaw(x_half=300, n=2), DarkEnergyLaw(x_half=5, n=1)]:
    t0 = time.time()
    r, lm = evaluate(law)
    d_exp = r["chi2_exp"] - BASE["chi2_exp"]
    d_vel = r["chi2_vel"] - BASE["chi2_vel"]
    print(f"x_half={law.x_half:6.1f} n={law.n:.1f} cut={law.x_cut:.0f}: Om={r['Om']:.3f} d_exp={d_exp:7.2f} d_vel={d_vel:6.2f} "
          f"k={r['k_pred']:.3f} kappa0={r['kappa0']:.3f} w0={r['w0']:.3f} w1={r['w1']:.3f} DE(z=1)/DE0={r['de_z1']:.3f} meanS0={lm.mean_S0:.3f} ({time.time()-t0:.2f}s)")
