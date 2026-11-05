# Prototype of the notebook engine (cells are assembled from these blocks).
# ---------------------------------------------------------------- setup
import hashlib
import math
import shutil
import urllib.request
import zipfile
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import integrate, optimize, special

DATA = Path("data")
OUT = Path("lab_outputs")
DATA.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)

C_KMS = 299792.458
G_SI = 6.67430e-11
MPC_M = 3.0856776e22
KPC_M = 3.0856776e19
M_P = 1.67262192e-27
LN10_5 = math.log(10) / 5
rng = np.random.default_rng(2026)

# ---------------------------------------------------------------- data
RAW_PP = "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/"
RAW_BAO = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/"
SOURCES = {
    "pantheon": dict(path=DATA / "pantheon" / "Pantheon+SH0ES.dat",
                     urls=[RAW_PP + "Pantheon%2BSH0ES.dat"],
                     sha256="1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8"),
    "pantheon_cov": dict(path=DATA / "pantheon" / "Pantheon+SH0ES_STAT+SYS.cov",
                         urls=[RAW_PP + "Pantheon%2BSH0ES_STAT%2BSYS.cov"],
                         sha256="abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc"),
    "desi_mean": dict(path=DATA / "desi" / "desi_gaussian_bao_ALL_GCcomb_mean.txt",
                      urls=[RAW_BAO + "desi_gaussian_bao_ALL_GCcomb_mean.txt"],
                      sha256="9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585"),
    "desi_cov": dict(path=DATA / "desi" / "desi_gaussian_bao_ALL_GCcomb_cov.txt",
                     urls=[RAW_BAO + "desi_gaussian_bao_ALL_GCcomb_cov.txt"],
                     sha256="252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509"),
    "sparc_zip": dict(path=DATA / "sparc" / "Rotmod_LTG.zip",
                      urls=["https://astroweb.case.edu/SPARC/Rotmod_LTG.zip",
                            "https://zenodo.org/records/16284118/files/Rotmod_LTG.zip?download=1"],
                      sha256="0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588"),
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(key, required=True):
    """Download one data file into data/ unless it is already there, then check its fingerprint."""
    src = SOURCES[key]
    path = src["path"]
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        for url in src["urls"]:
            try:
                print(f"downloading {path.name} from {url.split('/')[2]} ...")
                req = urllib.request.Request(url, headers={"User-Agent": "dark-energy-law-lab"})
                with urllib.request.urlopen(req, timeout=300) as r, open(path, "wb") as fh:
                    shutil.copyfileobj(r, fh)
                break
            except Exception as err:  # try the next mirror
                print(f"  failed: {err}")
                if path.exists():
                    path.unlink()
        else:
            msg = f"Could not download {path.name}. Save it by hand as {path} and re-run this cell."
            if required:
                raise RuntimeError(msg)
            print(msg)
            return None
    same = sha256(path) == src["sha256"]
    print(f"{path.name}: {'fingerprint matches the reference copy' if same else 'present (fingerprint differs from the reference copy)'}")
    return path


def get_sparc():
    """SPARC rotation curves (optional, only used for the galaxy check)."""
    folder = DATA / "sparc" / "Rotmod_LTG"
    files = list(folder.glob("*_rotmod.dat")) if folder.exists() else []
    if len(files) >= 170:
        print(f"SPARC rotation curves: {len(files)} galaxies found in {folder}")
        return folder
    z = fetch("sparc_zip", required=False)
    if z is None:
        return None
    folder.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(z) as zz:
        for member in zz.namelist():
            if member.endswith("_rotmod.dat"):
                (folder / Path(member).name).write_bytes(zz.read(member))
    print(f"SPARC rotation curves: unpacked {len(list(folder.glob('*_rotmod.dat')))} galaxies")
    return folder


# ---------------------------------------------------------------- the law
@dataclass
class DarkEnergyLaw:
    """How much dark energy a patch of space carries, as a function of its matter density.

    Density x is measured in units of TODAY's average cosmic matter density
    (x = 1 is an average patch today; deep voids 0.1-0.3; walls and filaments 3-30;
    galaxy clusters ~200-1000; a galaxy's disc ~1e6). In the past the whole universe was
    denser, x_mean = (1 + z)^3, so the law also decides how dark energy switched on over time.

    shape(x) is the fraction of the maximum (the cap) that a patch carries:
      * shape(0) = 1: perfectly empty space carries the most dark energy (the cap),
      * it falls as density rises (the inverse part, controlled by x_half and n),
      * above x_cut it stops altogether (the cutoff).
    p is the compounding exponent of your growth rule dV/dt ∝ V^p (p = 1: in proportion to volume).
    """
    name: str = "capped inverse law"
    x_half: float = 30.0         # density where dark energy has fallen to half its cap
    n: float = 1.0               # sharpness of the fall: 1 = plain inverse, larger = sharper
    x_cut: float = 200.0         # above this density dark energy stops completely
    cut_width_dex: float = 0.15  # softness of the cutoff in decades of density (0 = hard switch)
    p: float = 1.0               # compounding exponent of the volume-growth rule

    def shape(self, x):
        x = np.maximum(np.asarray(x, dtype=float), 0.0)
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            s = 1.0 / (1.0 + (x / self.x_half) ** self.n)      # x_half = inf gives 1 everywhere
        if np.isfinite(self.x_cut):
            if self.cut_width_dex <= 0:
                s = np.where(x < self.x_cut, s, 0.0)
            else:
                t = np.log10(np.maximum(x, 1e-300) / self.x_cut) / self.cut_width_dex
                s = s * 0.5 * special.erfc(t / math.sqrt(2.0))
        return s


UNIFORM = DarkEnergyLaw(name="uniform dark energy (standard model)", x_half=math.inf, x_cut=math.inf)


# ---------------------------------------------------------------- the cosmic web
def _tophat(y):
    y = np.asarray(y, float)
    small = y < 1e-2
    out = np.empty_like(y)
    out[small] = 1 - y[small] ** 2 / 10
    ys = y[~small]
    out[~small] = 3 * (np.sin(ys) - ys * np.cos(ys)) / ys ** 3
    return out


def sigma_tophat(R_mpc, sigma8=0.81, Om=0.315, h=0.674, ns=0.965, Ob=0.049):
    """Linear rms density contrast in spheres of radius R (BBKS transfer function)."""
    gamma = Om * h * math.exp(-Ob * (1 + math.sqrt(2 * h) / Om))
    k = np.logspace(-4, 3, 8000)                                    # h/Mpc
    q = k / gamma
    T = np.log(1 + 2.34 * q) / (2.34 * q) * (1 + 3.89 * q + (16.1 * q) ** 2 + (5.46 * q) ** 3 + (6.71 * q) ** 4) ** -0.25
    P = k ** ns * T ** 2

    def var(Rh):
        return integrate.trapezoid(k ** 2 * P * _tophat(k * Rh) ** 2, k) / (2 * math.pi ** 2)
    return sigma8 * math.sqrt(var(R_mpc * h) / var(8.0))


def growth_lcdm(z, Om=0.315):
    """Standard-model linear growth factor, D(z=0) = 1."""
    def E(a):
        return math.sqrt(Om * a ** -3 + 1 - Om)

    def D(a):
        return 2.5 * Om * E(a) * integrate.quad(lambda aa: 1 / (aa * E(aa)) ** 3, 0, a, limit=200)[0]
    d0 = D(1.0)
    return np.array([D(1 / (1 + zi)) for zi in np.atleast_1d(z)]) / d0


class CosmicWeb:
    """How the universe's volume is shared out between densities, at each redshift.

    Matter is smoothed in cells of radius R_mpc and the cell densities follow a lognormal
    distribution whose width grows with standard structure growth (Coles & Jones 1991).
    Averages are volume averages, which is what sets the mean dark-energy density.
    """
    def __init__(self, R_mpc=2.0, sigma8=0.81, nodes=120):
        self.R_mpc = R_mpc
        self.z = np.concatenate([np.linspace(0, 3, 151), np.geomspace(3, 1100, 80)[1:]])
        self.sigma = sigma_tophat(R_mpc, sigma8) * growth_lcdm(self.z)
        s2 = np.log1p(self.sigma ** 2)
        mu, s = -0.5 * s2, np.sqrt(s2)
        t, w = np.polynomial.hermite.hermgauss(nodes)
        self.wts = w / math.sqrt(math.pi)
        self.lnx_rel = mu[:, None] + math.sqrt(2.0) * s[:, None] * t[None, :]   # ln(density / mean at that time)
        self.x = (1 + self.z)[:, None] ** 3 * np.exp(self.lnx_rel)             # density / today's mean

    def averages(self, law):
        """Volume averages at each redshift: <S> (mean dark energy relative to the cap) and
        -<x dS/dx> (how much the mean responds when a large region is a little denser)."""
        S = law.shape(self.x)
        mean_S = S @ self.wts
        h = 1e-3
        xdS = (law.shape(self.x * math.exp(h)) - law.shape(self.x * math.exp(-h))) / (2 * h)
        return mean_S, -(xdS @ self.wts)

    def volume_fraction(self, lo, hi, z=0.0):
        """Share of the volume with density (relative to today's mean) between lo and hi."""
        i = int(np.argmin(np.abs(self.z - z)))
        x = self.x[i]
        return float(self.wts[(x >= lo) & (x < hi)].sum())


# ---------------------------------------------------------------- background expansion
Z_GRID = np.concatenate([np.linspace(0, 0.2, 801), np.linspace(0.2, 2.6, 2401)[1:]])


class LawModel:
    """Everything about a law that does not depend on the matter density Om."""
    def __init__(self, law, web):
        self.law = law
        mean_S, resp = web.averages(law)
        if not mean_S[0] > 0:
            raise ValueError("this law leaves no dark energy anywhere today; raise x_half or x_cut")
        self.mean_S0 = float(mean_S[0])
        cap = (1 + web.z) ** (6 * (1 - law.p))                 # how the cap itself evolves (p = 1: constant)
        self.web_z = web.z
        self.f_web = cap * mean_S / self.mean_S0              # mean dark energy relative to today
        self.g_web = cap * resp / self.mean_S0                 # kappa(z) * f(z), the term that feeds growth
        with np.errstate(invalid="ignore", divide="ignore"):
            self.kappa_web = np.where(mean_S > 1e-12, resp / mean_S, 0.0)
        self._lz = np.log1p(web.z)
        self.f_grid = self.f(Z_GRID)

    def f(self, z):
        return np.interp(np.log1p(z), self._lz, self.f_web)

    def g(self, z):
        return np.interp(np.log1p(z), self._lz, self.g_web)

    def kappa(self, z):
        """How strongly large-scale dark energy follows large-scale matter: d rho_DE/rho_DE = -kappa d rho_m/rho_m."""
        return np.interp(np.log1p(z), self._lz, self.kappa_web)

    def w_eff(self, z):
        """Equation of state a constant-w fluid would need to mimic this dark energy at redshift z."""
        lf = np.log(np.maximum(self.f_web, 1e-300))
        dlnf = np.gradient(lf, self._lz)
        return -1 + np.interp(np.log1p(z), self._lz, dlnf) / 3


def background(lm, Om):
    E = np.sqrt(Om * (1 + Z_GRID) ** 3 + (1 - Om) * lm.f_grid)
    Dc = integrate.cumulative_trapezoid(1 / E, Z_GRID, initial=0.0)   # comoving distance in units of c/H0
    return E, Dc


# ---------------------------------------------------------------- supernovae
class Supernovae:
    """Pantheon+SH0ES likelihood with its full covariance. The supernova absolute magnitude M and
    H0 enter linearly, so they are solved for exactly at every step (Cepheid hosts anchor M)."""
    def __init__(self, table, cov):
        use = ((table.IS_CALIBRATOR == 1) | (table.zHD > 0.01)).to_numpy()
        self.idx = np.where(use)[0]
        d = table.iloc[self.idx]
        self.zhd, self.zhel = d.zHD.to_numpy(), d.zHEL.to_numpy()
        self.mb, self.ceph = d.m_b_corr.to_numpy(), d.CEPH_DIST.to_numpy()
        self.cal = (d.IS_CALIBRATOR == 1).to_numpy()
        self.W = np.linalg.inv(cov[np.ix_(self.idx, self.idx)])
        self.A = np.column_stack([np.ones(len(self.idx)), np.where(self.cal, 0.0, -1.0)])
        self.Finv = np.linalg.inv(self.A.T @ self.W @ self.A)
        self.n = len(self.idx)

    def chi2(self, lm, Om, details=False):
        _, Dc = background(lm, Om)
        dc = np.interp(self.zhd, Z_GRID, Dc)
        mu100 = 5 * np.log10((1 + self.zhel) * (C_KMS / 100) * dc) + 25     # distance modulus for H0 = 100
        y = self.mb - np.where(self.cal, self.ceph, mu100)
        Wy = self.W @ y
        b = self.A.T @ Wy
        theta = self.Finv @ b
        chi2 = float(y @ Wy - b @ theta)
        if details:
            return chi2, dict(M=float(theta[0]), H0=float(100 * 10 ** (theta[1] / 5)))
        return chi2


# ---------------------------------------------------------------- DESI BAO
class BAOData:
    """DESI DR2 baryon acoustic oscillations. Their absolute scale c/(H0 r_d) is solved for exactly,
    so BAO tests only the SHAPE of the expansion history (no conflict with the H0 tension)."""
    def __init__(self, table, cov):
        self.z, self.d, self.kind = table.z.to_numpy(), table.value.to_numpy(), table.kind.to_numpy()
        self.W = np.linalg.inv(cov)

    def u(self, lm, Om):
        E, Dc = background(lm, Om)
        dm = np.interp(self.z, Z_GRID, Dc)
        dh = 1 / np.interp(self.z, Z_GRID, E)
        dv = (self.z * dm ** 2 * dh) ** (1 / 3)
        return np.select([self.kind == "DM_over_rs", self.kind == "DH_over_rs"], [dm, dh], dv)

    def chi2(self, lm, Om, details=False):
        u = self.u(lm, Om)
        Wu = self.W @ u
        amp = float(Wu @ self.d) / float(Wu @ u)
        r = self.d - amp * u
        chi2 = float(r @ self.W @ r)
        if details:
            return chi2, dict(H0rd=C_KMS / amp, bao_ratio=self.d / (amp * u))
        return chi2


# ---------------------------------------------------------------- growth of structure
def growth(lm, Om):
    """Linear growth with the law's background plus the extra pull of dark-energy lumps.
    Vacuum energy gravitates with -2 x its density, so dark energy that follows matter
    (kappa > 0) adds 2 kappa rho_DE to the source of gravity for density contrasts.
    Returns (f, D) today, with D = a early on so it can be compared between models."""
    lna = np.linspace(math.log(1e-3), 0.0, 1200)
    a = np.exp(lna)
    z = 1 / a - 1
    fz = lm.f(z)
    E2 = Om * a ** -3 + (1 - Om) * fz
    dlnE = np.gradient(0.5 * np.log(E2), lna)
    Om_a = Om * a ** -3 / E2
    geff = 1 + 2 * lm.g(z) * (1 - Om) / (Om * a ** -3)

    def rhs(x, y):
        dl = np.interp(x, lna, dlnE)
        src = 1.5 * np.interp(x, lna, Om_a) * np.interp(x, lna, geff)
        return [y[1], -(2 + dl) * y[1] + src * y[0]]
    sol = integrate.solve_ivp(rhs, (lna[0], 0.0), [a[0], a[0]], rtol=1e-8, atol=1e-11)
    D, dD = sol.y[0, -1], sol.y[1, -1]
    return dD / D, D
