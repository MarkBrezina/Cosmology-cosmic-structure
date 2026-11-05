#!/usr/bin/env python3
"""
Expanding-voids hypothesis: quantitative tests
===============================================

Hypothesis (Mark): every unit of empty space expands at a rate proportional to the amount of
empty space already there (dV/dt = alpha * V). Inside galaxies this should change how fast the
outer, low-density regions rotate; between galaxies it should make the local Hubble rate vary
from place to place while averaging to the global value.

Part A  Cosmology   - fit flat LCDM / wCDM to Pantheon+ supernovae (SH0ES-calibrated) and turn
                      the cosmological constant into "volume growth of empty space" (alpha).
Part B  Galaxies    - test expansion-type extra accelerations against 175 SPARC rotation curves,
                      next to an acceleration law (RAR) and dark-matter haloes.
Part C  Hubble flow - how much the local expansion rate varies between supernova host galaxies,
                      and how much of that ordinary gravity (2M++ density field) already predicts.

Data (public):
  SPARC rotation curves, Lelli, McGaugh & Schombert (2016), AJ 152, 157   (Rotmod_LTG/*.dat)
  SPARC BTFR table,      Lelli, McGaugh & Schombert (2016), ApJL 816, L14 (BTFR_Lelli2016a.mrt)
  Pantheon+SH0ES,        Scolnic et al. (2022), Brout et al. (2022)        (Pantheon+SH0ES.dat + cov)

Run:
  python -I expanding_voids_analysis.py --sparc DIR --btfr FILE --pantheon DIR --out results.json
"""
import argparse
import glob
import json
import math
import os

import numpy as np
import pandas as pd
from scipy import integrate, optimize, stats

# ----------------------------------------------------------------------------- constants
C_KMS = 299792.458
G_SI = 6.67430e-11
KPC_M = 3.0856776e19
MPC_M = 3.0856776e22
MSUN = 1.98847e30
M_P = 1.67262192e-27
YR_S = 3.15576e7
GYR_S = 3.15576e16
GEV_CM3 = 1.78266192e-21          # 1 GeV/cm^3 in kg/m^3
CONV = 1e6 / KPC_M                # (km/s)^2/kpc  ->  m/s^2
UPS_DISK, UPS_BUL = 0.5, 0.7      # 3.6 micron mass-to-light ratios (Lelli+16, McGaugh+16)
LN10_5 = math.log(10) / 5.0

rng = np.random.default_rng(20261026)


def h_si(H_kms_mpc):
    """km/s/Mpc -> 1/s"""
    return H_kms_mpc * 1e3 / MPC_M


def r3(x, n=3):
    return float(f"{x:.{n}g}")


# ============================================================================= PART A
def comoving_grid(Om, w=-1.0, zmax=2.5):
    zg = np.concatenate([np.linspace(0, 0.2, 4001), np.linspace(0.2, zmax, 3001)[1:]])
    E = np.sqrt(Om * (1 + zg) ** 3 + (1 - Om) * (1 + zg) ** (3 * (1 + w)))
    dc = integrate.cumulative_trapezoid(1.0 / E, zg, initial=0.0)   # D_C * H0 / c
    return zg, dc, E


def mu_theory(zhd, zhel, H0, Om, w=-1.0):
    zg, dc, _ = comoving_grid(Om, w)
    d_l = (1 + zhel) * np.interp(zhd, zg, dc) * C_KMS / H0
    return 5 * np.log10(d_l) + 25


def numeric_hessian(f, p, steps):
    p = np.asarray(p, float)
    n = len(p)
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            ei = np.zeros(n); ej = np.zeros(n)
            ei[i] = steps[i]; ej[j] = steps[j]
            fpp = f(p + ei + ej); fpm = f(p + ei - ej)
            fmp = f(p - ei + ej); fmm = f(p - ei - ej)
            H[i, j] = H[j, i] = (fpp - fpm - fmp + fmm) / (4 * steps[i] * steps[j])
    return H


def load_pantheon(pdir):
    dat = pd.read_csv(os.path.join(pdir, "PantheonPlusSH0ES.dat"), sep=r"\s+")
    cov_flat = np.loadtxt(os.path.join(pdir, "PantheonPlusSH0ES_STAT_SYS.cov"))
    n = int(cov_flat[0])
    cov = cov_flat[1:].reshape(n, n)
    assert n == len(dat)
    return dat, cov


def fit_cosmology(dat, cov):
    """Official Pantheon+SH0ES recipe: calibrators use Cepheid distances, Hubble-flow SNe
    (zHD > 0.01) use the model; M (SN absolute magnitude) is fitted jointly."""
    use = (dat["IS_CALIBRATOR"] == 1) | (dat["zHD"] > 0.01)
    idx = np.where(use)[0]
    C = cov[np.ix_(idx, idx)]
    cho = np.linalg.cholesky(C)
    d = dat.iloc[idx]
    mb = d["m_b_corr"].to_numpy()
    calib = d["IS_CALIBRATOR"].to_numpy() == 1
    ceph = d["CEPH_DIST"].to_numpy()
    zhd = d["zHD"].to_numpy(); zhel = d["zHEL"].to_numpy()

    def chi2(p, model):
        if model == "lcdm":
            M, H0, Om = p; w = -1.0
        else:
            M, H0, Om, w = p
        if not (0.01 < Om < 0.99) or not (40 < H0 < 110) or not (-3 < w < 0.5):
            return 1e12
        mu_th = np.where(calib, ceph, mu_theory(zhd, zhel, H0, Om, w))
        r = mb - M - mu_th
        y = np.linalg.solve(cho, r)
        return float(y @ y)

    out = {}
    for model, p0, steps in [("lcdm", [-19.25, 73.0, 0.33], [0.002, 0.05, 0.002]),
                             ("wcdm", [-19.25, 73.0, 0.33, -1.0], [0.002, 0.05, 0.002, 0.01])]:
        f = lambda p: chi2(p, model)
        res = optimize.minimize(f, p0, method="Nelder-Mead",
                                options=dict(xatol=1e-7, fatol=1e-6, maxiter=40000, maxfev=40000))
        res = optimize.minimize(f, res.x, method="Nelder-Mead",
                                options=dict(xatol=1e-8, fatol=1e-7, maxiter=40000, maxfev=40000))
        Hs = numeric_hessian(f, res.x, steps)
        pcov = 2.0 * np.linalg.inv(Hs)
        out[model] = dict(p=res.x.tolist(), cov=pcov.tolist(), chi2=res.fun, n=len(idx))
    return out


def cosmology_summary(fit):
    p = np.array(fit["lcdm"]["p"]); pc = np.array(fit["lcdm"]["cov"])
    M, H0, Om = p
    eH0, eOm = np.sqrt(np.diag(pc))[1:]
    # propagate by sampling
    s = rng.multivariate_normal(p, pc, 20000)
    H0s, Oms = s[:, 1], s[:, 2]
    HLs = H0s * np.sqrt(1 - Oms)                      # km/s/Mpc  (de Sitter rate of empty space)
    alphas = 3 * h_si(HLs)                            # 1/s volume growth rate
    rhoL = 3 * h_si(HLs) ** 2 / (8 * np.pi * G_SI)    # kg/m^3
    pw = np.array(fit["wcdm"]["p"]); pcw = np.array(fit["wcdm"]["cov"])
    HL = float(np.median(HLs)); alpha = float(np.median(alphas))
    summ = dict(
        M=M, H0=H0, H0_err=eH0, Om=Om, Om_err=eOm, OL=1 - Om,
        chi2=fit["lcdm"]["chi2"], n=fit["lcdm"]["n"],
        w=pw[3], w_err=float(np.sqrt(pcw[3, 3])), w_Om=pw[2], w_H0=pw[1],
        HL=HL, HL_err=float(np.std(HLs)),
        alpha_per_s=alpha, alpha_per_s_err=float(np.std(alphas)),
        alpha_per_yr=alpha * YR_S, alpha_per_yr_err=float(np.std(alphas) * YR_S),
        mm3_per_m3_per_yr=alpha * YR_S * 1e9,          # 1 m^3 grows by this many mm^3 per year
        volume_doubling_Gyr=math.log(2) / alpha / GYR_S,
        length_efold_Gyr=1.0 / h_si(HL) / GYR_S,
        rhoL_kg_m3=float(np.median(rhoL)), rhoL_err=float(np.std(rhoL)),
        rhoL_protons_per_m3=float(np.median(rhoL)) / M_P,
        rhoL_J_m3=float(np.median(rhoL)) * (C_KMS * 1e3) ** 2,
        # the rate of the whole universe today (matter coasting + vacuum)
        alpha_total_today_per_yr=3 * h_si(H0) * YR_S,
    )
    # what the volume rate would be with Planck's H0 (Hubble tension)
    H0p, Omp = 67.4, 0.315
    summ["alpha_per_yr_planck"] = 3 * h_si(H0p * math.sqrt(1 - Omp)) * YR_S
    return summ


# ============================================================================= PART B
def load_sparc(sdir):
    gals = []
    for fn in sorted(glob.glob(os.path.join(sdir, "*_rotmod.dat"))):
        name = os.path.basename(fn).replace("_rotmod.dat", "")
        dist = None
        with open(fn) as fh:
            for line in fh:
                if line.startswith("#") and "Distance" in line:
                    dist = float(line.split("=")[1].split()[0])
        a = np.loadtxt(fn, comments="#")
        if a.ndim == 1:
            a = a[None, :]
        r, vobs, ev, vgas, vdisk, vbul = (a[:, i] for i in range(6))
        vb2 = vgas * np.abs(vgas) + UPS_DISK * vdisk * np.abs(vdisk) + UPS_BUL * vbul * np.abs(vbul)
        ok = (r > 0) & (ev > 0) & (vobs > 0)
        gals.append(dict(name=name, D=dist, r=r[ok], v=vobs[ok], e=ev[ok], vb2=vb2[ok],
                         vgas=vgas[ok], vdisk=vdisk[ok], vbul=vbul[ok]))
    return gals


def load_btfr(fn):
    rows = []
    with open(fn) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) == 11 and not line.startswith(("Note", " Note")):
                try:
                    vals = [float(x) for x in parts[1:]]
                except ValueError:
                    continue
                rows.append([parts[0]] + vals)
    cols = ["name", "D", "eD", "fD", "logMb", "elogMb", "logVf", "elogVf", "Fg", "SBeff", "Reff"]
    return pd.DataFrame(rows, columns=cols)


def g_rar(gbar, a0):
    x = np.sqrt(np.clip(gbar, 1e-30, None) / a0)
    return gbar / (-np.expm1(-x))


def v_iso2(r, Vinf, rc):
    x = r / rc
    return Vinf ** 2 * (1 - np.arctan(x) / x)


def vmodel(kind, g, params=None, a0k=None, HLk=None):
    r, vb2 = g["r"], g["vb2"]
    if kind == "bar":
        v2 = vb2
    elif kind == "lit":                                   # literal: vacuum pushes outward
        v2 = vb2 - (HLk * r) ** 2
    elif kind == "free":                                  # generous: inward push, free strength
        v2 = vb2 + (params * r) ** 2
    elif kind == "univ":
        v2 = vb2 + (params * r) ** 2
    elif kind == "rar":
        gb = np.clip(vb2, 1e-6, None) / r
        v2 = g_rar(gb, a0k) * r
    elif kind == "iso":
        v2 = vb2 + v_iso2(r, *params)
    return np.sqrt(np.clip(v2, 0, None))


def chi2_of(g, vm):
    return float(np.sum(((g["v"] - vm) / g["e"]) ** 2))


def fit_free_h(g):
    """Best inward 'expansion push' strength h (km/s/kpc) for one galaxy: v^2 = vbar^2 + (h r)^2."""
    f = lambda h: chi2_of(g, vmodel("free", g, h))
    hs = np.logspace(-3, 3.3, 800)
    c2 = np.array([f(h) for h in hs])
    j = int(np.argmin(c2))
    lo, hi = hs[max(j - 1, 0)], hs[min(j + 1, len(hs) - 1)]
    res = optimize.minimize_scalar(f, bounds=(lo, hi), method="bounded",
                                   options=dict(xatol=1e-9))
    hb, cb = float(res.x), float(res.fun)
    c0 = f(0.0)
    if c0 <= cb:
        return 0.0, c0, (0.0, float(hs[np.argmax(c2 > c0 + 1)]) if np.any(c2 > c0 + 1) else np.nan)
    # 1-sigma interval from delta chi2 = 1
    target = cb + 1.0
    def bound(a, b):
        try:
            return optimize.brentq(lambda h: f(h) - target, a, b)
        except ValueError:
            return np.nan
    lo_b = bound(1e-6, hb) if f(1e-6) > target else 0.0
    hi_b = bound(hb, hs[-1] * 2) if f(hs[-1] * 2) > target else np.nan
    return hb, cb, (lo_b, hi_b)


def fit_iso(g):
    best = None
    for Vi0 in (20, 60, 120, 200, 300):
        for rc0 in (0.2, 1.0, 3.0, 10.0, 30.0):
            def resid(p):
                return (g["v"] - vmodel("iso", g, p)) / g["e"]
            try:
                res = optimize.least_squares(resid, [Vi0, rc0], bounds=([0.0, 0.01], [2000.0, 2000.0]))
            except Exception:
                continue
            c2 = float(np.sum(res.fun ** 2))
            if best is None or c2 < best[0]:
                best = (c2, res.x.tolist())
    return best


def rar_points(gals):
    """Points used for the radial acceleration relation: dV/V < 10% and positive baryonic gbar."""
    rows = []
    for gi, g in enumerate(gals):
        gb = g["vb2"] / g["r"]
        go = g["v"] ** 2 / g["r"]
        m = (g["e"] / g["v"] < 0.1) & (gb > 0)
        for k in np.where(m)[0]:
            rows.append((gi, g["r"][k], go[k] * CONV, gb[k] * CONV,
                         2 * g["e"][k] / (g["v"][k] * math.log(10))))
    a = np.array(rows)
    return dict(gi=a[:, 0].astype(int), r=a[:, 1], gobs=a[:, 2], gbar=a[:, 3], slog=a[:, 4])


def global_loglike_fit(P, kind, HL_si=None):
    """Max-likelihood fit in log g space with free intrinsic scatter.
    kind: 'rar' (a0), 'univ' (h), 'bar', 'lit'."""
    lgo = np.log10(P["gobs"])
    r_m = P["r"] * KPC_M

    def model(theta):
        if kind == "rar":
            return np.log10(g_rar(P["gbar"], 10 ** theta[0]))
        if kind == "univ":
            h = 10 ** theta[0]
            return np.log10(P["gbar"] + h ** 2 * r_m)
        if kind == "bar":
            return np.log10(P["gbar"])
        if kind == "lit":
            return np.log10(np.clip(P["gbar"] - HL_si ** 2 * r_m, 1e-30, None))

    def nll(theta):
        ls = theta[-1]
        s2 = P["slog"] ** 2 + (10 ** ls) ** 2
        res = lgo - model(theta)
        return 0.5 * float(np.sum(res ** 2 / s2 + np.log(2 * np.pi * s2)))

    if kind == "rar":
        x0 = [math.log10(1.2e-10), -1.0]
    elif kind == "univ":
        x0 = [math.log10(h_si(5000.0)), -0.7]
    else:
        x0 = [-0.5]
    res = optimize.minimize(nll, x0, method="Nelder-Mead",
                            options=dict(xatol=1e-8, fatol=1e-8, maxiter=20000))
    return res.x, res.fun


def binned_scatter(x, y, nb=12):
    edges = np.quantile(x, np.linspace(0, 1, nb + 1))
    sds, ns = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (x >= a) & (x <= b)
        if m.sum() > 10:
            sds.append(np.std(y[m])); ns.append(m.sum())
    sds, ns = np.array(sds), np.array(ns)
    return float(np.sqrt(np.sum(ns * sds ** 2) / np.sum(ns)))


def galaxy_analysis(gals, btfr, cosmo):
    HL = cosmo["HL"]                       # km/s/Mpc
    HLk = HL / 1000.0                      # km/s/kpc
    HL_si = h_si(HL)
    H0 = cosmo["H0"]

    # ---------- global fits on the RAR point set
    P = rar_points(gals)
    th_rar, nll_rar = global_loglike_fit(P, "rar")
    th_univ, nll_univ = global_loglike_fit(P, "univ")
    th_bar, nll_bar = global_loglike_fit(P, "bar")
    th_lit, nll_lit = global_loglike_fit(P, "lit", HL_si=HL_si)
    a0 = 10 ** th_rar[0]
    h_univ_si = 10 ** th_univ[0]
    H_univ = h_univ_si * MPC_M / 1e3       # km/s/Mpc
    npts = len(P["gobs"])

    # bootstrap a0 and H_univ over galaxies
    ug = np.unique(P["gi"])
    a0_bs, hu_bs = [], []
    for _ in range(150):
        pick = rng.choice(ug, size=len(ug), replace=True)
        idx = np.concatenate([np.where(P["gi"] == k)[0] for k in pick])
        Pb = {k: v[idx] for k, v in P.items()}
        a0_bs.append(10 ** global_loglike_fit(Pb, "rar")[0][0])
        hu_bs.append(10 ** global_loglike_fit(Pb, "univ")[0][0] * MPC_M / 1e3)
    a0k = a0 / CONV
    h_univ_k = H_univ / 1000.0

    # ---------- per-galaxy fits
    per = []
    for g in gals:
        n = len(g["r"])
        if n < 5:
            continue
        hb, c_free, (hlo, hhi) = fit_free_h(g)
        iso = fit_iso(g)
        c_bar = chi2_of(g, vmodel("bar", g))
        c_lit = chi2_of(g, vmodel("lit", g, HLk=HLk))
        c_rar = chi2_of(g, vmodel("rar", g, a0k=a0k))
        c_univ = chi2_of(g, vmodel("univ", g, h_univ_k))
        rmax = g["r"].max()
        # outer slope of missing acceleration
        gmiss = (g["v"] ** 2 - g["vb2"]) / g["r"]
        gobs = g["v"] ** 2 / g["r"]
        m = (g["r"] >= 0.5 * rmax) & (gmiss > 0) & (gmiss / gobs > 0.3)
        slope = np.nan
        if m.sum() >= 4:
            slope = float(np.polyfit(np.log(g["r"][m]), np.log(gmiss[m]), 1)[0])
        # dark-energy push vs missing pull at the last point
        k = int(np.argmax(g["r"]))
        a_lambda = HL_si ** 2 * g["r"][k] * KPC_M
        g_missing = gmiss[k] * CONV
        per.append(dict(
            name=g["name"], n=n, rmax=float(rmax), vlast=float(g["v"][k]),
            h_free=hb, h_lo=hlo, h_hi=hhi, H_free=hb * 1000.0,
            iso=iso[1] if iso else None,
            chi2=dict(bar=c_bar, lit=c_lit, free=c_free, univ=c_univ, rar=c_rar,
                      iso=iso[0] if iso else np.nan),
            slope=slope, a_lambda=a_lambda, g_missing=g_missing,
        ))

    # model comparison table
    kpar = dict(bar=0, lit=0, free=1, univ=0, rar=0, iso=2)
    table = {}
    for mk in ["bar", "lit", "free", "univ", "rar", "iso"]:
        c2 = np.array([p["chi2"][mk] for p in per])
        dof = np.array([max(p["n"] - kpar[mk], 1) for p in per])
        table[mk] = dict(median_chi2nu=float(np.median(c2 / dof)),
                         frac_chi2nu_lt2=float(np.mean(c2 / dof < 2)),
                         total_chi2=float(np.sum(c2)),
                         n_params_total=kpar[mk] * len(per) + (1 if mk in ("rar", "univ") else 0))
    # which model wins per galaxy (BIC) among free / rar / iso
    wins = dict(free=0, rar=0, iso=0)
    bic_free_vs_iso, bic_free_vs_rar = [], []
    for p in per:
        n = p["n"]
        b = {mk: p["chi2"][mk] + kpar[mk] * math.log(n) for mk in wins}
        wins[min(b, key=b.get)] += 1
        bic_free_vs_iso.append(b["free"] - b["iso"])
        bic_free_vs_rar.append(b["free"] - b["rar"])

    # fitted H_loc distribution (galaxies with a meaningful detection)
    Hf = np.array([p["H_free"] for p in per])
    Hpos = Hf[Hf > 0]
    # is it one universal value? compare per-galaxy chi2 at the universal value
    c_univ_tot = sum(p["chi2"]["univ"] for p in per)
    c_free_tot = sum(p["chi2"]["free"] for p in per)

    # correlations of H_loc with galaxy properties (BTFR table)
    bt = btfr.set_index("name")
    xs_m, xs_sb, ys = [], [], []
    for p in per:
        if p["name"] in bt.index and p["H_free"] > 0:
            xs_m.append(bt.loc[p["name"], "logMb"]); xs_sb.append(math.log10(bt.loc[p["name"], "SBeff"]))
            ys.append(math.log10(p["H_free"]))
    rho_m = stats.spearmanr(xs_m, ys)
    rho_sb = stats.spearmanr(xs_sb, ys)

    # slopes
    sl = np.array([p["slope"] for p in per if np.isfinite(p["slope"])])

    # vacuum push vs missing pull at last point
    ratio = np.array([p["a_lambda"] / p["g_missing"] for p in per if p["g_missing"] > 0])

    # residual pattern vs normalised radius (stacked)
    bins = np.linspace(0, 1, 11)
    respat = {mk: [] for mk in ("free", "rar", "iso")}
    gal_by_name = {g["name"]: g for g in gals}
    allx = {mk: [] for mk in respat}
    ally = {mk: [] for mk in respat}
    for p in per:
        g = gal_by_name[p["name"]]
        x = g["r"] / g["r"].max()
        vms = dict(free=vmodel("free", g, p["h_free"]), rar=vmodel("rar", g, a0k=a0k),
                   iso=vmodel("iso", g, p["iso"]) if p["iso"] else None)
        for mk in respat:
            if vms[mk] is None:
                continue
            allx[mk].append(x); ally[mk].append((g["v"] - vms[mk]) / g["v"])
    for mk in respat:
        x = np.concatenate(allx[mk]); y = np.concatenate(ally[mk])
        for a, b in zip(bins[:-1], bins[1:]):
            m = (x > a) & (x <= b)
            respat[mk].append(float(np.median(y[m])) if m.sum() else None)

    # discrepancy vs radius and vs baryonic acceleration
    lD = np.log10(P["gobs"] / P["gbar"])
    scat_r = binned_scatter(np.log10(P["r"]), lD)
    scat_g = binned_scatter(np.log10(P["gbar"]), lD)
    rho_r = stats.spearmanr(np.log10(P["r"]), lD)[0]
    rho_g = stats.spearmanr(np.log10(P["gbar"]), lD)[0]

    # RAR residual scatter (total, in dex)
    res_rar = np.log10(P["gobs"]) - np.log10(g_rar(P["gbar"], a0))
    res_univ = np.log10(P["gobs"]) - np.log10(P["gbar"] + h_univ_si ** 2 * P["r"] * KPC_M)
    res_bar = np.log10(P["gobs"]) - np.log10(P["gbar"])

    # BTFR
    bt2 = btfr.copy()
    slope_btfr, icpt_btfr = np.polyfit(bt2["logVf"], bt2["logMb"], 1)
    a0_impl = ((10 ** bt2["logVf"]) * 1e3) ** 4 / (G_SI * (10 ** bt2["logMb"]) * MSUN)
    # errors-in-both-variables fit for the slope: chi2 = sum (y - a x - b)^2 / (sy^2 + a^2 sx^2)
    xv_, yv_ = bt2["logVf"].to_numpy(), bt2["logMb"].to_numpy()
    sx_, sy_ = bt2["elogVf"].to_numpy(), bt2["elogMb"].to_numpy()
    def c2_btfr(p):
        a, b = p
        return float(np.sum((yv_ - a * xv_ - b) ** 2 / (sy_ ** 2 + a ** 2 * sx_ ** 2)))
    rb = optimize.minimize(c2_btfr, [4.0, 2.0], method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-9))
    bs_sl = []
    for _ in range(1000):
        i = rng.integers(0, len(xv_), len(xv_))
        def c2b(p, i=i):
            a, b = p
            return float(np.sum((yv_[i] - a * xv_[i] - b) ** 2 / (sy_[i] ** 2 + a ** 2 * sx_[i] ** 2)))
        bs_sl.append(optimize.minimize(c2b, rb.x, method="Nelder-Mead").x[0])
    class _O:  # small shim with the fields used below
        beta = rb.x
        sd_beta = [float(np.std(bs_sl)), 0.0]
    o = _O()

    gsum = dict(
        n_gal=len(gals), n_gal_fit=len(per), n_points_total=int(sum(len(g["r"]) for g in gals)),
        n_rar_points=npts, n_rar_gal=int(len(ug)),
        a0=a0, a0_err=float(np.std(a0_bs)), sig_int_rar=10 ** th_rar[1],
        H_univ=H_univ, H_univ_err=float(np.std(hu_bs)), sig_int_univ=10 ** th_univ[1],
        sig_int_bar=10 ** th_bar[0], sig_int_lit=10 ** th_lit[0],
        nll=dict(rar=nll_rar, univ=nll_univ, bar=nll_bar, lit=nll_lit),
        dBIC_univ_minus_rar=2 * (nll_univ - nll_rar),
        rms_dex=dict(rar=float(np.std(res_rar)), univ=float(np.std(res_univ)), bar=float(np.std(res_bar))),
        mean_dex=dict(rar=float(np.mean(res_rar)), univ=float(np.mean(res_univ)), bar=float(np.mean(res_bar))),
        table=table, wins=wins,
        frac_free_beats_iso=float(np.mean(np.array(bic_free_vs_iso) < 0)),
        frac_free_beats_rar=float(np.mean(np.array(bic_free_vs_rar) < 0)),
        H_free_median=float(np.median(Hpos)), H_free_p16=float(np.percentile(Hpos, 16)),
        H_free_p84=float(np.percentile(Hpos, 84)), H_free_min=float(Hpos.min()), H_free_max=float(Hpos.max()),
        H_free_n_pos=int(len(Hpos)), H_free_n_zero=int(np.sum(Hf == 0)),
        H_free_dex_scatter=float(np.std(np.log10(Hpos))),
        chi2_univ_total=c_univ_tot, chi2_free_total=c_free_tot,
        spearman_H_vs_logMb=[float(rho_m[0]), float(rho_m[1])],
        spearman_H_vs_logSB=[float(rho_sb[0]), float(rho_sb[1])],
        n_spearman=len(ys),
        slope_n=int(len(sl)), slope_median=float(np.median(sl)),
        slope_p16=float(np.percentile(sl, 16)), slope_p84=float(np.percentile(sl, 84)),
        slope_frac_gt0=float(np.mean(sl > 0)), slope_frac_gt05=float(np.mean(sl > 0.5)),
        lambda_ratio_median=float(np.median(ratio)), lambda_ratio_p16=float(np.percentile(ratio, 16)),
        lambda_ratio_p84=float(np.percentile(ratio, 84)), lambda_ratio_n=int(len(ratio)),
        g_missing_last_median=float(np.median([p["g_missing"] for p in per if p["g_missing"] > 0])),
        r_last_median_kpc=float(np.median([p["rmax"] for p in per])),
        a_lambda_last_median=float(np.median([p["a_lambda"] for p in per])),
        disc_scatter_by_r=scat_r, disc_scatter_by_g=scat_g, rho_disc_r=float(rho_r), rho_disc_g=float(rho_g),
        btfr_n=len(bt2), btfr_slope_ols=float(slope_btfr), btfr_slope_odr=float(o.beta[0]),
        btfr_slope_odr_err=float(o.sd_beta[0]),
        btfr_a0_median=float(np.median(a0_impl)), btfr_a0_scatter_dex=float(np.std(np.log10(a0_impl))),
        btfr_scatter_dex=float(np.std(bt2["logMb"] - (o.beta[0] * bt2["logVf"] + o.beta[1]))),
        respat=respat, respat_bins=bins.tolist(),
        HL=HL, H0=H0,
    )
    return gsum, per, P, (a0, a0k, H_univ, h_univ_k, HLk)


# ============================================================================= PART C
def unique_sne(dat, M_fit):
    """One row per supernova (some were observed by several surveys); distance modulus uses the
    absolute magnitude M from our own cosmology fit so everything is internally consistent."""
    d = dat.copy()
    d["MU"] = d["m_b_corr"] - M_fit
    w = 1 / d["m_b_corr_err_RAW"] ** 2
    d["w"] = w
    agg = []
    for cid, grp in d.groupby("CID", sort=False):
        ww = grp["w"].to_numpy()
        row = grp.iloc[0].copy()
        row["MU"] = float(np.sum(ww * grp["MU"]) / np.sum(ww))
        row["err_raw"] = float(1 / np.sqrt(np.sum(ww)))
        agg.append(row)
    return pd.DataFrame(agg).reset_index(drop=True)


def fit_noise(lnH, sig_raw_mag, cz, sig_int_mag=None):
    """Max-likelihood split of the scatter in ln(H_i) into supernova noise
    [ (ln10/5)^2 (raw^2 + sig_int^2) ] and random peculiar velocities [ (sig_v / cz)^2 ]."""
    def nll(th):
        off = th[0]
        si = sig_int_mag if sig_int_mag is not None else 10 ** th[2]
        sv = 10 ** th[1]
        s2 = LN10_5 ** 2 * (sig_raw_mag ** 2 + si ** 2) + (sv / cz) ** 2
        return 0.5 * float(np.sum((lnH - off) ** 2 / s2 + np.log(s2)))
    x0 = [0.0, math.log10(250.0)] + ([] if sig_int_mag is not None else [math.log10(0.1)])
    res = optimize.minimize(nll, x0, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-7, maxiter=20000))
    Hs = numeric_hessian(nll, res.x, [1e-3] * len(res.x))
    cov = np.linalg.inv(Hs)
    sv = 10 ** res.x[1]
    sv_err = sv * math.log(10) * math.sqrt(cov[1, 1])
    si = sig_int_mag if sig_int_mag is not None else 10 ** res.x[2]
    return dict(offset=float(res.x[0]), sig_v=float(sv), sig_v_err=float(sv_err), sig_int=float(si))


def radec_to_gal(ra, dec):
    R = np.array([[-0.0548755604, -0.8734370902, -0.4838350155],
                  [0.4941094279, -0.4448296300, 0.7469822445],
                  [-0.8676661490, -0.1980763734, 0.4559837762]])
    a, d = np.radians(ra), np.radians(dec)
    x = np.vstack([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)])
    xg = R @ x
    l = np.degrees(np.arctan2(xg[1], xg[0])) % 360
    b = np.degrees(np.arcsin(np.clip(xg[2], -1, 1)))
    return l, b, xg.T


def mollweide(l_deg, b_deg):
    lon = np.radians(((l_deg + 180) % 360) - 180)
    lat = np.radians(b_deg)
    th = lat.copy()
    for _ in range(60):
        den = 2 + 2 * np.cos(2 * th)
        den = np.where(np.abs(den) < 1e-9, 1e-9, den)
        th = th - (2 * th + np.sin(2 * th) - np.pi * np.sin(lat)) / den
    x = 2 * np.sqrt(2) / np.pi * lon * np.cos(th)
    y = np.sqrt(2) * np.sin(th)
    return -x, y    # longitude increasing to the left (sky convention)


def sigma_R(R_mpc, h=0.70, Om=0.31, Ob=0.049, ns=0.965, s8=0.81):
    """Linear rms density contrast in a top-hat sphere (BBKS transfer function)."""
    Gamma = Om * h * math.exp(-Ob * (1 + math.sqrt(2 * h) / Om))
    k = np.logspace(-4, 2, 6000)
    q = k / Gamma
    T = np.log(1 + 2.34 * q) / (2.34 * q) * (1 + 3.89 * q + (16.1 * q) ** 2 + (5.46 * q) ** 3 + (6.71 * q) ** 4) ** -0.25
    P = k ** ns * T ** 2

    def s2(Rh):
        x = k * Rh
        W = 3 * (np.sin(x) - x * np.cos(x)) / x ** 3
        return integrate.trapezoid(k ** 2 * P * W ** 2, k) / (2 * np.pi ** 2)

    norm = s8 ** 2 / s2(8.0)
    return np.array([math.sqrt(norm * s2(R * h)) for R in np.atleast_1d(R_mpc)])


def hubble_flow_analysis(dat, cosmo):
    H0, Om = cosmo["H0"], cosmo["Om"]
    u = unique_sne(dat, cosmo["M"])
    u = u[(u["zCMB"] > 0.0035) & (u["zHD"] < 0.8)].copy()
    zg, dc, E = comoving_grid(Om)
    dl_obs = 10 ** ((u["MU"] - 25) / 5)                      # Mpc
    chi_obs = dl_obs / (1 + u["zHEL"])
    # cosmological redshift implied by the measured distance
    zcos = np.interp(chi_obs * H0 / C_KMS, dc, zg)
    v_obs = C_KMS * ((1 + u["zCMB"]) / (1 + zcos) - 1)
    u["d"] = dl_obs / (1 + zcos)                             # ~ comoving distance today (Mpc)
    u["zcos"] = zcos; u["v_obs"] = v_obs
    # per-supernova "Hubble constant": H such that model distance at the observed z equals measured
    mu_raw = mu_theory(u["zCMB"].to_numpy(), u["zHEL"].to_numpy(), H0, Om)
    mu_cor = mu_theory(u["zHD"].to_numpy(), u["zHEL"].to_numpy(), H0, Om)
    u["H_raw"] = H0 * 10 ** ((mu_raw - u["MU"]) / 5)
    u["H_cor"] = H0 * 10 ** ((mu_cor - u["MU"]) / 5)

    # --- noise model: SN scatter (raw error + fitted intrinsic) and random peculiar velocities
    fitset = u[u["zHD"] < 0.15]
    cz_f = C_KMS * fitset["zcos"].to_numpy()
    nz_cor = fit_noise(np.log(fitset["H_cor"] / H0).to_numpy(), fitset["err_raw"].to_numpy(), cz_f)
    nz_raw = fit_noise(np.log(fitset["H_raw"] / H0).to_numpy(), fitset["err_raw"].to_numpy(), cz_f,
                       sig_int_mag=nz_cor["sig_int"])
    sig_int = nz_cor["sig_int"]
    u["sig_mu_sn"] = np.sqrt(u["err_raw"] ** 2 + sig_int ** 2)
    u["sig_lnH_sn"] = LN10_5 * u["sig_mu_sn"]
    Ez = np.interp(zcos, zg, E)
    dz = chi_obs * LN10_5 * u["sig_mu_sn"] * H0 * Ez / C_KMS
    u["sv"] = C_KMS * (1 + u["zCMB"]) / (1 + zcos) ** 2 * dz   # velocity error from distance error

    loc = u[u["zHD"] < 0.15].copy()

    # binned scatter of per-SN Hubble rates
    edges = [15, 30, 50, 80, 130, 200, 320, 650]
    bins = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (loc["d"] >= a) & (loc["d"] < b)
        if m.sum() < 8:
            continue
        sub = loc[m]
        lr = np.log(sub["H_raw"] / H0); lc = np.log(sub["H_cor"] / H0)
        s_sn = np.sqrt(np.mean(sub["sig_lnH_sn"] ** 2))
        sd_r, sd_c = np.std(lr), np.std(lc)
        cz = np.median(C_KMS * sub["zcos"])
        sv_r = cz * math.sqrt(max(sd_r ** 2 - s_sn ** 2, 0))     # implied peculiar-velocity dispersion
        sv_c = cz * math.sqrt(max(sd_c ** 2 - s_sn ** 2, 0))
        bins.append(dict(d_lo=a, d_hi=b, n=int(m.sum()), d_med=float(np.median(sub["d"])),
                         sd_raw=float(sd_r), sd_cor=float(sd_c), sd_sn=float(s_sn),
                         p16_raw=float(np.percentile(sub["H_raw"], 16)), p84_raw=float(np.percentile(sub["H_raw"], 84)),
                         p16_cor=float(np.percentile(sub["H_cor"], 16)), p84_cor=float(np.percentile(sub["H_cor"], 84)),
                         min_raw=float(sub["H_raw"].min()), max_raw=float(sub["H_raw"].max()),
                         sigv_raw=float(sv_r), sigv_cor=float(sv_c)))

    # observed peculiar velocity vs gravity-predicted (2M++ density field) peculiar velocity
    reg = loc[(loc["zCMB"] > 0.01) & (loc["zCMB"] < 0.06)].copy()
    sig_star = nz_cor["sig_v"]

    def wls(x, y, s):
        w = 1 / s ** 2
        X = np.vstack([np.ones_like(x), x]).T
        A = X.T @ (X * w[:, None])
        bvec = X.T @ (w * y)
        beta = np.linalg.solve(A, bvec)
        return beta

    s_tot = np.sqrt(reg["sv"] ** 2 + sig_star ** 2).to_numpy()
    xv, yv = reg["VPEC"].to_numpy(), reg["v_obs"].to_numpy()
    beta = wls(xv, yv, s_tot)
    bs = []
    for _ in range(2000):
        i = rng.integers(0, len(reg), len(reg))
        bs.append(wls(xv[i], yv[i], s_tot[i])[1])
    corr = stats.pearsonr(xv, yv)

    # bulk flow (dipole) fit, galactic coords
    ra = np.where(reg["HOST_RA"] > -900, reg["HOST_RA"], reg["RA"])
    de = np.where(reg["HOST_DEC"] > -900, reg["HOST_DEC"], reg["DEC"])
    l, b, nvec = radec_to_gal(ra, de)

    def dipole(v, s):
        w = 1 / s ** 2
        X = np.hstack([np.ones((len(v), 1)), nvec])
        A = X.T @ (X * w[:, None])
        beta_ = np.linalg.solve(A, X.T @ (w * v))
        cov = np.linalg.inv(A)
        B = beta_[1:]
        amp = float(np.linalg.norm(B))
        lb = (float(np.degrees(np.arctan2(B[1], B[0])) % 360), float(np.degrees(np.arcsin(B[2] / amp))))
        # amplitude uncertainty via sampling
        smp = rng.multivariate_normal(beta_, cov, 4000)[:, 1:]
        return dict(amp=amp, amp_err=float(np.std(np.linalg.norm(smp, axis=1))), l=lb[0], b=lb[1],
                    monopole=float(beta_[0]))
    s_r = np.sqrt(reg["sv"] ** 2 + nz_raw["sig_v"] ** 2).to_numpy()
    s_c = np.sqrt(reg["sv"] ** 2 + sig_star ** 2).to_numpy()
    dip_raw = dipole(yv, s_r)
    dip_cor = dipole(yv - xv, s_c)

    # shells: local Hubble rate in distance shells (raw redshifts vs gravity-corrected)
    sh_edges = [15, 40, 80, 140, 220, 330, 480, 650]
    shells = []
    for a, b_ in zip(sh_edges[:-1], sh_edges[1:]):
        m = (loc["d"] >= a) & (loc["d"] < b_)
        if m.sum() < 8:
            continue
        sub = loc[m]
        cz = C_KMS * sub["zcos"]
        out = dict(d_lo=a, d_hi=b_, n=int(m.sum()), d_mid=float(np.median(sub["d"])))
        for key, col, sv_ in (("raw", "H_raw", nz_raw["sig_v"]), ("cor", "H_cor", nz_cor["sig_v"])):
            s2 = sub["sig_lnH_sn"] ** 2 + (sv_ / cz) ** 2
            w = 1 / s2
            mu = float(np.sum(w * np.log(sub[col] / H0)) / np.sum(w))
            err = float(1 / math.sqrt(np.sum(w)))
            out[key] = H0 * math.exp(mu)
            out[key + "_err"] = H0 * math.exp(mu) * err
        shells.append(out)
    # LCDM cosmic-variance band for the mean expansion inside radius R (linear theory)
    f_growth = Om ** 0.55
    Rs = np.array([20, 30, 45, 70, 100, 150, 220, 330, 500, 650.0])
    sigH = f_growth / 3 * sigma_R(Rs, h=H0 / 100, Om=Om)
    band = dict(R=Rs.tolist(), frac=sigH.tolist())

    # sky map data (z < 0.04)
    sky = loc[(loc["zCMB"] < 0.04)].copy()
    ra_s = np.where(sky["HOST_RA"] > -900, sky["HOST_RA"], sky["RA"])
    de_s = np.where(sky["HOST_DEC"] > -900, sky["HOST_DEC"], sky["DEC"])
    ls_, bs_, _ = radec_to_gal(ra_s, de_s)
    mx, my = mollweide(ls_, bs_)
    sky_pts = [[round(float(a), 4), round(float(b), 4), round(float(c), 4), round(float(d), 4), round(float(e), 1)]
               for a, b, c, d, e in zip(mx, my, sky["H_raw"] / H0 - 1, sky["H_cor"] / H0 - 1, sky["d"])]

    # inside vs outside ~300 Mpc (local-void test)
    def wmean_H(sub, col, sv_):
        cz = C_KMS * sub["zcos"]
        w = 1 / (sub["sig_lnH_sn"] ** 2 + (sv_ / cz) ** 2)
        mu = float(np.sum(w * np.log(sub[col] / H0)) / np.sum(w))
        return H0 * math.exp(mu), H0 * math.exp(mu) / math.sqrt(np.sum(w))
    inner = u[(u["zHD"] > 0.023) & (u["d"] < 300)]
    outer = u[(u["d"] >= 300) & (u["zHD"] < 0.15)]
    Hin, eHin = wmean_H(inner, "H_cor", nz_cor["sig_v"])
    Hout, eHout = wmean_H(outer, "H_cor", nz_cor["sig_v"])
    Hin_r, eHin_r = wmean_H(inner, "H_raw", nz_raw["sig_v"])
    Hout_r, eHout_r = wmean_H(outer, "H_raw", nz_raw["sig_v"])
    # expected spread of individual "Hubble constants" vs distance from the fitted noise model
    dgrid = np.geomspace(15, 650, 60)
    zgrid = np.interp(dgrid * H0 / C_KMS, dc, zg)
    sn_typ = float(np.median(loc["sig_lnH_sn"]))
    exp_raw = np.sqrt(sn_typ ** 2 + (nz_raw["sig_v"] / (C_KMS * zgrid)) ** 2)
    exp_sn = np.full_like(dgrid, sn_typ)
    expected = dict(d=dgrid.round(2).tolist(), raw=exp_raw.tolist(), sn=exp_sn.tolist())

    pts = loc[["d", "H_raw", "H_cor", "sig_lnH_sn", "zcos"]].to_numpy()
    sn_pts = [[round(float(a), 2), round(float(b), 2), round(float(c), 2), round(float(dd), 4)]
              for a, b, c, dd, _ in pts]
    reg_pts = [[round(float(a), 0), round(float(b), 0), round(float(c), 0)] for a, b, c in
               zip(reg["VPEC"], reg["v_obs"], reg["sv"])]

    return dict(
        n_unique=int(len(u)), n_local=int(len(loc)), n_reg=int(len(reg)),
        bins=bins, reg_slope=float(beta[1]), reg_slope_err=float(np.std(bs)), reg_icpt=float(beta[0]),
        reg_pearson=[float(corr[0]), float(corr[1])],
        noise_raw=nz_raw, noise_cor=nz_cor, sig_star=sig_star,
        dipole_raw=dip_raw, dipole_cor=dip_cor, shells=shells, band=band, expected=expected,
        H_in300=Hin, H_in300_err=eHin, H_out300=Hout, H_out300_err=eHout,
        H_in300_raw=Hin_r, H_in300_raw_err=eHin_r, H_out300_raw=Hout_r, H_out300_raw_err=eHout_r,
        f_growth=f_growth, sn_pts=sn_pts, reg_pts=reg_pts, sky_pts=sky_pts,
    )


# ============================================================================= order-of-magnitude
def scale_numbers(cosmo, gsum):
    HL_si = h_si(cosmo["HL"]); H0_si = h_si(cosmo["H0"])
    r30 = 30 * KPC_M
    a_lambda_30 = HL_si ** 2 * r30
    cH0 = C_KMS * 1e3 * H0_si
    # stellar volume filling factor, solar neighbourhood: ~0.1 stars/pc^3 with solar size (generous)
    v_sun = 4 / 3 * math.pi * (6.957e8) ** 3
    pc3 = (KPC_M / 1000) ** 3
    fill_local = 0.1 * v_sun / pc3
    fill_core = 1e6 * v_sun / pc3          # ~1e6 stars/pc^3 near the Galactic Centre
    # local dark matter density ~0.4 GeV/cm^3
    rho_dm_local = 0.4 * GEV_CM3
    # where vacuum repulsion balances the gravity of the Local Group (M ~ 2e12 Msun)
    M_lg = 2e12 * MSUN
    R_zg = (3 * M_lg / (8 * math.pi * cosmo["rhoL_kg_m3"])) ** (1 / 3) / MPC_M
    # Sun's orbit
    a_sun = (233e3) ** 2 / (8.2 * KPC_M)
    ladder = [
        dict(label="Dark-energy push at 30 kpc (your rule, measured strength)", a=a_lambda_30, kind="push"),
        dict(label="Missing pull at the edge of a typical galaxy (SPARC median)", a=gsum["g_missing_last_median"], kind="need"),
        dict(label="Acceleration scale where the 'missing mass' appears (fit, a₀)", a=gsum["a0"], kind="scale"),
        dict(label="c·H₀ / 2π  (speed of light × Hubble rate)", a=cH0 / (2 * math.pi), kind="scale"),
        dict(label="Sun's acceleration around the Milky Way", a=a_sun, kind="ref"),
        dict(label="Earth's acceleration around the Sun", a=5.93e-3, kind="ref"),
        dict(label="Gravity at Earth's surface", a=9.81, kind="ref"),
    ]
    return dict(
        a_lambda_30=a_lambda_30, cH0=cH0, cH0_2pi=cH0 / (2 * math.pi),
        horizon_over_galaxy=C_KMS / (cosmo["H0"] * 0.030),
        hubble_velocity_30kpc=cosmo["H0"] * 0.030, hubble_velocity_30kpc_HL=cosmo["HL"] * 0.030,
        fill_local=fill_local, fill_core=fill_core,
        rho_dm_local=rho_dm_local, rho_dm_local_protons=rho_dm_local / M_P,
        rho_ratio_dm_to_lambda=rho_dm_local / cosmo["rhoL_kg_m3"],
        R_zero_gravity_LG_Mpc=R_zg, ladder=ladder,
    )


# ============================================================================= main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sparc", required=True)
    ap.add_argument("--btfr", required=True)
    ap.add_argument("--pantheon", required=True)
    ap.add_argument("--out", default="results.json")
    ap.add_argument("--explorer", default="NGC3198,NGC2403,DDO154,UGC02885,F568-3,NGC7814,IC2574")
    args = ap.parse_args()

    print("A) Cosmology fit (Pantheon+SH0ES) ...")
    dat, cov = load_pantheon(args.pantheon)
    fit = fit_cosmology(dat, cov)
    cosmo = cosmology_summary(fit)
    for k in ("H0", "H0_err", "Om", "Om_err", "w", "w_err", "HL", "alpha_per_yr", "mm3_per_m3_per_yr",
              "volume_doubling_Gyr", "rhoL_kg_m3", "rhoL_protons_per_m3", "chi2", "n"):
        print(f"   {k:22s} {cosmo[k]:.6g}")

    print("B) Galaxies (SPARC) ...")
    gals = load_sparc(args.sparc)
    btfr = load_btfr(args.btfr)
    gsum, per, P, (a0, a0k, H_univ, h_univ_k, HLk) = galaxy_analysis(gals, btfr, cosmo)
    for k, v in gsum.items():
        if k not in ("table", "respat", "respat_bins", "nll"):
            print(f"   {k:28s} {v}")
    print("   model table:")
    for mk, row in gsum["table"].items():
        print(f"     {mk:5s} {row}")

    print("C) Hubble flow (Pantheon+) ...")
    hf = hubble_flow_analysis(dat, cosmo)
    for k, v in hf.items():
        if k not in ("sn_pts", "reg_pts", "sky_pts"):
            print(f"   {k}: {v}")

    scales = scale_numbers(cosmo, gsum)
    print("Scales:", {k: v for k, v in scales.items() if k != "ladder"})

    # ------------------------------------------------------------- export for the report
    gal_by_name = {g["name"]: g for g in gals}
    per_by_name = {p["name"]: p for p in per}
    explorer = []
    for nm in args.explorer.split(","):
        g = gal_by_name[nm]; p = per_by_name[nm]
        explorer.append(dict(
            name=nm, D=g["D"], r=[round(float(x), 3) for x in g["r"]], v=[round(float(x), 2) for x in g["v"]],
            e=[round(float(x), 2) for x in g["e"]], vb2=[round(float(x), 2) for x in g["vb2"]],
            h_free=p["h_free"], iso=p["iso"],
            chi2nu={k: (v / max(p["n"] - dict(free=1, iso=2).get(k, 0), 1)) for k, v in p["chi2"].items()},
        ))
    rar_export = [[round(float(math.log10(a)), 4), round(float(math.log10(b)), 4), int(c), round(float(d), 3)]
                  for a, b, c, d in zip(P["gbar"], P["gobs"], P["gi"], P["r"])]
    per_export = [dict(name=p["name"], n=p["n"], H=p["H_free"], slope=p["slope"],
                       ratio=(p["a_lambda"] / p["g_missing"]) if p["g_missing"] > 0 else None,
                       c2nu_free=p["chi2"]["free"] / max(p["n"] - 1, 1),
                       c2nu_rar=p["chi2"]["rar"] / max(p["n"], 1),
                       c2nu_iso=p["chi2"]["iso"] / max(p["n"] - 2, 1),
                       c2nu_bar=p["chi2"]["bar"] / max(p["n"], 1),
                       rmax=p["rmax"]) for p in per]
    btfr_export = [[round(float(a), 3), round(float(b), 3), nm] for a, b, nm in
                   zip(btfr["logVf"], btfr["logMb"], btfr["name"])]
    out = dict(cosmo=cosmo, galaxies=gsum, explorer=explorer, rar=rar_export, per=per_export,
               gal_names=[g["name"] for g in gals], btfr=btfr_export, hubble=hf, scales=scales)

    def clean(o):
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, (np.floating, float)):
            return None if not np.isfinite(o) else float(o)
        if isinstance(o, np.integer):
            return int(o)
        return o
    with open(args.out, "w") as fh:
        json.dump(clean(out), fh, separators=(",", ":"))
    print("wrote", args.out, os.path.getsize(args.out), "bytes")


if __name__ == "__main__":
    main()
