"""
Compare the quasi-static mu (Rebouças 2026 eq. 17, use_qsa=True) against the
scale-dependent mu of Cataneo+2024 eq. (31a) (use_qsa=False), for the
alpha_K proportional to Omega_DE parametrization (alpha_K_parametrization=1).

Usage: python docs/qsa_comparison.py   ->  writes docs/qsa_comparison.png
"""

import os

import matplotlib.pyplot as plt
import numpy as np

import camb

# Background / DE settings shared by every run
H0, ombh2, omch2 = 67.5, 0.022, 0.12
w0, wa, cs2 = -0.838, -0.6, 1.0
alpha_K_0_values = [1.0, 3.0, 10.0]  # mu_p > 0 needs cs2*alpha_K_0 >~ 1 here
lmax = 2000
z_mu = 0.5  # redshift at which mu(k) is shown

# Ordered values -> one-hue ramp, light to dark (blue 300/450/650)
colors = ["#6da7ec", "#2a78d6", "#104281"]
ink, muted, grid = "#1f1f1e", "#6b6a64", "#e4e3dd"


def get_pars(**kw):
    pars = camb.CAMBparams()
    pars.DarkEnergy = camb.dark_energy.DarkEnergyPPF()
    pars.DarkEnergy.set_params(w=w0, wa=wa, cs2_0=cs2)
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, alpha_K_parametrization=1, **kw)
    pars.InitPower.set_params(As=2e-9, ns=0.965)
    pars.set_for_lmax(lmax, lens_potential_accuracy=1)
    pars.set_matter_power(redshifts=[0.0], kmax=2.0)
    return pars


def spectra(results):
    kh, _, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints=300)
    cls = results.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)
    return kh, pk[0], cls["total"][2 : lmax + 1, 0], cls["lens_potential"][2 : lmax + 1, 0]


def mu_of_k(results, k_mpc, z):
    """mu(k, z) rebuilt from the Fortran tables, as in equations.f90 (MG_interp_tables)."""
    P = results.Params
    log_a, alpha_B = np.array(P.log_a), np.array(P.alpha_B)
    c_sN2, mu_p = np.array(P.c_sN2), np.array(P.mu_p)
    r1 = np.interp(-np.log1p(z), log_a, mu_p / c_sN2)
    mu_inf = np.interp(-np.log1p(z), log_a, 1 + alpha_B**2 / (2 * c_sN2))
    calH = results.hubble_parameter(z) / (1 + z) / 299792.458  # aH in 1/Mpc
    x = (k_mpc / calH) ** 2
    return (r1 + x * mu_inf) / (r1 + x), mu_inf


# GR reference (same dark energy background, no MG)
kh, pk_gr, cl_gr, clpp_gr = spectra(camb.get_results(get_pars()))
ell = np.arange(2, lmax + 1)

fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
ax_mu, ax_pk, ax_tt, ax_pp = axes.flat

for aK0, color in zip(alpha_K_0_values, colors):
    for use_qsa, ls in [(True, "--"), (False, "-")]:
        res = camb.get_results(get_pars(use_cs2=True, use_qsa=use_qsa, alpha_K_0=aK0))
        _, pk, cl, clpp = spectra(res)
        ax_pk.plot(kh, pk / pk_gr, color=color, ls=ls, lw=2)
        ax_tt.plot(ell, cl / cl_gr, color=color, ls=ls, lw=2)
        ax_pp.plot(ell, clpp / clpp_gr, color=color, ls=ls, lw=2)
        if not use_qsa:
            mu, mu_inf = mu_of_k(res, kh * H0 / 100, z_mu)
            ax_mu.plot(kh, mu, color=color, ls="-", lw=2)
            ax_mu.plot(kh, np.full_like(kh, mu_inf), color=color, ls="--", lw=2)

# Axes styling
panels = [
    (ax_mu, r"$\mu(k,\,z=%g)$" % z_mu, r"$k\ [h\,\mathrm{Mpc}^{-1}]$", "log"),
    (ax_pk, r"$P(k,\,z=0)\,/\,P_\mathrm{GR}$", r"$k\ [h\,\mathrm{Mpc}^{-1}]$", "log"),
    (ax_tt, r"$C_\ell^{TT}\,/\,C_{\ell,\mathrm{GR}}^{TT}$", r"$\ell$", "log"),
    (ax_pp, r"$C_\ell^{\phi\phi}\,/\,C_{\ell,\mathrm{GR}}^{\phi\phi}$", r"$\ell$", "log"),
]
for ax, ylabel, xlabel, xscale in panels:
    ax.set_xscale(xscale)
    ax.axhline(1, color=muted, lw=1, ls=":", zorder=0)
    ax.set_xlabel(xlabel, color=ink)
    ax.set_ylabel(ylabel, color=ink)
    ax.grid(True, color=grid, lw=0.8)
    ax.tick_params(colors=muted, labelcolor=ink)
    for spine in ax.spines.values():
        spine.set_color(grid)
ax_tt.set_xlim(2, lmax)
ax_pp.set_xlim(2, lmax)
# Legend: color = alpha_K_0, line style = mu prescription
from matplotlib.lines import Line2D

handles = [Line2D([], [], color=c, lw=2, label=r"$\alpha_{K,0} = %g$" % a) for a, c in zip(alpha_K_0_values, colors)]
handles += [
    Line2D([], [], color=ink, lw=2, ls="--", label="QSA (eq. 17)"),
    Line2D([], [], color=ink, lw=2, ls="-", label="scale-dependent (Cataneo+24 eq. 31a)"),
]
fig.legend(handles=handles, loc="outside lower center", ncol=5, frameon=False)
fig.suptitle(
    r"QSA vs scale-dependent $\mu$:  $\alpha_K \propto \Omega_\mathrm{DE}$,  PPF DE with $w_0=%g$, $w_a=%g$, $c_s^2=%g$" % (w0, wa, cs2),
    color=ink,
)

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "qsa_comparison.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print("saved", out)
