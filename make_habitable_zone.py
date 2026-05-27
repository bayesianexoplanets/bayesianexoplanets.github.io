"""Habitable-zone diagram for the TESS catalog.

Generates ``habitable_zone.png`` in this directory. Two side-by-side panels show
the known TOIs (left) and the new candidates (right) in the
(instellation flux, T_eff) plane, with the Kopparapu (2014) conservative and
optimistic HZ regions shaded. Points are colored by ``1 - NST p-value`` (NST
column ``log10(p value)`` in ``tois.csv`` / ``tois_new.csv``). Candidates that
fail the systematics tests (``passed_all_tests == False``) get a dashed red
ring around them; the three newly-identified HZ candidates are highlighted
with gold stars.

Style follows the Steve Bryson DR25 insolation figure:
https://github.com/stevepur/DR25-occurrence-public/blob/main/insolation/insolation_figures.ipynb
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "habitable_zone.png")

# Kopparapu (2014) HZ-flux polynomial.
# i: 0 = Recent Venus, 1 = Runaway Greenhouse, 2 = Maximum Greenhouse, 3 = Early Mars.
SEFFSUN = [1.776, 1.107, 0.356, 0.320]
A = [2.136e-4, 1.332e-4, 6.171e-5, 5.547e-5]
B = [2.533e-8, 1.580e-8, 1.698e-9, 1.526e-9]
C = [-1.332e-11, -8.308e-12, -3.198e-12, -2.874e-12]
D = [-3.097e-15, -1.931e-15, -5.575e-16, -5.011e-16]


def hz_edge(teff, i):
    Ts = teff - 5780.0
    return SEFFSUN[i] + A[i] * Ts + B[i] * Ts**2 + C[i] * Ts**3 + D[i] * Ts**4


def insolation(rstar, logg, teff, period):
    """Bolometric flux in Earth units (S/S_earth)."""
    mstar = 10.0**logg * rstar**2.0 / 10.0**4.437
    semia = mstar ** (1.0 / 3.0) * (period / 365.25) ** (2.0 / 3.0)
    lum = rstar**2.0 * (teff / 5778.0) ** 4.0
    return lum / semia**2.0


def prep(df):
    need = ["Radius", "logg", "Teff", "Period"]
    d = df.dropna(subset=need).copy()
    d["S"] = insolation(d.Radius, d.logg, d.Teff, d.Period)
    p = 10.0 ** d["log10(p value)"].astype(float)
    d["conf"] = (1.0 - p).clip(0.0, 1.0)
    return d


def in_hz(df, mode):
    j_in, j_out = (0, 3) if mode == "optimistic" else (1, 2)
    return (df.S <= hz_edge(df.Teff, j_in)) & (df.S >= hz_edge(df.Teff, j_out))


known = prep(pd.read_csv(os.path.join(HERE, "tois.csv")))
new = prep(pd.read_csv(os.path.join(HERE, "tois_new.csv")))

# Newly-identified HZ candidates, surfaced from tois_new.csv.
HL_TICS = sorted(new.loc[in_hz(new, "optimistic"), "TIC"].astype(int).unique())

plt.rcParams["lines.linewidth"] = 2
fig, axes = plt.subplots(1, 2, figsize=(30, 10))

T_GRID = np.linspace(2800, 7300, 1000)
XLIM = (2.2, 0)
YLIM = (3000, 6900)


def panel(ax, df, title):
    opt_in, opt_out = hz_edge(T_GRID, 0), hz_edge(T_GRID, 3)
    con_in, con_out = hz_edge(T_GRID, 1), hz_edge(T_GRID, 2)
    ax.fill_betweenx(T_GRID, opt_in, opt_out, color="tab:green", alpha=0.2, zorder=11)
    ax.fill_betweenx(T_GRID, con_in, con_out, color="tab:green", alpha=0.4, zorder=11)

    d = df[
        (df.Teff >= YLIM[0]) & (df.Teff <= YLIM[1])
        & (df.S >= XLIM[1]) & (df.S <= XLIM[0])
    ].copy()

    rad = d["Radius_planet"].fillna(1.0).clip(0.3, 4).to_numpy()
    order = np.flip(np.argsort(rad))
    size = 120 * rad[order] ** 2

    scp = ax.scatter(
        d["S"].to_numpy()[order], d["Teff"].to_numpy()[order],
        c=d["conf"].to_numpy()[order],
        cmap="RdBu", vmin=0, vmax=1, edgecolors="k",
        s=size, alpha=0.95, zorder=100,
    )

    fail = d[d["passed_all_tests"] == False]
    if len(fail):
        rf = fail["Radius_planet"].fillna(1.0).clip(0.3, 4).to_numpy()
        ax.scatter(
            fail["S"], fail["Teff"],
            facecolors="none", edgecolors="red",
            linestyles="--", linewidths=2.0,
            s=120 * rf**2 * 1.9, zorder=150,
        )

    for y, ls, a in [
        (6000, "--", 0.3), (5300, "--", 0.3), (3900, "--", 0.3),
        (4800, "-", 0.5), (6300, "-", 0.5),
    ]:
        ax.plot([XLIM[0], XLIM[1]], [y, y], "k", linestyle=ls, alpha=a, zorder=10)

    ax.text(-0.01, 6040, "F", alpha=0.5, fontsize=24)
    ax.text(-0.01, 5340, "G", alpha=0.5, fontsize=24)
    ax.text(-0.01, 3940, "K", alpha=0.5, fontsize=24)
    ax.text(-0.01, YLIM[0] + 50, "M", alpha=0.5, fontsize=24)
    ax.text(0.12, 4840, "4800K", alpha=0.5, fontsize=16)
    ax.text(0.12, 6340, "6300K", alpha=0.5, fontsize=16)

    ax.text(1.1, 6700, "Conservative HZ", color="tab:green", fontsize=22)
    ax.text(1.8, 6700, "Optimistic HZ", color="tab:green", fontsize=22)

    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_xlabel(r"Instellation Flux [$I_\oplus$]", fontsize=32)
    ax.set_title(title, fontsize=28)
    ax.tick_params(labelsize=18)
    return scp


scp_left = panel(axes[0], known, "Known TOIs")
scp_right = panel(axes[1], new, "New Candidates")
axes[0].set_ylabel("Effective Temperature [K]", fontsize=32)

hl = new[new.TIC.isin(HL_TICS)]
axes[1].scatter(
    hl["S"], hl["Teff"], marker="*", s=1600, color="gold",
    edgecolors="k", linewidths=1.5, zorder=200, label="Newly-identified HZ",
)
for _, r in hl.iterrows():
    axes[1].annotate(
        f"TIC {int(r.TIC)}", (r.S, r.Teff),
        xytext=(10, 10), textcoords="offset points", fontsize=15,
        color="black", zorder=210,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="gold", alpha=0.85),
    )

kw = dict(
    prop="sizes", num=[1, 2, 4], color="grey",
    fmt=r"{x:.0f} $R_\oplus$", markeredgecolor="k",
    func=lambda s: np.sqrt(s / 120),
)
axes[1].legend(*scp_right.legend_elements(**kw), loc=[0.72, 0.78], fontsize=18)

cbh = fig.colorbar(scp_right, ax=axes, pad=0.02, fraction=0.04)
cbh.ax.set_ylabel("1 - NST p-value", fontsize=24)
cbh.ax.tick_params(labelsize=18)

plt.savefig(OUT, dpi=110, bbox_inches="tight")
plt.close()
print(f"Wrote {OUT}")


def _count(df, label):
    opt = in_hz(df, "optimistic"); con = in_hz(df, "conservative")
    pa = df["passed_all_tests"] == True
    print(
        f"{label}: n={len(df)}, "
        f"optimistic HZ={int(opt.sum())} (passed={int((opt & pa).sum())}), "
        f"conservative HZ={int(con.sum())} (passed={int((con & pa).sum())})"
    )


_count(known, "Known")
_count(new, "New")
print(f"Highlighted (new in optimistic HZ): {HL_TICS}")
