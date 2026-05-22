"""Regenerate new_outliers.png with the Teff axes capped at 10000 K.
OoD scores/flags are read from tois_new.csv (unchanged)."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/global/u2/j/julius/TESS corrected"
TEFF_MAX = 10000.0

known = pd.read_csv(f"{BASE}/tois.csv")
new = pd.read_csv(f"{BASE}/tois_new.csv")

# known-population features
mult = known.groupby("TIC")["TIC"].transform("size").clip(lower=1)
kR, kP, kT, kM = known["Radius_planet"], known["Period"], known["Teff"], mult

# new candidates (host Teff mapped from the known catalog where missing, as the OoD used it)
teff_map = known.dropna(subset=["Teff"]).groupby("TIC")["Teff"].mean()
nR, nP, nM = new["Radius_planet"], new["Period"], new["Multiplicity"]
nT = new["Teff"].fillna(new["TIC"].map(teff_map))
interesting = new["Interesting"].fillna(False).astype(bool)
errp, errm = new["Radius_planet_errp"], new["Radius_planet_errm"]

rng = np.random.default_rng(0)
def jit(x, w=0.10):
    return x + rng.uniform(-w, w, size=len(x))

GREY = dict(s=8, c="0.6", alpha=0.30, linewidths=0, label="Known TESS planets")
NEW = dict(s=16, c="#ff7f0e", alpha=0.85, linewidths=0, label="New candidates")
STAR = dict(marker="*", s=150, c="#a11212", edgecolors="k", linewidths=0.4, zorder=5,
            label="Interesting (OoD)")

fig, ax = plt.subplots(2, 2, figsize=(11, 8))
fig.suptitle("New candidates vs. known TESS planets — density-based OoD outliers (★)",
             fontsize=12)

# (0,0) Period vs radius
a = ax[0, 0]
a.scatter(kP, kR, **GREY); a.scatter(nP, nR, **NEW)
a.scatter(nP[interesting], nR[interesting], **STAR)
a.set_xscale("log"); a.set_yscale("log")
a.set_xlabel("Period [d]"); a.set_ylabel("Planet radius [R⊕]")
a.legend(loc="upper left", fontsize=7, framealpha=0.9)

# (0,1) Teff vs radius  (Teff capped)
a = ax[0, 1]
a.scatter(kT, kR, **GREY); a.scatter(nT, nR, **NEW)
a.scatter(nT[interesting], nR[interesting], **STAR)
a.set_yscale("log"); a.set_xlim(2500, TEFF_MAX)
a.set_xlabel("Host T$_{eff}$ [K]"); a.set_ylabel("Planet radius [R⊕]")

# (1,0) Period vs Teff  (Teff capped)
a = ax[1, 0]
a.scatter(kP, kT, **GREY); a.scatter(nP, nT, **NEW)
a.scatter(nP[interesting], nT[interesting], **STAR)
a.set_xscale("log"); a.set_ylim(2500, TEFF_MAX)
a.set_xlabel("Period [d]"); a.set_ylabel("Host T$_{eff}$ [K]")

# (1,1) Multiplicity vs radius (with asymmetric radius error bars on interesting)
a = ax[1, 1]
a.scatter(jit(kM), kR, **GREY)
a.scatter(jit(nM), nR, **NEW)
mi = interesting.values
a.errorbar(nM[mi], nR[mi], yerr=[errm[mi].fillna(0), errp[mi].fillna(0)],
           fmt="*", ms=13, mfc="#a11212", mec="k", mew=0.4, ecolor="#a11212",
           elinewidth=1.0, capsize=2, zorder=5)
a.set_yscale("log")
a.set_xlabel("System multiplicity"); a.set_ylabel("Planet radius [R⊕]")

fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(f"{BASE}/new_outliers.png", dpi=110)
print("wrote new_outliers.png  |  interesting:", int(mi.sum()),
      " Teff capped at", TEFF_MAX, "K")
