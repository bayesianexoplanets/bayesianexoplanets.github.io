"""apply_rerun.py — merge the KnownRun rerun results into the website catalog.

For the rerun TICs (CleanUp.md "old catalog values" list), replace the
pipeline-recomputed columns of ``TESS corrected/tois.csv`` with the values from
``KnownRun/known_candidates/{tic}.csv`` (produced by
``pipeline/known_planets_vetting.py``), and refresh the per-planet plot the
website shows (``plots/{TIC}/{known_idx}.jpg``).

Only columns the rerun authoritatively recomputes are touched; NST columns
(μ, σ, nst_samples), stellar columns (Mass, Radius, logg, FEH, Teff) and
``Has Visible TTVs`` are kept. ``log10(p value)`` is recomputed from the new SNR
and the kept μ/σ. ``passed_all_tests`` / ``failed_tests`` are recomputed from the
new direct diagnostics (same rule as recompute_failed_tests.py).

Matching (robust to mid-list vetting failures): each website row of a TIC, taken
in file order (that order *is* the frontend ``known_idx`` — index.html L573), is
matched to a rerun planet by TOI. A rerun row is tied to its TOI and to its
Stage-A plot file via ``catalog_period`` -> the read_known_planets_tess planet
order (which names the plot files and indexes the ExoFOP TOI) — NOT via event_id,
which shifts when a middle planet's vetting throws.

Verified conventions: Duration = 2*Tau, Epoch = t_start + Phase, Radius_planet_err{p,m} = radius{p,m},
log10(p) = interp of the per-star Singh-Maddala posterior-mean-SF grid (sm_sf_grid) at the SNR.

Dry-run by default (prints planned changes); pass --apply to write tois.csv + plots.
"""
import os
import sys
import glob

os.chdir("/global/u2/j/julius")
sys.path.insert(0, "/global/u2/j/julius")
import numpy as np
import pandas as pd
sys.path.insert(0, "/global/u2/j/julius/TESS corrected")
import sm_pvalue
from scipy.stats import norm
from TESS.load_tess import StarInfo_tess, read_known_planets_tess

HERE = "/global/u2/j/julius/TESS corrected"
TOIS = os.path.join(HERE, "tois.csv")
KNOWN = "/pscratch/sd/j/julius/exoprob/KnownRun/known_candidates"
PLOTS_SRC = "/pscratch/sd/j/julius/exoprob/KnownRun/plots"
PLOTS_DST = os.path.join(HERE, "plots")
STARS = "/pscratch/sd/j/julius/exoprob/FullRun/stars"
BULK = "/pscratch/sd/j/julius/Bulk Download/Data"
LN10 = np.log(10.0)

# reconstructed vetting rule (recompute_failed_tests.py)
SPURIOUS_MAX, SNRD_MIN, SNR_OVERRIDE, NTRANSITS_MIN, HARMONICS_SNRD = 0.75, 1e-2, 50.0, 3, 0.95


def has_sharp_peak(tic):
    try:
        df = pd.read_csv(os.path.join(STARS, f"{tic}.csv"), sep="\t")
        return bool(df["has_sharp_peak"].iloc[0]) if not df.empty else False
    except Exception:
        return False


def planet_order(tic):
    """(list of (plot_index, period, TOI) in read_known order, t_start) — the
    order Stage A names the plot files and the ExoFOP TOI each planet maps to."""
    sd = np.load(os.path.join(BULK, f"StellarData_{tic}.npy"))
    sdc = sd.copy()
    if np.isnan(sdc[3]):
        sdc[3] = 1
    if np.isnan(sdc[-1]):
        sdc[-1] = 1
        sdc[-2] = 1
    star = StarInfo_tess(*sdc[1:])
    pp = read_known_planets_tess(star, sd[0], False)
    e = pd.read_csv("TESS/tois.csv")
    e = e[e["TIC ID"] == tic]
    e = e[e["TFOPWG Disposition"].isin(["PC", "CP", "KP", "APC"])]
    eper = e["Period (days)"].to_numpy(float)
    etoi = e["TOI"].to_numpy()
    order = []
    for i, p in enumerate(pp):
        per = float(p.params[0])
        toi = str(etoi[int(np.nanargmin(np.abs(eper - per)))]) if len(eper) else ""
        order.append((i, per, toi))
    return order, float(sd[0])


def failures(spurious1, snrd, snr, ntr, sharp):
    f = []
    if spurious1 >= SPURIOUS_MAX:
        f.append("spurious")
    if snrd <= SNRD_MIN and snr <= SNR_OVERRIDE:
        f.append("snrd")
    if ntr < NTRANSITS_MIN:
        f.append("ntransits")
    if sharp and snrd > HARMONICS_SNRD and snr <= SNR_OVERRIDE:
        f.append("sharp_freq")
    return f


def main(apply):
    w = pd.read_csv(TOIS)
    rerun_tics = sorted({int(os.path.basename(p).split(".")[0]) for p in glob.glob(os.path.join(KNOWN, "*.csv"))})
    print(f"website rows={len(w)} | rerun TICs with a CSV={len(rerun_tics)}")

    n_rows = n_plots = 0
    skipped, unmatched = [], []

    for tic in rerun_tics:
        rr = pd.read_csv(os.path.join(KNOWN, f"{tic}.csv"), sep="\t")
        wsub = w[w["TIC"] == tic]
        if len(rr) == 0:
            skipped.append((tic, "rerun produced 0 vetted planets"))
            continue
        if len(wsub) == 0:
            skipped.append((tic, "TIC not in website catalog"))
            continue
        sharp = has_sharp_peak(tic)
        order, ts = planet_order(tic)

        def lookup(cp):  # rerun catalog_period -> (plot_index, TOI)
            if not order:
                return None, None
            k = int(np.argmin([abs(per - cp) for (_, per, _) in order]))
            idx, per, toi = order[k]
            return (idx, toi) if abs(per - cp) < 1e-4 * max(per, 1.0) else (idx, toi)

        # rerun row keyed by its true TOI, carrying its Stage-A plot index
        by_toi = {}
        for _, r in rr.iterrows():
            plot_idx, toi = lookup(float(r["catalog_period"]))
            if toi is not None:
                by_toi[toi] = (r, plot_idx)

        for known_idx, (wi, wrow) in enumerate(wsub.iterrows()):
            toi = str(wrow["TOI"])
            if toi not in by_toi:
                unmatched.append((tic, toi))
                continue
            r, plot_idx = by_toi[toi]
            snr = float(r["SNR"]); snrd = float(r["snrd_pvalue"])
            sp1 = float(r["spurious1"]); ntr = int(r["num_available_transits"])
            fails = failures(sp1, snrd, snr, ntr, sharp)
            passed = len(fails) == 0
            # Singh-Maddala posterior-mean-SF p-value: interpolate the star's tabulated log10-SF grid.
            log10p = sm_pvalue.log10p_from_grid(wrow["sm_sf_grid"], snr)
            if not np.isfinite(log10p):
                log10p = wrow["log10(p value)"]

            upd = {
                "Period": float(r["period"]), "Phase": float(r["phase"]), "Tau": float(r["tau"]),
                "Duration": 2.0 * float(r["tau"]), "Epoch": ts + float(r["phase"]), "SNR": snr,
                "err_Period": float(r["err_period"]), "err_Epoch": float(r["err_phase"]),
                "err_Duration": 2.0 * float(r["err_tau"]),
                "Radius_planet": float(r["radius"]),
                "Radius_planet_errp": float(r["radiusp"]), "Radius_planet_errm": float(r["radiusm"]),
                "Number of Valid Transits": ntr, "log10(p value)": log10p,
                "passed_all_tests": bool(passed), "failed_tests": "|".join(fails),
            }
            if apply:
                for k, v in upd.items():
                    w.at[wi, k] = v
            else:
                print(f"  TIC {tic} idx{known_idx} TOI {toi} (plot {plot_idx}): "
                      f"P {wrow['Period']:.4f}->{upd['Period']:.4f}  SNR {wrow['SNR']:.2f}->{snr:.2f}  "
                      f"ntr {wrow['Number of Valid Transits']}->{ntr}  "
                      f"pass {wrow['passed_all_tests']}->{passed} [{upd['failed_tests']}]")
            n_rows += 1

            src = os.path.join(PLOTS_SRC, str(tic), f"{plot_idx}.png")
            dst = os.path.join(PLOTS_DST, str(tic), f"{known_idx}.jpg")
            if os.path.exists(src):
                if apply:
                    _to_jpg(src, dst)
                n_plots += 1
            else:
                print(f"    WARN missing plot {src}")

    print(f"\nrows to update: {n_rows} | plots to refresh: {n_plots}")
    if skipped:
        print("skipped TICs:", skipped)
    if unmatched:
        print("unmatched website rows (kept as-is — rerun didn't vet this planet):", unmatched)
    if apply:
        w.to_csv(TOIS, index=False)
        print(f"\nWROTE {TOIS}  ({len(w)} rows)")
    else:
        print("\nDRY RUN — pass --apply to write tois.csv and plots")


def _to_jpg(src, dst, max_w=900, quality=88):
    from PIL import Image
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    im = Image.open(src).convert("RGB")
    if im.width > max_w:
        im = im.resize((max_w, int(im.height * max_w / im.width)), Image.LANCZOS)
    im.save(dst, "JPEG", quality=quality)


if __name__ == "__main__":
    main(apply="--apply" in sys.argv)
