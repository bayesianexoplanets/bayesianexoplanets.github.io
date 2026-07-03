"""merge_ephem68.py -- merge the RerunEphem68 "corrected ephemeris" re-vetting into tois.csv.

For the 68 TICs whose catalog (phase, tau) sat on a local matched-filter maximum (found by
the TESS TTVs catalog run, visually confirmed in TESS TTVs/corrections_review/), replace the
pipeline-recomputed columns of TESS corrected/tois.csv with the values from
RerunEphem68/known_candidates/{tic}.csv (matched by TOI; sibling planets of the same star are
also refreshed since the star was re-vetted as a whole), refresh the per-planet website plot,
and recompute log10(p) from the kept per-star Singh-Maddala sm_sf_grid + passed/failed tests.

Dry-run by default; pass --apply to write tois.csv + plots.
"""
import os, sys, glob
os.chdir("/global/u2/j/julius"); sys.path.insert(0, "/global/u2/j/julius")
import numpy as np, pandas as pd
sys.path.insert(0, "/global/u2/j/julius/TESS corrected")
import sm_pvalue

HERE = "/global/u2/j/julius/TESS corrected"
TOIS = os.path.join(HERE, "tois.csv")
KNOWN = "/pscratch/sd/j/julius/exoprob/RerunEphem68/known_candidates"
PLOTS_SRC = "/pscratch/sd/j/julius/exoprob/RerunEphem68/plots"
PLOTS_DST = os.path.join(HERE, "plots")
STARS = "/pscratch/sd/j/julius/exoprob/FullRun/stars"
BULK = "/pscratch/sd/j/julius/Bulk Download/Data"
SPURIOUS_MAX, SNRD_MIN, SNR_OVERRIDE, NTRANSITS_MIN, HARMONICS_SNRD = 0.75, 1e-2, 50.0, 3, 0.95


def has_sharp_peak(tic):
    try:
        df = pd.read_csv(os.path.join(STARS, f"{tic}.csv"), sep="\t")
        return bool(df["has_sharp_peak"].iloc[0]) if not df.empty else False
    except Exception:
        return False


def failures(sp1, snrd, snr, ntr, sharp):
    f = []
    if sp1 >= SPURIOUS_MAX: f.append("spurious")
    if snrd <= SNRD_MIN and snr <= SNR_OVERRIDE: f.append("snrd")
    if ntr < NTRANSITS_MIN: f.append("ntransits")
    if sharp and snrd > HARMONICS_SNRD and snr <= SNR_OVERRIDE: f.append("sharp_freq")
    return f


def toi_key(x):
    try: return f"{float(x):.2f}"
    except Exception: return str(x)


def main(apply):
    w = pd.read_csv(TOIS)
    tics = sorted({int(os.path.basename(p).split(".")[0]) for p in glob.glob(os.path.join(KNOWN, "*.csv"))})
    print(f"website rows={len(w)} | RerunEphem68 TICs={len(tics)}")
    n_rows = n_plots = 0; unmatched = []

    for tic in tics:
        rr = pd.read_csv(os.path.join(KNOWN, f"{tic}.csv"), sep="\t")
        if len(rr) == 0:
            print(f"  TIC {tic}: 0 vetted planets -- SKIP"); continue
        sd = np.load(os.path.join(BULK, f"StellarData_{tic}.npy")); t_start = float(sd[0])
        sharp = has_sharp_peak(tic)
        by_toi = {toi_key(r["TOI"]): r for _, r in rr.iterrows()}

        wsub = w[w["TIC"] == tic]
        for known_idx, (wi, wrow) in enumerate(wsub.iterrows()):
            k = toi_key(wrow["TOI"])
            if k not in by_toi:
                unmatched.append((tic, k)); continue
            r = by_toi[k]
            snr = float(r["SNR"]); snrd = float(r["snrd_pvalue"]); sp1 = float(r["spurious1"])
            ntr = int(r["num_available_transits"])
            fails = failures(sp1, snrd, snr, ntr, sharp); passed = len(fails) == 0
            upd = {"Period": float(r["period"]), "Phase": float(r["phase"]), "Tau": float(r["tau"]),
                   "Duration": 2.0 * float(r["tau"]), "Epoch": t_start + float(r["phase"]), "SNR": snr,
                   "err_Period": float(r["err_period"]), "err_Epoch": float(r["err_phase"]),
                   "err_Duration": 2.0 * float(r["err_tau"]), "Radius_planet": float(r["radius"]),
                   "Radius_planet_errp": float(r["radiusp"]), "Radius_planet_errm": float(r["radiusm"]),
                   "Number of Valid Transits": ntr,
                   "passed_all_tests": bool(passed), "failed_tests": "|".join(fails)}
            lp = sm_pvalue.log10p_from_grid(wrow["sm_sf_grid"], snr)
            upd["log10(p value)"] = lp if np.isfinite(lp) else wrow["log10(p value)"]

            if apply:
                for c, v in upd.items(): w.at[wi, c] = v
            else:
                print(f"  TIC {tic} TOI {k} (idx{known_idx}, plot {int(r['plot_index'])}): "
                      f"P {wrow['Period']:.4f}->{upd['Period']:.4f}  tau {wrow['Tau']:.4f}->{upd['Tau']:.4f}  "
                      f"SNR {wrow['SNR']:.2f}->{snr:.2f}  ntr {wrow['Number of Valid Transits']}->{ntr}  "
                      f"pass {wrow['passed_all_tests']}->{passed} [{upd['failed_tests']}]")
            n_rows += 1
            src = os.path.join(PLOTS_SRC, str(tic), f"{int(r['plot_index'])}.png")
            dst = os.path.join(PLOTS_DST, str(tic), f"{known_idx}.jpg")
            if os.path.exists(src):
                if apply: _to_jpg(src, dst)
                n_plots += 1
            else:
                print(f"    WARN missing plot {src}")

    print(f"\nrows updated: {n_rows} | plots refreshed: {n_plots}")
    if unmatched: print("unmatched website rows (left as-is):", unmatched)
    if apply:
        w.to_csv(TOIS, index=False); print(f"WROTE {TOIS} ({len(w)} rows)")
    else:
        print("DRY RUN -- pass --apply to write")


def _to_jpg(src, dst, max_w=900, quality=88):
    from PIL import Image
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    im = Image.open(src).convert("RGB")
    if im.width > max_w:
        im = im.resize((max_w, int(im.height * max_w / im.width)), Image.LANCZOS)
    im.save(dst, "JPEG", quality=quality)


if __name__ == "__main__":
    main(apply="--apply" in sys.argv)
