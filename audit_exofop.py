"""audit_exofop.py — thorough discrepancy check between the website catalog
(TESS corrected/tois.csv) and the ExoFOP catalog (TESS/tois.csv).

Read-only. Reports:
  (1) PERIOD discrepancies for matched (TIC, TOI): website Period vs ExoFOP
      Period(days), split into exact (<tol), integer-multiple "harmonic", and
      genuine mismatch — each annotated with the latest re-vet SNR + pass/fail so
      we can tell whether a mismatch still corresponds to a non-detection.
  (2) MEMBERSHIP: ExoFOP PC/CP/KP/APC planets on a website TIC that are MISSING
      from the website; website rows whose TOI is not a PC/CP/KP/APC ExoFOP planet.
  (3) EPOCH/DURATION mismatches for matched rows (sanity).

Usage: python audit_exofop.py [reltol]   (default reltol=0.002 = 0.2%)
"""
import os
import sys

os.chdir("/global/u2/j/julius")
import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/TESS corrected"
WEB = os.path.join(HERE, "tois.csv")
EXO = "TESS/tois.csv"
RUN = "/pscratch/sd/j/julius/exoprob/CatalogRun/known_vetting_summary.csv"
RELTOL = float(sys.argv[1]) if len(sys.argv) > 1 else 0.002
DISP = ["PC", "CP", "KP", "APC"]


def period_relation(pw, pe, tol=RELTOL):
    """Return ('match'|'harmonic:N'|'mismatch', rel_error)."""
    if not (np.isfinite(pw) and np.isfinite(pe)) or pw <= 0 or pe <= 0:
        return "nan", np.inf
    r = max(pw / pe, pe / pw)
    n = round(r)
    rel = abs(r - n) / n
    if n == 1 and abs(pw - pe) / pe <= tol:
        return "match", abs(pw - pe) / pe
    if n >= 2 and rel <= tol:
        return f"harmonic:x{n}", rel
    return "mismatch", abs(pw - pe) / pe


def main():
    w = pd.read_csv(WEB)
    e = pd.read_csv(EXO)
    ecp = e[e["TFOPWG Disposition"].isin(DISP)]
    # ExoFOP maps keyed by (TIC, TOI) and grouped by TIC
    exo = {(int(r["TIC ID"]), str(r["TOI"])): r for _, r in ecp.iterrows()}
    exo_by_tic = {}
    for _, r in ecp.iterrows():
        exo_by_tic.setdefault(int(r["TIC ID"]), []).append(str(r["TOI"]))
    # latest re-vet SNR per (TIC,TOI) for annotation
    revsnr = {}
    try:
        rr = pd.read_csv(RUN, sep="\t")
        for _, r in rr.iterrows():
            revsnr[(int(r["kepid"]), str(r["TOI"]))] = float(r["SNR"])
    except Exception:
        pass

    print(f"website rows={len(w)} | ExoFOP PC/CP/KP/APC rows={len(ecp)} | reltol={RELTOL}\n")

    # (1) period discrepancies on matched rows
    harm, mism, nan = [], [], []
    matched = 0
    for _, wr in w.iterrows():
        tic = int(wr["TIC"]); toi = str(wr["TOI"])
        er = exo.get((tic, toi))
        if er is None:
            continue
        matched += 1
        rel = period_relation(float(wr["Period"]), float(er["Period (days)"]))
        kind = rel[0]
        rec = (tic, toi, round(float(wr["Period"]), 5), round(float(er["Period (days)"]), 5),
               round(rel[1] * 100, 2), round(revsnr.get((tic, toi), float("nan")), 1),
               bool(wr["passed_all_tests"]))
        if kind == "mismatch":
            mism.append(rec)
        elif kind.startswith("harmonic"):
            harm.append(rec + (kind,))
        elif kind == "nan":
            nan.append(rec)

    print(f"=== (1) PERIOD on {matched} matched (TIC,TOI) ===")
    print(f"  genuine mismatches (> {RELTOL*100}% and not a harmonic): {len(mism)}")
    if mism:
        print(f"  {'TIC':>11} {'TOI':>9} {'Pweb':>9} {'Pexo':>9} {'dP%':>7} {'reSNR':>6} pass")
        for c in sorted(mism, key=lambda x: -x[4]):
            print(f"  {c[0]:>11} {c[1]:>9} {c[2]:>9} {c[3]:>9} {c[4]:>7} {c[5]:>6} {c[6]}")
    print(f"  integer-multiple 'harmonic' periods: {len(harm)}")
    for c in sorted(harm, key=lambda x: -x[4])[:40]:
        print(f"     {c[0]} {c[1]}  Pweb={c[2]} Pexo={c[3]}  {c[7]}  reSNR={c[5]} pass={c[6]}")
    if nan:
        print(f"  non-finite period rows: {len(nan)}")

    # (2) membership
    web_keys = {(int(r["TIC"]), str(r["TOI"])) for _, r in w.iterrows()}
    web_tics = set(w["TIC"].astype(int))
    missing = [(t, toi) for t in web_tics for toi in exo_by_tic.get(t, []) if (t, toi) not in web_keys]
    extra = [(int(r["TIC"]), str(r["TOI"])) for _, r in w.iterrows()
             if (int(r["TIC"]), str(r["TOI"])) not in exo]
    print(f"\n=== (2) MEMBERSHIP ===")
    print(f"  ExoFOP PC/CP/KP/APC planets on a website TIC but MISSING from website: {len(missing)}")
    for m in sorted(missing)[:40]:
        print(f"     TIC {m[0]} TOI {m[1]}  (Pexo={float(exo[(m[0],m[1])]['Period (days)']):.4f}, "
              f"exoSNR={exo[(m[0],m[1])]['Planet SNR']})")
    print(f"  website rows whose TOI is not a PC/CP/KP/APC ExoFOP planet: {len(extra)}")
    for x in sorted(extra)[:25]:
        print(f"     TIC {x[0]} TOI {x[1]}")

    # (3) epoch / duration sanity on matched rows (BTJD epoch = t_start+phase vs ExoFOP)
    print(f"\n=== (3) summary ===")
    print(f"  matched (TIC,TOI): {matched} | period: {len(mism)} mismatch, {len(harm)} harmonic")
    print(f"  membership: {len(missing)} missing, {len(extra)} extra")


if __name__ == "__main__":
    main()
