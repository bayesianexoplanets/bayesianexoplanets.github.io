"""refresh_masked_figures.py — re-render the figures of every website row whose star carries another catalogued signal.

`pipeline/pretty_plot.py:load_data` now masks the star's other signals (curated catalog, few-transit rows, new
candidates) before folding, as the blind search masks known planets. Only rows with at least one other signal change.
Each figure is rendered from the website row's own Period/Phase/Tau, so table and figure agree, then written in place:
``plots/{TIC}/{known_idx}.jpg`` (catalog) and ``plots_new/{TIC}/{_cand_idx}.jpg`` (new candidates).

Usage:
    python refresh_masked_figures.py render [workers]   # PNGs to scratch, resumable (an existing PNG is skipped)
    python refresh_masked_figures.py install [--apply]  # PNG -> JPG into the site tree (dry run by default)
"""
import os
import sys
import time
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/exoplanets/TESS corrected"
sys.path.insert(0, HERE); sys.path.insert(0, "/global/u2/j/julius/exoplanets")
from constants import scratch
from merge_known_run import _to_jpg

OUT = scratch + "MaskedFigures_20261001/"
SUMMARY = OUT + "summary.csv"
PAGES = (("tess", "tois.csv", "plots"), ("new", "tois_new.csv", "plots_new"))


def tasks():
    """One task per website row whose star has another catalogued signal: (page, TIC, slot, period, phase, tau, snr)."""
    from pipeline.pretty_plot import other_signals
    rows = []
    for page, table_name, _ in PAGES:
        table = pd.read_csv(os.path.join(HERE, table_name))
        slots = table.groupby("TIC").cumcount() if page == "tess" else table["_cand_idx"]
        for slot, r in zip(slots, table.itertuples()):
            params = np.array([r.Period, r.Phase, r.Tau], dtype=float)
            if np.all(np.isfinite(params)) and other_signals(int(r.TIC), params[0]):
                rows.append((page, int(r.TIC), int(slot), *params, float(r.SNR)))
    return rows


def png_path(page, tic, slot):
    return OUT + f"{page}/{tic}_{slot}.png"


def one(task):
    """Render one figure; never raises."""
    warnings.filterwarnings("ignore")
    import matplotlib
    matplotlib.use("Agg")
    from pipeline.pretty_plot import pretty_transit_plot
    page, tic, slot, period, phase, tau, snr = task
    t0 = time.time()
    try:
        summary = pretty_transit_plot(tic, period, phase, tau, png_path(page, tic, slot), snr=snr)
        return dict(page=page, TIC=tic, slot=slot, model=str(summary.get("model", "none")), seconds=time.time() - t0, error="")
    except Exception as e:
        return dict(page=page, TIC=tic, slot=slot, model="none", seconds=time.time() - t0, error=repr(e)[:200])


def render(workers):
    for page, _, _ in PAGES:
        os.makedirs(OUT + page, exist_ok=True)
    todo = [t for t in tasks() if not os.path.exists(png_path(*t[:3]))]
    print(len(todo), "figures to render,", workers, "workers", flush=True)
    t_start = time.time()
    with Pool(workers, maxtasksperchild=20) as pool:
        for i, result in enumerate(pool.imap_unordered(one, todo, chunksize=1)):
            pd.DataFrame([result]).to_csv(SUMMARY, mode="a", header=not os.path.exists(SUMMARY), index=False)
            if (i + 1) % 50 == 0 or result["error"]:
                print(f"[{i + 1}/{len(todo)}] {result['page']} TIC {result['TIC']} #{result['slot']} {result['seconds']:.0f}s "
                      f"{result['error']} ({(time.time() - t_start) / 60:.1f} min)", flush=True)
    print("done in %.1f min" % ((time.time() - t_start) / 60), flush=True)


def install(apply):
    destination = {page: os.path.join(HERE, folder) for page, _, folder in PAGES}
    all_tasks = tasks()
    ready = [t for t in all_tasks if os.path.exists(png_path(*t[:3]))]
    print(f"{len(ready)} of {len(all_tasks)} figures rendered", flush=True)
    for page, _, _ in PAGES:
        print(f"  {page}: {sum(t[0] == page for t in ready)}")
    if not apply:
        print("dry run; pass --apply to write the JPGs")
        return
    for page, tic, slot, *_ in ready:
        _to_jpg(png_path(page, tic, slot), os.path.join(destination[page], str(tic), f"{slot}.jpg"))
    print("installed", len(ready))


if __name__ == "__main__":
    stage = sys.argv[1]
    if stage == "render":
        render(int(sys.argv[2]) if len(sys.argv) > 2 else 4)
    elif stage == "install":
        install("--apply" in sys.argv)
    else:
        sys.exit(__doc__)
