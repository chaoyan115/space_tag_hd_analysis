#!/usr/bin/env python3
"""Fragments-per-capturing-spot QC for SPACE-Tag-HD binned fragment files.

Reads a binned fragments BED (gzip; col4 = spot barcode, col5 = count) and the
matching SpaceRanger tissue_positions.parquet, restricts to in_tissue spots, and
writes (a) summary stats and (b) a spatial map coloured by log10 fragments/spot.

Fragment count per spot = sum of col5 over fragment records with that barcode.
"""
import argparse
import gzip
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def count_frags_per_spot(frag_path):
    counts = defaultdict(int)
    opener = gzip.open if frag_path.endswith(".gz") else open
    with opener(frag_path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                if len(parts) == 4:
                    counts[parts[3]] += 1
                continue
            bc = parts[3]
            try:
                counts[bc] += int(parts[4])
            except ValueError:
                counts[bc] += 1
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frag", required=True, help="binned fragments bed.gz")
    ap.add_argument("--positions", required=True, help="tissue_positions.parquet")
    ap.add_argument("--outprefix", required=True, help="output prefix for .pdf/.png")
    ap.add_argument("--stats", required=True, help="output stats txt path")
    ap.add_argument("--sample", default="sample")
    args = ap.parse_args()

    pos = pd.read_parquet(args.positions)
    # normalize column names (SpaceRanger HD parquet uses these)
    bc_col = "barcode" if "barcode" in pos.columns else pos.columns[0]
    pos = pos[pos["in_tissue"] == 1].copy()
    pos = pos.set_index(bc_col)

    counts = count_frags_per_spot(args.frag)
    pos["frags"] = pos.index.map(lambda b: counts.get(b, 0)).astype(float)

    in_tissue_frags = pos["frags"].values
    n_tissue = len(pos)
    n_with = int((in_tissue_frags > 0).sum())
    total = float(in_tissue_frags.sum())
    median = float(np.median(in_tissue_frags))
    mean = float(np.mean(in_tissue_frags))
    # fraction of all fragments that fall on in-tissue spots
    total_all = float(sum(counts.values()))
    frac_in_tissue = total / total_all if total_all > 0 else float("nan")

    with open(args.stats, "w") as fh:
        fh.write(f"SAMPLE\t{args.sample}\n")
        fh.write(f"N_TISSUE_SPOTS\t{n_tissue}\n")
        fh.write(f"N_SPOTS_WITH_FRAGS\t{n_with}\n")
        fh.write(f"FRAC_SPOTS_WITH_FRAGS\t{n_with / n_tissue if n_tissue else float('nan'):.4f}\n")
        fh.write(f"TOTAL_FRAGS_IN_TISSUE\t{total:.0f}\n")
        fh.write(f"TOTAL_FRAGS_ALL_SPOTS\t{total_all:.0f}\n")
        fh.write(f"FRAC_FRAGS_IN_TISSUE\t{frac_in_tissue:.4f}\n")
        fh.write(f"MEDIAN_FRAGS_PER_SPOT\t{median:.1f}\n")
        fh.write(f"MEAN_FRAGS_PER_SPOT\t{mean:.1f}\n")

    # spatial map: pixel coords coloured by log10(frags+1)
    xcol = "pxl_col_in_fullres"
    ycol = "pxl_row_in_fullres"
    if xcol not in pos.columns or ycol not in pos.columns:
        # fall back to array grid
        xcol, ycol = "array_col", "array_row"

    logf = np.log10(pos["frags"].values + 1.0)
    fig, ax = plt.subplots(figsize=(7, 7))
    sc = ax.scatter(pos[xcol].values, -pos[ycol].values, c=logf, s=1.0,
                    cmap="viridis", linewidths=0, rasterized=True)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"{args.sample}  fragments/spot (log10)\n"
                 f"n={n_tissue} in-tissue spots, median={median:.0f}, total={total:.0f}")
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("log10(fragments + 1)")
    fig.tight_layout()
    fig.savefig(args.outprefix + ".pdf", dpi=200)
    fig.savefig(args.outprefix + ".png", dpi=150)
    print(f"[qc_fragspot] {args.sample}: n_tissue={n_tissue} median={median:.0f} "
          f"total={total:.0f} frac_in_tissue={frac_in_tissue:.3f}", file=sys.stderr)


if __name__ == "__main__":
    main()
