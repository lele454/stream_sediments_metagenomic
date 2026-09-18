#!/usr/bin/env python3
"""
Figure: alpha diversity of the stream communities
---------------------------------------------------
Boxplots of Hill numbers q = 0, 1 and 2 with Shannon entropy, compared between
preserved and degraded streams.

Each panel carries the Mann-Whitney p value for that index. Shannon entropy H and the Hill number of
order one are the same quantity on different scales, q1 = exp(H), so they will
always tell the same story and are shown together only because both are
conventionally reported. The information is in the series across orders: q0
counts genera regardless of abundance, q1 weights them by frequency, q2
weights the common ones more heavily still. A difference at q0 but not at q2
means rare genera differ while the dominant ones do not, which is a different
result from a difference across all three.

Richness is sensitive to sequencing depth, and depth differs severalfold here.
Counts are therefore rarefied to the smallest sample total, averaged over
repeated draws, before q0 is computed. q1 and q2 are computed on proportions,
where the effect of depth is weaker but not absent.

Inputs:
  --bracken   directory of Bracken genus reports, one per sample
  --outdir    where to write (default: current directory)
  --n-rare    rarefaction draws (default 100)

Outputs:
  Fig_alpha_diversity.tif / .jpg
  alpha_diversity.csv
  alpha_diversity_tests.csv

Usage:
  python3 fig_alpha_diversity.py --bracken bracken/
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SAMPLES = [f"Pres0{i}" for i in range(1, 8)] + [f"Deg0{i}" for i in range(1, 8)]

COL_PRES = "#2c7bb6"
COL_DEG = "#d7191c"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "savefig.bbox": "tight",
})


def sample_of(path):
    name = os.path.basename(path)
    for s in SAMPLES:
        if s.lower() in name.lower():
            return s
    return None


def read_bracken(directory):
    counts = {}
    for path in sorted(glob.glob(os.path.join(directory, "*"))):
        if not os.path.isfile(path):
            continue
        s = sample_of(path)
        if s is None:
            continue
        d = pd.read_csv(path, sep="\t")
        if not {"name", "taxonomy_lvl", "new_est_reads"}.issubset(d.columns):
            print(f"[skipped] {os.path.basename(path)}: unexpected columns")
            continue
        d = d[d.taxonomy_lvl == "G"]
        counts[s] = d.set_index("name").new_est_reads
    if not counts:
        sys.exit(f"no Bracken reports recognised in {directory}")
    m = pd.DataFrame(counts).fillna(0)
    return m.reindex(columns=[s for s in SAMPLES if s in m.columns])


def rarefied_richness(counts, depth, n_draws, rng):
    """Mean genera observed in repeated draws of `depth` reads."""
    counts = counts[counts > 0].astype(int)
    total = counts.sum()
    if total <= depth:
        return float((counts > 0).sum())
    p = counts.values / total
    obs = [int((rng.multinomial(depth, p) > 0).sum()) for _ in range(n_draws)]
    return float(np.mean(obs))


def hill(p, q):
    p = p[p > 0]
    if q == 1:
        return float(np.exp(-(p * np.log(p)).sum()))
    return float((p ** q).sum() ** (1 / (1 - q)))


# ---------------------------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument("--bracken", required=True)
ap.add_argument("--outdir", default=".")
ap.add_argument("--n-rare", type=int, default=100)
args = ap.parse_args()

os.makedirs(args.outdir, exist_ok=True)
rng = np.random.default_rng(42)

m = read_bracken(args.bracken)
depth = int(m.sum().min())
print(f"{m.shape[0]} genera, {m.shape[1]} samples")
print(f"reads per sample: {int(m.sum().min())} to {int(m.sum().max())}")
print(f"rarefying to {depth} reads, {args.n_rare} draws\n")

rows = []
for s in m.columns:
    col = m[s]
    p = (col / col.sum()).values
    rows.append(dict(
        sample=s,
        status="preserved" if s.startswith("Pres") else "degraded",
        reads=int(col.sum()),
        q0_rarefied=round(rarefied_richness(col, depth, args.n_rare, rng), 1),
        q0_observed=int((col > 0).sum()),
        q1=round(hill(p, 1), 2),
        q2=round(hill(p, 2), 2),
        shannon=round(float(np.log(hill(p, 1))), 3)))

d = pd.DataFrame(rows)
d.to_csv(os.path.join(args.outdir, "alpha_diversity.csv"), index=False)
print(d.to_string(index=False))

# ---------------------------------------------------------------------
INDICES = [("q0_rarefied", "$^{0}D$ (rarefied)"),
           ("q1", "$^{1}D$"),
           ("q2", "$^{2}D$")]

print("\n--- preserved vs degraded ---")
tests = {}
for key in [k for k, _ in INDICES] + ["shannon"]:
    a = d.loc[d.status == "preserved", key]
    b = d.loc[d.status == "degraded", key]
    u, p = mannwhitneyu(a, b)
    tests[key] = dict(index=key, preserved_median=round(a.median(), 3),
                      degraded_median=round(b.median(), 3),
                      U=float(u), p=round(p, 4))
    note = "" if key in dict(INDICES) else "   (table only)"
    print(f"  {key:14s} preserved {a.median():9.2f}  degraded {b.median():9.2f}"
          f"   U = {u:5.1f}  p = {p:.4f}{note}")

pd.DataFrame(tests.values()).to_csv(
    os.path.join(args.outdir, "alpha_diversity_tests.csv"), index=False)

print("\nShannon is in the table but not plotted: H = ln(q1), so its test is"
      "\nidentical to that of q1 by construction.")

# ---------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(6.4, 3.0))

for panel_letter, ax, (key, ylabel) in zip("abc", axes, INDICES):
    groups = [d.loc[d.status == g, key].values
              for g in ("preserved", "degraded")]

    bp = ax.boxplot(groups, widths=0.55, patch_artist=True,
                    medianprops=dict(color="#1a1a1a", linewidth=1.2),
                    whiskerprops=dict(color="#595959", linewidth=0.8),
                    capprops=dict(color="#595959", linewidth=0.8),
                    flierprops=dict(marker="", linestyle="none"))
    for patch, colour in zip(bp["boxes"], (COL_PRES, COL_DEG)):
        patch.set_facecolor(colour)
        patch.set_alpha(0.28)
        patch.set_edgecolor(colour)
        patch.set_linewidth(1.0)

    # the points matter more than the box with seven per group
    for i, (vals, colour) in enumerate(zip(groups, (COL_PRES, COL_DEG)), start=1):
        x = rng.normal(i, 0.055, size=len(vals))
        ax.scatter(x, vals, s=24, c=colour, edgecolors="white",
                   linewidths=0.6, zorder=3)

    # the p value sits above a bracket spanning the two groups, placed in
    # axis coordinates so it clears the highest point whatever the scale
    pv = tests[key]["p"]
    lo, hi = ax.get_ylim()
    span = hi - lo
    ax.set_ylim(lo, hi + 0.16 * span)
    y = hi + 0.03 * span
    ax.plot([1, 1, 2, 2], [y, y + 0.02 * span, y + 0.02 * span, y],
            color="#595959", linewidth=0.8, clip_on=False)
    ax.text(1.5, y + 0.045 * span,
            "$p$ < 0.001" if pv < 0.001 else f"$p$ = {pv:.3f}",
            ha="center", va="bottom", fontsize=7.5, color="#1a1a1a")

    ax.text(-0.30, 1.04, f"({panel_letter})", transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="bottom", ha="left")

    ax.set_xticks([1, 2])
    ax.set_xticklabels(["P", "D"])
    ax.set_ylabel(ylabel)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

fig.tight_layout()
base = os.path.join(args.outdir, "Fig_alpha_diversity")
fig.savefig(f"{base}.tif", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(f"{base}.jpg", dpi=300, pil_kwargs={"quality": 95})
plt.close(fig)

print(f"\nWritten: {base}.tif, {base}.jpg")
print("         alpha_diversity.csv, alpha_diversity_tests.csv")
print("\nEach panel carries its Mann-Whitney p value: (a) q0, (b) q1, (c) q2.")
print("P is preserved and D is degraded; medians and U statistics are in the")
print("tests file.")
