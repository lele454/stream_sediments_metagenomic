#!/usr/bin/env python3
"""
Figure: NMDS of taxonomic and functional composition, side by side
--------------------------------------------------------------------
Two panels on a shared layout: genus composition from the Bracken reports, and
KEGG orthologue composition from the KO table. Both on Bray-Curtis distances.

The two panels answer different questions and are worth showing together. A
taxonomic difference between the stream groups says the assemblages differ; a
functional difference says the gene content differs, which need not follow.
Where taxonomy separates and function does not, different organisms are doing
comparable things.

Each panel carries its PERMANOVA result and NMDS stress, and a dashed
ellipse at 1.96 standard deviations around each group. The ellipse helps the
eye group the points; it is not a test, and with seven samples per group it
should not be read as one. The within-group
dispersion is printed to the console instead: a significant PERMANOVA with
unequal dispersion can reflect heterogeneity rather than a shift in location,
so the two belong together in the caption even though only one fits on the
plot.

Abundances are converted to proportions within each sample before the
distances are computed. Sequencing depth differs severalfold between samples
and Bray-Curtis on raw counts would track depth.

Inputs:
  --bracken   directory of Bracken genus reports, one per sample
  --ko        KO abundance table (.xlsx or .csv), KEGG_KO in the first column
  --outdir    where to write (default: current directory)

Outputs:
  Fig_NMDS_taxonomic_functional.tif / .jpg
  nmds_scores.csv        sample scores of both ordinations
  permanova_results.csv

Usage:
  python3 fig_nmds.py --bracken bracken/ --ko Table_S4_KO_abundance.xlsx
"""

import os
import re
import sys
import glob
import argparse
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

try:
    from sklearn.manifold import MDS
except ImportError:
    sys.exit("scikit-learn is required for NMDS:\n"
             "  conda install -c conda-forge scikit-learn")

SAMPLES = [f"Pres0{i}" for i in range(1, 8)] + [f"Deg0{i}" for i in range(1, 8)]

COL_PRES = "#2c7bb6"
COL_DEG = "#d7191c"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "savefig.bbox": "tight",
})


# ---------------------------------------------------------------------
def sample_of(path):
    name = os.path.basename(path)
    for s in SAMPLES:
        if s.lower() in name.lower():
            return s
    return None


def read_bracken(directory):
    """Genus-by-sample counts from Bracken reports."""
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


def read_ko(path):
    d = (pd.read_excel(path) if path.lower().endswith((".xlsx", ".xls"))
         else pd.read_csv(path))
    d = d.set_index(d.columns[0])
    cols = [c for c in SAMPLES if c in d.columns]
    if not cols:
        sys.exit(f"no sample columns of the expected names found in {path}")
    return d[cols]


def permanova(D, groups, n_perm=9999, seed=42):
    """
    One-way PERMANOVA on a distance matrix. Implemented here rather than
    imported so the script has no dependency beyond scipy: the statistic is
    the ratio of between- to within-group sums of squares from the distances,
    and its null distribution comes from permuting the labels.
    """
    g = np.asarray(groups)
    n = len(g)
    levels = np.unique(g)
    a = len(levels)

    def ss_within(labels):
        w = 0.0
        for lv in levels:
            idx = np.where(labels == lv)[0]
            if len(idx) == 0:
                return np.nan
            sub = D[np.ix_(idx, idx)]
            w += (sub ** 2).sum() / (2 * len(idx))
        return w

    total = (D ** 2).sum() / (2 * n)

    def F_of(labels):
        w = ss_within(labels)
        return ((total - w) / (a - 1)) / (w / (n - a))

    obs = F_of(g)
    r2 = 1 - ss_within(g) / total

    rng = np.random.default_rng(seed)
    null = np.array([F_of(rng.permutation(g)) for _ in range(n_perm)])
    p = (np.sum(null >= obs) + 1) / (n_perm + 1)
    return dict(R2=r2, F=obs, p=p, n_perm=n_perm)


def dispersion(D, groups):
    """Mean within-group distance, to be read alongside the PERMANOVA."""
    g = np.asarray(groups)
    out = {}
    for lv in np.unique(g):
        idx = np.where(g == lv)[0]
        sub = D[np.ix_(idx, idx)]
        out[lv] = float(sub[np.triu_indices(len(idx), 1)].mean())
    return out


def confidence_ellipse(ax, x, y, colour, n_std=1.96):
    """
    Ellipse covering the group at the stated multiple of the standard
    deviation, drawn from the covariance of the two ordination axes. With
    seven points per group this is a visual aid to reading the panel, not an
    inference: the PERMANOVA is what tests the grouping.
    """
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    if np.any(vals <= 0):
        return
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    w, h = 2 * n_std * np.sqrt(vals)
    e = Ellipse((x.mean(), y.mean()), width=w, height=h, angle=angle,
                facecolor=colour, alpha=0.12, edgecolor=colour,
                linewidth=1.0, linestyle="--", zorder=2)
    ax.add_patch(e)


def nmds(D, seed=42):
    """Non-metric MDS on a precomputed distance matrix, two dimensions."""
    # scikit-learn is renaming these arguments; the call is built to suit the
    # installed version so the script does not emit deprecation warnings or
    # break on the rename
    import inspect
    sig = inspect.signature(MDS.__init__).parameters
    kw = dict(n_components=2, random_state=seed, n_init=50, max_iter=1000,
              normalized_stress=True)
    kw["metric_mds" if "metric_mds" in sig else "metric"] = False
    if "dissimilarity" in sig:
        kw["dissimilarity"] = "precomputed"
    else:
        kw["metric"] = "precomputed"
    if "init" in sig:
        # pin it: the default is scheduled to change from random to
        # classical_mds, which would silently alter the ordination
        kw["init"] = "random"
    mds = MDS(**kw)
    coords = mds.fit_transform(D)
    return coords, mds.stress_


# ---------------------------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument("--bracken", required=True)
ap.add_argument("--ko", required=True)
ap.add_argument("--outdir", default=".")
args = ap.parse_args()

os.makedirs(args.outdir, exist_ok=True)

tax = read_bracken(args.bracken)
fun = read_ko(args.ko)

shared = [s for s in SAMPLES if s in tax.columns and s in fun.columns]
if len(shared) < 6:
    sys.exit(f"only {len(shared)} samples present in both tables")
tax, fun = tax[shared], fun[shared]
tax = tax.loc[tax.sum(axis=1) > 0]
fun = fun.loc[fun.sum(axis=1) > 0]

status = np.array(["preserved" if s.startswith("Pres") else "degraded"
                   for s in shared])

print(f"{len(shared)} samples: {int((status=='preserved').sum())} preserved, "
      f"{int((status=='degraded').sum())} degraded")
print(f"taxonomy: {tax.shape[0]} genera")
print(f"function: {fun.shape[0]} KEGG orthologues\n")

panels = []
records = []
for label, table in [("Taxonomic (genus)", tax), ("Functional (KEGG KO)", fun)]:
    rel = table / table.sum()
    D = squareform(pdist(rel.T.values, metric="braycurtis"))
    res = permanova(D, status)
    disp = dispersion(D, status)
    coords, stress = nmds(D)

    print(f"--- {label} ---")
    print(f"  PERMANOVA  R2 = {res['R2']:.3f}  F = {res['F']:.2f}  "
          f"p = {res['p']:.4f}  ({res['n_perm']} permutations)")
    print(f"  dispersion preserved {disp['preserved']:.3f}, "
          f"degraded {disp['degraded']:.3f}")
    print(f"  NMDS stress {stress:.4f}\n")

    panels.append((label, coords, stress, res))
    records.append(dict(dataset=label, n_features=table.shape[0],
                        R2=round(res["R2"], 4), F=round(res["F"], 3),
                        p=res["p"], stress=round(float(stress), 4),
                        disp_preserved=round(disp["preserved"], 4),
                        disp_degraded=round(disp["degraded"], 4)))

pd.DataFrame(records).to_csv(
    os.path.join(args.outdir, "permanova_results.csv"), index=False)

if any(r["stress"] > 0.2 for r in records):
    print("[note] a stress above 0.2 means the two dimensions represent the")
    print("       distances poorly; read such an ordination with caution.\n")

# ---------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5))

scores = []
for panel_letter, ax, (label, coords, stress, res) in zip("ab", axes, panels):
    for grp, colour in (("preserved", COL_PRES), ("degraded", COL_DEG)):
        sel = status == grp
        confidence_ellipse(ax, coords[sel, 0], coords[sel, 1], colour)
        ax.scatter(coords[sel, 0], coords[sel, 1], s=46, c=colour,
                   edgecolors="white", linewidths=0.8, zorder=3)

    ax.axhline(0, color="#d9d9d9", linewidth=0.6, zorder=1)
    ax.axvline(0, color="#d9d9d9", linewidth=0.6, zorder=1)

    # p is reported as an inequality at the permutation floor: with 9999
    # permutations no smaller value is resolvable, and printing one would
    # imply precision the test does not have
    pstr = (f"$p$ < {1 / (res['n_perm'] + 1):.4f}"
            if res["p"] <= 1 / (res["n_perm"] + 1)
            else f"$p$ = {res['p']:.3f}")
    ax.text(0.97, 0.97,
            f"$R^2$ = {res['R2']:.3f}\n{pstr}\nstress = {stress:.3f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            linespacing=1.5, color="#1a1a1a")

    # panel letter outside the axes, top left, so it never sits on a point
    ax.text(-0.18, 1.04, f"({panel_letter})", transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="bottom", ha="left")

    ax.set_xlabel("NMDS1")
    ax.set_ylabel("NMDS2")
    ax.set_aspect("equal", adjustable="datalim")
    # headroom so the statistics block does not sit on top of a point
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.18 * (hi - lo))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    for s, (x, y) in zip(shared, coords):
        scores.append(dict(dataset=label, sample=s,
                           status="preserved" if s.startswith("Pres") else "degraded",
                           NMDS1=round(float(x), 4), NMDS2=round(float(y), 4)))

pd.DataFrame(scores).to_csv(
    os.path.join(args.outdir, "nmds_scores.csv"), index=False)

fig.tight_layout()
base = os.path.join(args.outdir, "Fig_NMDS_taxonomic_functional")
fig.savefig(f"{base}.tif", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(f"{base}.jpg", dpi=300, pil_kwargs={"quality": 95})
plt.close(fig)

print(f"Written: {base}.tif, {base}.jpg")
print("         nmds_scores.csv, permanova_results.csv")
print("\nEach panel carries its PERMANOVA and stress. (a) is taxonomic,")
print("(b) is functional, blue is preserved and red is degraded; those go")
print("in the caption, along with the dispersion values above.")
