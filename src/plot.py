import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import MultipleLocator

# ============================================================
# Setup
# ============================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE_DIR)

# ============================================================
# Helper functions
# ============================================================
def cdf(values):
    x = np.sort(values)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y


def beautify_axes(
    ax,
    x_major=None, x_minor=None,
    y_major=None, y_minor=None
):
    if x_major is not None:
        ax.xaxis.set_major_locator(MultipleLocator(x_major))
    if x_minor is not None:
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor))

    if y_major is not None:
        ax.yaxis.set_major_locator(MultipleLocator(y_major))
    if y_minor is not None:
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor))

    ax.grid(True, which="major", alpha=0.4)
    ax.grid(True, which="minor", alpha=0.15)


def logP(df):
    return -np.log10(df[df["pvalue"] > 0]["pvalue"])


# ============================================================
# Load data
# ============================================================
avg_real  = pd.read_csv("segments_real_avg.csv")
avg_null  = pd.read_csv("segments_null_avg.csv")
pv_real   = pd.read_csv("segments_real_pvalue.csv")
pv_null   = pd.read_csv("segments_null_pvalue.csv")

# keep valid p-values
avg_real = avg_real[avg_real["pvalue"] > 0]
avg_null = avg_null[avg_null["pvalue"] > 0]
pv_real  = pv_real[pv_real["pvalue"] > 0]
pv_null  = pv_null[pv_null["pvalue"] > 0]

# ============================================================
# 1. CDF — AVG pipeline (Real vs Null)
# ============================================================
x_r, y_r = cdf(avg_real["pvalue"])
x_n, y_n = cdf(avg_null["pvalue"])

plt.figure(figsize=(8,6))
plt.plot(x_r, y_r, linewidth=2, label="Real segments (AVG)")
plt.plot(x_n, y_n, "--", linewidth=2, label="Null segments (AVG)", color="black")

plt.xlabel("P-value")
plt.ylabel("CDF")
plt.title("CDF of P-values — AVG-based segmentation")

ax = plt.gca()
beautify_axes(ax, x_major=0.05, x_minor=0.01, y_major=0.1, y_minor=0.05)

plt.legend()
plt.tight_layout()
plt.show()

# ============================================================
# 2. CDF — P-value pipeline (Real vs Null)
# ============================================================
x_r, y_r = cdf(pv_real["pvalue"])
x_n, y_n = cdf(pv_null["pvalue"])

plt.figure(figsize=(8,6))
plt.plot(x_r, y_r, linewidth=2, label="Real segments (P-value)")
plt.plot(x_n, y_n, "--", linewidth=2, label="Null segments (P-value)", color="black")

plt.xlabel("P-value")
plt.ylabel("CDF")
plt.title("CDF of P-values — P-value–based segmentation")

ax = plt.gca()
beautify_axes(ax, x_major=0.05, x_minor=0.01, y_major=0.1, y_minor=0.05)

plt.legend()
plt.tight_layout()
plt.show()

# ============================================================
# 3. CDF — Comparison of segmentation strategies (Real only)
# ============================================================
x_avg, y_avg = cdf(avg_real["pvalue"])
x_pv,  y_pv  = cdf(pv_real["pvalue"])

plt.figure(figsize=(8,6))
plt.plot(x_avg, y_avg, linewidth=2, label="AVG-based segmentation")
plt.plot(x_pv,  y_pv,  linewidth=2, label="P-value–based segmentation")

plt.xlabel("P-value")
plt.ylabel("CDF")
plt.title("Comparison of segmentation strategies (Real segments)")

ax = plt.gca()
beautify_axes(ax, x_major=0.05, x_minor=0.01, y_major=0.1, y_minor=0.05)

plt.legend()
plt.tight_layout()
plt.show()

# ============================================================
# 4. Histogram of significance (-log10 P)
# ============================================================
avg_real_lp = logP(avg_real)
avg_null_lp = logP(avg_null)
pv_real_lp  = logP(pv_real)
pv_null_lp  = logP(pv_null)

bins = np.linspace(
    0,
    max(avg_real_lp.max(), pv_real_lp.max()),
    40
)

plt.figure(figsize=(12,5))

# AVG pipeline
plt.subplot(1,2,1)
plt.hist(avg_real_lp, bins=bins, alpha=0.8, label="Real", edgecolor="black")
plt.hist(avg_null_lp, bins=bins, alpha=0.6, label="Null", edgecolor="black")

plt.xlabel("-log10(P-value)")
plt.ylabel("Number of segments")
plt.title("AVG-based segmentation")

ax = plt.gca()
beautify_axes(ax, x_major=2, x_minor=1, y_major=5, y_minor=1)
plt.legend()

# P-value pipeline
plt.subplot(1,2,2)
plt.hist(pv_real_lp, bins=bins, alpha=0.8, label="Real", edgecolor="black")
plt.hist(pv_null_lp, bins=bins, alpha=0.6, label="Null", edgecolor="black")

plt.xlabel("-log10(P-value)")
plt.ylabel("Number of segments")
plt.title("P-value–based segmentation")

ax = plt.gca()
beautify_axes(ax, x_major=2, x_minor=1, y_major=5, y_minor=1)
plt.legend()

plt.tight_layout()
plt.show()

# ============================================================
# 5. Scatter — AVG vs Significance
# ============================================================
avg_real["logP"] = -np.log10(avg_real["pvalue"])
pv_real["logP"]  = -np.log10(pv_real["pvalue"])

plt.figure(figsize=(8,6))

plt.scatter(
    avg_real["AVG"], avg_real["logP"],
    s=35, alpha=0.7, label="AVG-based pipeline"
)

plt.scatter(
    pv_real["AVG"], pv_real["logP"],
    s=35, alpha=0.7, label="P-value–based pipeline"
)

plt.xlabel("AVG match (%)")
plt.ylabel("-log10(P-value)")
plt.title("Relationship between AVG match and statistical significance")

ax = plt.gca()
beautify_axes(ax, x_major=2, x_minor=1, y_major=2, y_minor=1)

plt.legend()
plt.tight_layout()
plt.show()

print("\nAll figures generated successfully.")
# ============================================================
# 6. Segment length vs statistical significance (STRONG FIGURE)
# ============================================================

plt.figure(figsize=(8,6))

# compute logP if not exists
avg_real["logP"] = -np.log10(avg_real["pvalue"])
pv_real["logP"]  = -np.log10(pv_real["pvalue"])

plt.scatter(
    avg_real["length"], avg_real["logP"],
    s=35, alpha=0.7,
    label="AVG-based pipeline"
)

plt.scatter(
    pv_real["length"], pv_real["logP"],
    s=35, alpha=0.7,
    label="P-value–based pipeline"
)

plt.xlabel("Segment length (bp)")
plt.ylabel("-log10(P-value)")
plt.title("Segment length vs statistical significance")

ax = plt.gca()

# Axes control (קטן, ברור, אקדמי)
beautify_axes(
    ax,
    x_major=500, x_minor=250,   # תוכל לשנות לפי L
    y_major=2,   y_minor=1
)

plt.legend()
plt.tight_layout()
plt.show()
# ============================================================
# FIXED CDF DIFFERENCE — meaningful range only
# ============================================================

p_grid = np.logspace(-6, -1, 400)  # 1e-6 → 0.1

def cdf_interp(values, grid):
    x, y = cdf(values)
    return np.interp(grid, x, y, left=0, right=1)

diff_avg = (
    cdf_interp(avg_real["pvalue"], p_grid)
    - cdf_interp(avg_null["pvalue"], p_grid)
)

diff_pv = (
    cdf_interp(pv_real["pvalue"], p_grid)
    - cdf_interp(pv_null["pvalue"], p_grid)
)

plt.figure(figsize=(8,6))

plt.plot(p_grid, diff_avg, label="AVG-based pipeline")
plt.plot(p_grid, diff_pv, label="P-value–based pipeline")

plt.axhline(0, color="black", linestyle="--", linewidth=1)

plt.xscale("log")
plt.xlabel("P-value (log scale)")
plt.ylabel("CDF(real) − CDF(null)")
plt.title("Separation between real and null segments (significant range)")

ax = plt.gca()
beautify_axes(ax, y_major=0.05, y_minor=0.01)

plt.legend()
plt.tight_layout()
plt.show()
