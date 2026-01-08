import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import MultipleLocator
from scipy.stats import ks_2samp

# ============================================================
# Setup
# ============================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE_DIR)

OUT_DIR = "figures"
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# Helper functions
# ============================================================
def cdf(values):
    x = np.sort(values)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y

def beautify(ax, xM=None, xm=None, yM=None, ym=None):
    if xM: ax.xaxis.set_major_locator(MultipleLocator(xM))
    if xm: ax.xaxis.set_minor_locator(MultipleLocator(xm))
    if yM: ax.yaxis.set_major_locator(MultipleLocator(yM))
    if ym: ax.yaxis.set_minor_locator(MultipleLocator(ym))
    ax.grid(True, which="major", alpha=0.4)
    ax.grid(True, which="minor", alpha=0.15)

def logP(df):
    return -np.log10(df[df["pvalue"] > 0]["pvalue"])

# ============================================================
# Load data
# ============================================================
avg_real = pd.read_csv("segments_real_avg.csv")
avg_null = pd.read_csv("segments_null_avg.csv")
pv_real  = pd.read_csv("segments_real_pvalue.csv")
pv_null  = pd.read_csv("segments_null_pvalue.csv")

avg_real = avg_real[avg_real["pvalue"] > 0].copy()
avg_null = avg_null[avg_null["pvalue"] > 0].copy()
pv_real  = pv_real[pv_real["pvalue"] > 0].copy()
pv_null  = pv_null[pv_null["pvalue"] > 0].copy()

avg_real["logP"] = -np.log10(avg_real["pvalue"])
pv_real["logP"]  = -np.log10(pv_real["pvalue"])

classes = ["strong", "weak", "noise"]
colors  = {"strong":"red", "weak":"orange", "noise":"blue"}

# ============================================================
# 1. CDF — AVG (Real vs Null)
# ============================================================
plt.figure(figsize=(8,6))
x,y = cdf(avg_real["pvalue"])
xn,yn = cdf(avg_null["pvalue"])
plt.plot(x,y,label="Real")
plt.plot(xn,yn,"--",color="black",label="Null")
plt.xlabel("P-value"); plt.ylabel("CDF")
plt.title("CDF of P-values — AVG-based")
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig01_cdf_avg.pdf", dpi=300)
plt.savefig(f"{OUT_DIR}/fig01_cdf_avg.png", dpi=300)

# ============================================================
# 2. CDF — P-value (Real vs Null)
# ============================================================
plt.figure(figsize=(8,6))
x,y = cdf(pv_real["pvalue"])
xn,yn = cdf(pv_null["pvalue"])
plt.plot(x,y,label="Real")
plt.plot(xn,yn,"--",color="black",label="Null")
plt.xlabel("P-value"); plt.ylabel("CDF")
plt.title("CDF of P-values — P-value-based")
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig02_cdf_pvalue.pdf", dpi=300)
plt.savefig(f"{OUT_DIR}/fig02_cdf_pvalue.png", dpi=300)

# ============================================================
# 3. CDF — AVG vs P-value (Real only)
# ============================================================
plt.figure(figsize=(8,6))
x1,y1 = cdf(avg_real["pvalue"])
x2,y2 = cdf(pv_real["pvalue"])
plt.plot(x1,y1,label="AVG-based")
plt.plot(x2,y2,label="P-value-based")
plt.xlabel("P-value"); plt.ylabel("CDF")
plt.title("CDF comparison (REAL segments)")
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig03_cdf_compare.pdf", dpi=300)
plt.savefig(f"{OUT_DIR}/fig03_cdf_compare.png", dpi=300)

# ============================================================
# 4. Histogram of significance
# ============================================================
plt.figure(figsize=(9,5))
plt.hist(avg_real["logP"], bins=40, alpha=0.7, label="AVG")
plt.hist(pv_real["logP"],  bins=40, alpha=0.7, label="P-value")
plt.xlabel("-log10(P-value)")
plt.ylabel("Number of segments")
plt.title("Distribution of statistical significance")
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig04_hist_logp.pdf", dpi=300)
plt.savefig(f"{OUT_DIR}/fig04_hist_logp.png", dpi=300)

# ============================================================
# 5. CDF separation (Real − Null)
# ============================================================
p_grid = np.logspace(-6,-1,400)

def interp_cdf(v):
    x,y = cdf(v)
    return np.interp(p_grid,x,y,left=0,right=1)

plt.figure(figsize=(8,6))
plt.plot(p_grid, interp_cdf(avg_real["pvalue"]) - interp_cdf(avg_null["pvalue"]), label="AVG")
plt.plot(p_grid, interp_cdf(pv_real["pvalue"])  - interp_cdf(pv_null["pvalue"]),  label="P-value")
plt.axhline(0,color="black",ls="--")
plt.xscale("log")
plt.xlabel("P-value")
plt.ylabel("CDF(real) − CDF(null)")
plt.title("Separation between real and null")
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig05_cdf_separation.pdf", dpi=300)
plt.savefig(f"{OUT_DIR}/fig05_cdf_separation.png", dpi=300)

# ============================================================
# 6. KS test + bar plot
# ============================================================
ks_avg = ks_2samp(avg_real["pvalue"], avg_null["pvalue"])
ks_pv  = ks_2samp(pv_real["pvalue"],  pv_null["pvalue"])

plt.figure(figsize=(6,4))
plt.bar(["AVG","P-value"], [ks_avg.statistic, ks_pv.statistic], edgecolor="black")
plt.ylabel("KS distance")
plt.title("KS test: Real vs Null")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig06_ks.pdf", dpi=300)
plt.savefig(f"{OUT_DIR}/fig06_ks.png", dpi=300)

# ============================================================
# 7. Class distribution (counts)
# ============================================================
cnt_avg = avg_real["class"].value_counts().reindex(classes, fill_value=0)
cnt_pv  = pv_real["class"].value_counts().reindex(classes, fill_value=0)

x = np.arange(len(classes))
w = 0.35

plt.figure(figsize=(7,5))
plt.bar(x-w/2, cnt_avg, w, label="AVG")
plt.bar(x+w/2, cnt_pv,  w, label="P-value")
plt.xticks(x,[c.upper() for c in classes])
plt.ylabel("Number of segments")
plt.title("Class distribution (REAL)")
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig07_class_counts.pdf", dpi=300)
plt.savefig(f"{OUT_DIR}/fig07_class_counts.png", dpi=300)

# ============================================================
# 8. AVG per class
# ============================================================
plt.figure(figsize=(8,5))
data = [avg_real[avg_real["class"]==c]["AVG"] for c in classes]
plt.boxplot(data, labels=[c.upper() for c in classes], showfliers=False)
plt.ylabel("AVG match (%)")
plt.title("AVG distribution per class (REAL)")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig08_avg_per_class.pdf", dpi=300)
plt.savefig(f"{OUT_DIR}/fig08_avg_per_class.png", dpi=300)

# ============================================================
# 9. Significance per class
# ============================================================
plt.figure(figsize=(8,5))
data = [avg_real[avg_real["class"]==c]["logP"] for c in classes]
plt.boxplot(data, labels=[c.upper() for c in classes], showfliers=False)
plt.ylabel("-log10(P-value)")
plt.title("Statistical significance per class (REAL)")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig09_logp_per_class.pdf", dpi=300)
plt.savefig(f"{OUT_DIR}/fig09_logp_per_class.png", dpi=300)

# ============================================================
# 10. AVG vs significance (colored)
# ============================================================
plt.figure(figsize=(8,6))
for c in classes:
    d = avg_real[avg_real["class"]==c]
    plt.scatter(d["AVG"], d["logP"], label=c.upper(), color=colors[c], alpha=0.6)

plt.xlabel("AVG match (%)")
plt.ylabel("-log10(P-value)")
plt.title("AVG vs statistical significance (REAL)")
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig10_avg_vs_logp.pdf", dpi=300)
plt.savefig(f"{OUT_DIR}/fig10_avg_vs_logp.png", dpi=300)

# ============================================================
# Show all figures together
# ============================================================
plt.show()