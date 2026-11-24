# ======================================================================
# figS2_cohort_distributions.py — Figure S2
# #
# Distributions, proportions and correlations of clinical variables in the four analytic cohorts.
# #
# Inputs : clinical files
# Outputs: figures/png, figures/pdf
# Paper  : Figure S2
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
import os, sys
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

OLD = "--old-labels" in sys.argv
args = [a for a in sys.argv[1:] if not a.startswith("--")]
BASE = os.environ.get("HIPBMD_ROOT", os.path.abspath("."))
AR, R11 = BASE + "/derived/analysis_sets/", BASE + "/results/"
OUT = args[0] if args else BASE + "/figures"

mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 8, "axes.titlesize": 8,
    "axes.labelsize": 8, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7, "axes.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 3, "ytick.major.size": 3, "lines.linewidth": 1.2, "axes.titlelocation": "left",
    "axes.titlepad": 6, "axes.labelpad": 4, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.1})
def panel_letter(ax, letter, dx=-0.18, dy=1.02):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontweight="bold", fontsize=9, va="bottom", ha="left")

DISC, HOLD = "Modeling discovery", "Internal holdout"
LC, GZ = "LC subset", "Guangzhou cohort"
TITLE_A = "The LC subset and the Guangzhou cohort differ most in age, weight and height"

b = pd.read_excel(BASE + "/data/phenotypes/los_dxa_bmd.xlsx")
ph = pd.DataFrame({"IID": b.Sample_ID.astype(str).str.strip(), "y": pd.to_numeric(b.HTOT_BMD, errors="coerce")}).dropna().drop_duplicates("IID")
clB = pd.read_csv(BASE + "/data/phenotypes/los_clinical_imputed.csv"); clB["IID"] = clB.Subject_ID.astype(str).str.strip()
clL = pd.read_csv(BASE + "/data/phenotypes/los_lc_clinical_imputed.csv"); clL["IID"] = clL.Subject_ID.astype(str).str.strip()
chC = pd.read_csv(BASE + "/data/phenotypes/guangzhou_clinical.csv"); chC["IID"] = chC.ID.astype(str).str.strip()
pcs = pd.read_csv(BASE + "/data/genotypes/los_qc_pca.eigenvec", sep=r"\s+", header=None)
pcmap = pcs.rename(columns={1: "IID", 2: "gPC1", 3: "gPC2"})[["IID", "gPC1", "gPC2"]]; pcmap["IID"] = pcmap.IID.astype(str)
ids = {DISC: [l.strip() for l in open(AR + "ids_modeling_discovery.txt") if l.strip()],
       HOLD: pd.read_csv(R11 + "holdout_predictions.csv").IID.astype(str).tolist(),
       LC:   pd.read_csv(R11 + "lc_predictions.csv").IID.astype(str).tolist(),
       GZ:   pd.read_csv(R11 + "guangzhou_predictions.csv").IID.astype(str).tolist()}
CATS = {"Female": None, "Regular exercise": "Exercise_Regular_cat", "Smoking": "Smoke_cat", "Sun exposure": "Sun_Exposure_cat",
        "Alcohol": "Alcohol_Drinking_cat", "Cheese": "Cheese_Eat_cat", "Milk": "Milk_Drink_cat", "Yogurt": "Yogurt_Eat_cat",
        "Caucasian": "Race_Caucasian", "African American": "Race_African", "Asian or other": "Race_AsianOther"}
RACE3 = "Asian or other"
if OLD:
    CATS = {(k if k != RACE3 else "Asian"): (v if k != RACE3 else "Race_Asian") for k, v in CATS.items()}; RACE3 = "Asian"
for _cl in (clB, clL):
    if {"Race_Asian", "Race_Other"} <= set(_cl.columns):
        _a = _cl.Race_Asian.astype(str).str.strip(); _o = _cl.Race_Other.astype(str).str.strip(); _si = ["Yes", "1", "1.0", "True", "TRUE"]
        _cl["Race_AsianOther"] = np.where(_cl.Race_Asian.isna() & _cl.Race_Other.isna(), None, np.where(_a.isin(_si) | _o.isin(_si), "Yes", "No"))
def yes(v):
    s = v.astype(str).str.strip(); out = s.isin(["Yes", "yes", "1", "1.0", "True", "TRUE"]).astype(float); out[v.isna().values] = np.nan; return out
rows = []
for coh, idl in ids.items():
    d = pd.DataFrame({"IID": [str(i).strip() for i in idl]})
    if coh != GZ:
        cl = clB if coh != LC else clL
        d = d.merge(ph, on="IID").merge(cl, on="IID", how="left").merge(pcmap, on="IID", how="left")
        d["Female"] = (d.gender_dummy.astype(str).str.lower() == "female").astype(float)
        for k, c in CATS.items():
            if c: d[k] = yes(d[c]) if c in d.columns else np.nan
    else:
        d = d.merge(chC, on="IID", how="left"); d["y"] = pd.to_numeric(d.HTOT_BMD); d["Female"] = 1.0; d["gPC1"] = np.nan; d["gPC2"] = np.nan
        m = {"Regular exercise": "Exercise_cat", "Smoking": "Smoke_cat", "Alcohol": "Alcohol_cat", "Milk": "Milk_cat", "Yogurt": "Yogurt_cat"}
        for k in CATS:
            if k in m: d[k] = yes(d[m[k]])
            elif k in ("Caucasian", "African American", RACE3, "Sun exposure", "Cheese"): d[k] = np.nan
    d["Cohort"] = coh; rows.append(d[["Cohort", "IID", "y", "Age", "Weight_kg", "Height_cm", "BMI", "gPC1", "gPC2"] + list(CATS)])
S1 = pd.concat(rows, ignore_index=True)
assert S1.groupby("Cohort").size().to_dict() == {DISC: 1430, HOLD: 351, LC: 495, GZ: 492}

COH = list(ids); COL = {DISC: "#2b6cb0", HOLD: "#8fb3e0", LC: "#d98c00", GZ: "#117a65"}
NN = S1.groupby("Cohort").size().to_dict()
fig = plt.figure(figsize=(7.2, 7.4))
gs = fig.add_gridspec(3, 1, height_ratios=[1, 1.05, 1.3], hspace=0.42)
g1 = gs[0].subgridspec(1, 5, wspace=0.35); g3 = gs[2].subgridspec(1, 4, wspace=0.12)
CONT = [("y", "Total hip BMD (g/cm²)"), ("Age", "Age (years)"), ("Weight_kg", "Weight (kg)"), ("Height_cm", "Height (cm)"), ("BMI", "BMI (kg/m²)")]
A = []
for j, (v, lab) in enumerate(CONT):
    ax = fig.add_subplot(g1[0, j]); A.append(ax)
    allv = S1[v].dropna(); xs = np.linspace(allv.quantile(0.002), allv.quantile(0.998), 300)
    for c in COH:
        x = S1.loc[S1.Cohort == c, v].dropna()
        ax.plot(xs, gaussian_kde(x)(xs), color=COL[c], lw=1.3 if c != HOLD else 1.1, ls="-" if c != HOLD else (0, (3, 1.5)))
    ax.set_xlabel(lab); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    if j == 0: ax.set_ylabel("Density")
    ax.margins(x=0.02)
A[0].set_title(TITLE_A, loc="left", x=0)
hand = [Line2D([0], [0], color=COL[c], lw=1.4, ls="-" if c != HOLD else (0, (3, 1.5))) for c in COH]
axb = fig.add_subplot(gs[1]); cats = list(CATS); off = np.linspace(-0.27, 0.27, 4)
for k, c in enumerate(COH):
    p = S1[S1.Cohort == c][cats].mean() * 100
    xs = np.arange(len(cats)) + off[k]; ok = p.notna().values
    axb.vlines(xs[ok], 0, p.values[ok], color=COL[c], lw=1.0, alpha=0.7)
    axb.scatter(xs[ok], p.values[ok], s=16, color=COL[c], zorder=3, edgecolor="white", lw=0.3)
    for x in xs[~ok]: axb.text(x, 2, "n.c.", rotation=90, fontsize=6, color=COL[c], ha="center", va="bottom")
axb.set_xticks(range(len(cats)))
axb.set_xticklabels([c.replace(" ", "\n") if c in ("Regular exercise", "African American", "Sun exposure") else (c.replace(" or ", "\nor ") if c == "Asian or other" else c) for c in cats], fontsize=6)
axb.set_ylim(0, 105); axb.set_ylabel("Participants answering yes (%)"); axb.set_xlim(-0.6, len(cats) - 0.4)
for x in np.arange(len(cats) - 1) + 0.5: axb.axvline(x, color="#e5e5e5", lw=0.5, zorder=0)
axb.set_title("Sex, smoking, alcohol and dairy intake differ most between cohorts", loc="left")
VARS = [("y", "BMD"), ("Age", "Age"), ("Weight_kg", "Weight"), ("Height_cm", "Height"), ("BMI", "BMI"), ("gPC1", "gPC1"), ("gPC2", "gPC2")]
cmap = plt.get_cmap("RdBu_r"); AX3 = []
for k, c in enumerate(COH):
    ax = fig.add_subplot(g3[0, k]); AX3.append(ax)
    M_ = S1[S1.Cohort == c][[v for v, _ in VARS]].corr().values; n = len(VARS)
    for i in range(n):
        for j in range(n):
            if j > i: continue
            val = M_[i, j]
            if np.isnan(val):
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="#f2f2f2", edgecolor="white", lw=0.6)); continue
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=cmap((val + 1) / 2), edgecolor="white", lw=0.6))
            if i != j: ax.text(j, i, f"{val:.2f}".replace("-", "−").replace("0.", "."), ha="center", va="center", fontsize=5.5,
                               color="white" if abs(val) > 0.55 else "#222222")
    ax.set_xlim(-0.5, n - 0.5); ax.set_ylim(n - 0.5, -0.5); ax.set_aspect("equal")
    ax.set_xticks(range(n)); ax.set_xticklabels([l for _, l in VARS], rotation=90, fontsize=6)
    ax.set_yticks(range(n)); ax.set_yticklabels([l for _, l in VARS] if k == 0 else [], fontsize=6)
    for s_ in ax.spines.values(): s_.set_visible(False)
    ax.tick_params(length=0)
    lab = c
    ax.text(n - 0.6, -0.2, lab, ha="right", va="top", fontsize=7, color=COL[c], fontweight="bold")
AX3[0].set_title("Pearson correlations within each cohort", loc="left")
sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-1, 1), cmap=cmap)
cax = fig.add_axes([0.93, 0.085, 0.012, 0.17]); cb = fig.colorbar(sm, cax=cax, ticks=[-1, 0, 1]); cb.ax.tick_params(labelsize=6); cb.outline.set_visible(False)
cb.set_label("r", fontsize=7)
panel_letter(A[0], "a"); panel_letter(axb, "b"); panel_letter(AX3[0], "c")
fig.subplots_adjust(left=0.08, right=0.91, top=0.925, bottom=0.075)
fig.legend(hand, [f"{c} (n = {NN[c]:,})" for c in COH], loc="upper center", ncol=4, frameon=False,
           bbox_to_anchor=(0.5, 0.985), bbox_transform=fig.transFigure, fontsize=7, handlelength=2.2)
os.makedirs(OUT + "/png", exist_ok=True); os.makedirs(OUT + "/pdf", exist_ok=True)
fig.savefig(OUT + "/png/figS4_cohort_distributions.png", dpi=300); fig.savefig(OUT + "/pdf/figS4_cohort_distributions.pdf")
print("saved", OUT)
