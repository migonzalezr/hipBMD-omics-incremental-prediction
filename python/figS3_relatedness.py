# ======================================================================
# figS3_relatedness.py — Figure S3: relatedness screening (PLINK PI_HAT and KING)
# Inputs : derived/qc/pihat_histogram.tsv, derived/qc/related_pairs_king_pihat.tsv, derived/qc/pihat_sample_pc1.tsv
#          (summary files exported from the PLINK QC and GWAS runs; no individual identifiers are plotted)
# Outputs: figures/png/figS3_relatedness.png, figures/pdf/figS3_relatedness.pdf
# Run from the repository root.
# ======================================================================
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
BASE = os.environ.get("HIPBMD_ROOT", ".")

def apply_figure_style(sizes=(8,7,6)):
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': sizes[2],
        'axes.titlesize': sizes[1],
        'axes.labelsize': sizes[2],
        'xtick.labelsize': sizes[2],
        'ytick.labelsize': sizes[2],
    })

def panel_letter(ax, letter):
    ax.text(-0.12, 1.05, letter, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', ha='left')

h = pd.read_csv(os.path.join(BASE, "derived/qc/pihat_histogram.tsv"), sep="\t", header=None, names=["bin","n"]).sort_values("bin")
kp = pd.read_csv(os.path.join(BASE, "derived/qc/related_pairs_king_pihat.tsv"), sep="\t", header=None, names=["IID1","IID2","PI_HAT","Z0","Z1","Z2","KING"])
s = pd.read_csv(os.path.join(BASE, "derived/qc/pihat_sample_pc1.tsv"), sep="\t", header=None, names=["PI_HAT","pc1a","pc1b"])

kp["deg"] = np.where(kp.KING>=0.0884,"2nd","3rd")

g = lambda x: np.where(x>0.0, "EUR-like (PC1>0)", "AFR-like (PC1<0)")
s["pair"] = [" / ".join(sorted([a,b])) for a,b in zip(g(s.pc1a), g(s.pc1b))]

apply_figure_style(sizes=(8,7,6))
C_AB, C_OT, C2, C3, THR = "#C2562B", "#4C78A8", "#7B3294", "#E08214", "#555555"
fig, axs = plt.subplots(1,3, figsize=(7.2,2.7), gridspec_kw=dict(width_ratios=[1.1,1,1]))

ax=axs[0]; hh=h.copy()
ax.bar(hh.bin+0.0025, hh.n, width=0.005, color="#9a9a9a", edgecolor="none")
ax.set_yscale("log"); ax.set_xlim(-0.005,0.25); ax.set_ylim(30, 4e7)
ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v,_: {1e2:"100",1e3:"1k",1e4:"10k",1e5:"100k",1e6:"1M",1e7:"10M"}.get(v,"")))
ax.axvline(0.2, color=THR, ls="--", lw=0.8); ax.text(0.203, 1.2e7, "QC threshold\n(PI_HAT = 0.2)", fontsize=6, color=THR, va="top")
ax.set_xlabel("PI_HAT (PLINK IBD)"); ax.set_ylabel("Participant pairs")
ax.set_title("Spurious PI_HAT mode near 0.16", loc="left")

ax=axs[1]; bins=np.arange(0,0.205,0.005)
grp = np.where((s.pc1a>0)&(s.pc1b>0), "both PC1 > 0", "other pairs")
for lab,col in [("other pairs",C_OT),("both PC1 > 0",C_AB)]:
    v = s.PI_HAT[grp==lab]; w = np.full(len(v), 1/len(v))
    ax.hist(v, bins=bins, weights=w, histtype="step", lw=1.2, color=col)
ax.set_yscale("log"); ax.set_ylim(1e-5, 1.5)
ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v,_: {1e-4:"0.01%",1e-3:"0.1%",1e-2:"1%",1e-1:"10%",1:"100%"}.get(v,"")))
ax.axvline(0.2, color=THR, ls="--", lw=0.8)
ax.text(0.022, 2.5e-3, "PC1 < 0 or\nmixed pairs", color=C_OT, fontsize=7, va="bottom")
ax.text(0.150, 0.35, "both PC1 > 0", color=C_AB, fontsize=7, va="bottom", ha="center")
ax.set_xlabel("PI_HAT (PLINK IBD)"); ax.set_ylabel("Share of pairs in group")
ax.set_title("Confined to one PC1 cluster", loc="left")

ax=axs[2]
for d,col,lab in [("3rd",C3,"third degree"),("2nd",C2,"second degree")]:
    q=kp[kp.deg==d]; ax.scatter(q.KING, q.PI_HAT, s=9, color=col, edgecolor="white", lw=0.3, zorder=3, label=f"{lab} (n = {len(q)})")
ax.axhline(0.2, color=THR, ls="--", lw=0.8); ax.axvline(0.0442, color="#999999", lw=0.6, ls=":"); ax.axvline(0.0884, color="#999999", lw=0.6, ls=":")
ax.set_xlim(0.040, 0.108); ax.set_ylim(0.07, 0.215)
ax.set_xlabel("KING-robust kinship coefficient"); ax.set_ylabel("PI_HAT (PLINK IBD)")
ax.legend(loc="lower right", frameon=False, fontsize=6, handletextpad=0.2, borderaxespad=0.1, labelspacing=0.3)
ax.set_title("KING pairs all below 0.2", loc="left")

for a,l in zip(axs,"abc"): panel_letter(a, l)
fig.tight_layout(w_pad=1.2)
os.makedirs("suppfig", exist_ok=True)
fig.savefig(os.path.join(BASE, "figures/png/figS3_relatedness.png"), dpi=300)
fig.savefig(os.path.join(BASE, "figures/pdf/figS3_relatedness.pdf"))
