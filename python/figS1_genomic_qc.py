# ======================================================================
# figS1_genomic_qc.py — Figure S1: genomic quality control, population structure and GWAS diagnostics
# Inputs : derived/qc/variant_maf.tsv, derived/qc/heterozygosity_residuals.csv, derived/qc/pca_sample_meta.txt,
#          data/genotypes/los_qc_pca.eigenvec/.eigenval, results/gwas_full_cohort_pvalues.txt, results/gwas_qq_thinned.tsv
#          (summary files exported from the PLINK QC and GWAS runs; no individual identifiers are plotted)
# Outputs: figures/png/figS1_genomic_qc.png, figures/pdf/figS1_genomic_qc.pdf
# Run from the repository root.
# ======================================================================
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import chi2
import re
import os

BASE = os.environ.get("HIPBMD_ROOT", ".")
G = BASE + "/"
R9 = BASE + "/"

def apply_figure_style(sizes=(8, 7, 6)):
    mpl.rcParams.update({
        'font.family': 'Arial',
        'font.size': sizes[2],
        'axes.labelsize': sizes[1],
        'axes.titlesize': sizes[1],
        'xtick.labelsize': sizes[2],
        'ytick.labelsize': sizes[2],
        'axes.spines.top': False,
        'axes.spines.right': False,
        'figure.dpi': 100,
    })

def panel_letter(ax, letter):
    ax.text(-0.12, 1.05, letter, transform=ax.transAxes, fontsize=9,
            fontweight='bold', va='top', ha='left')

maf = pd.read_csv(BASE + "/derived/qc/variant_maf.tsv", sep="\t", header=None, names=["CHR", "MAF"], dtype={"CHR": np.int8, "MAF": np.float32})
pv = np.loadtxt(BASE + "/results/gwas_full_cohort_pvalues.txt", dtype=np.float64)

lam = np.median(chi2.isf(pv, 1)) / chi2.ppf(0.5, 1)
print(len(maf), len(pv), "lambda:", round(lam, 4), " min P:", pv.min(), " n<5e-8:", (pv < 5e-8).sum(), " n<=1e-5:", (pv <= 1e-5).sum())

hr = pd.read_csv(G + "derived/qc/heterozygosity_residuals.csv")
pcF = pd.read_csv(G + "data/genotypes/los_qc_pca.eigenvec", sep=r"\s+", header=None)
meta = pd.read_csv(G + "derived/qc/pca_sample_meta.txt", sep="\t")
meta.columns = [c.strip() for c in meta.columns]
meta["race"] = meta.race.astype(str).str.strip()

ev = pd.read_csv(G + "data/genotypes/los_qc_pca.eigenval", header=None)[0].values
evp = 100 * ev / ev.sum()

pcd = pcF.rename(columns={1: "IID", 2: "PC1", 3: "PC2", 4: "PC3"})[["IID", "PC1", "PC2", "PC3"]].merge(meta, on="IID")

qqt = pd.read_csv(R9 + "results/gwas_qq_thinned.tsv", sep="\t", header=None, names=["rank", "p"])
Nq = int(qqt.iloc[0, 1])
qqt = qqt.iloc[1:].astype(float)
qq_lambda = 1.043

L38 = {1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259, 6: 170805979,
       7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422, 11: 135086622, 12: 133275309,
       13: 114364328, 14: 107043718, 15: 101991189, 16: 90338345, 17: 83257441, 18: 80373285,
       19: 58617616, 20: 64444167, 21: 46709983, 22: 50818468}

cnt = maf.CHR.value_counts().sort_index()

RCOL = {"EUR": "#9E9E9E", "AFR": "#AA3377", "ASN": "#EE7733", "OTH": "#33BBEE"}
RLAB = {"EUR": "European", "AFR": "African American", "ASN": "Asian", "OTH": "Other"}

apply_figure_style(sizes=(8, 7, 6))

def build_S2(qq=None):
    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.05], hspace=0.62, wspace=0.42)
    for j, (pcy) in enumerate(["PC2", "PC3"]):
        ax = fig.add_subplot(gs[0, j])
        for r in ["EUR", "AFR", "ASN", "OTH"]:
            q = pcd[pcd.race == r]
            ax.scatter(q.PC1, q[pcy], s=2.5 if r in ("EUR", "AFR") else 7,
                       color=RCOL[r], alpha=0.55 if r in ("EUR", "AFR") else 0.95,
                       lw=0, label=f"{RLAB[r]} (n = {len(q):,})", rasterized=True)
        ax.set_xlabel("PC1"); ax.set_ylabel(pcy); ax.set_xticks([]); ax.set_yticks([]); ax.margins(0.04)
        if j == 0:
            panel_letter(ax, "a"); ax.set_title("Genetic structure (n = 6,732)", loc="left")
            ax.legend(loc="lower left", frameon=False, fontsize=6, markerscale=2.5,
                      handletextpad=0.1, borderaxespad=0.1, labelspacing=0.2)
    ax = fig.add_subplot(gs[0, 2]); k = np.arange(1, 11)
    ax.bar(k, evp, color="#9E9E9E", width=0.7); ax.set_yscale("log"); ax.set_ylim(0.5, 150)
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: {1: "1%", 10: "10%", 100: "100%"}.get(v, "")))
    ax.set_xticks(k); ax.set_xlabel("Principal component"); ax.set_ylabel("Share of top-10 eigenvalues")
    ax.text(1.5, evp[0] * 0.95, f"{evp[0]:.1f}%", fontsize=6, va="top")
    ax.text(2.45, evp[1] * 1.15, f"{evp[1]:.1f}%", fontsize=6)
    ax.set_title("PC1 dominates", loc="left"); panel_letter(ax, "b")
    ax = fig.add_subplot(gs[1, 0]); out = hr.het_outlier.astype(bool)
    ax.scatter(hr.PC1[~out], hr.F[~out], s=2, color="#9E9E9E", alpha=0.5, lw=0, rasterized=True)
    ax.scatter(hr.PC1[out], hr.F[out], s=6, color="#C0392B", lw=0, label=f"excluded (n = {out.sum()})")
    ax.set_xlabel("PC1 (before sample QC)"); ax.set_ylabel("Heterozygosity (F)"); ax.set_xticks([]); ax.margins(0.04)
    ax.legend(loc="lower right", frameon=False, fontsize=6, markerscale=1.8, handletextpad=0.1, borderaxespad=0.1)
    ax.set_title("F depends on ancestry (n = 7,678)", loc="left"); panel_letter(ax, "c")
    ax = fig.add_subplot(gs[1, 1]); sd_ = hr.het_resid.std(); mu = hr.het_resid.mean()
    ax.hist(hr.het_resid, bins=120, color="#9E9E9E"); ax.set_yscale("log")
    for s_ in (-3, 3): ax.axvline(mu + s_ * sd_, color="#C0392B", ls="--", lw=0.8)
    ax.text(mu + 3 * sd_, ax.get_ylim()[1] * 0.5, " ±3 SD", color="#C0392B", fontsize=6, va="top")
    ax.set_xlabel("F residual after PC1–PC10"); ax.set_ylabel("Participants")
    ax.set_title("Outliers flagged on residuals", loc="left")
    ax = fig.add_subplot(gs[1, 2]); panel_letter(ax, "e")
    if qq is not None:
        N, rk, pp = qq
        exp_ = -np.log10((rk - 0.5) / N); obs = -np.log10(pp)
        ax.scatter(exp_, obs, s=2, color="#333333", lw=0, rasterized=True)
        mx = max(exp_.max(), obs.max()) * 1.03; ax.plot([0, mx], [0, mx], color="#C0392B", lw=0.8)
        ax.axhline(-np.log10(5e-8), color="#777777", ls=":", lw=0.6)
        ax.set_xlim(0, exp_.max() * 1.04); ax.set_ylim(0, mx)
        ax.text(0.04, 0.96, f"λ = {qq_lambda:.3f}\nn = 6,704", transform=ax.transAxes, va="top", fontsize=7)
    else:
        ax.text(0.5, 0.5, "QQ pending\n(Cypress P values)", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="#C0392B")
    ax.set_xlabel("Expected −log₁₀ P"); ax.set_ylabel("Observed −log₁₀ P"); ax.set_title("Full-cohort GWAS of hip BMD", loc="left")
    ax = fig.add_subplot(gs[2, 0]); ax.hist(maf.MAF, bins=np.arange(0, 0.5025, 0.0125), color="#9E9E9E")
    ax.axvline(0.01, color="#C0392B", ls="--", lw=0.8)
    ax.text(0.02, ax.get_ylim()[1] * 0.93, "MAF cutoff 0.01", color="#C0392B", fontsize=6, va="top")
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v / 1e6:.1f}M" if v else "0"))
    ax.set_xlabel("Minor allele frequency"); ax.set_ylabel("Variants"); ax.set_title("7.30M variants after QC", loc="left"); panel_letter(ax, "d")
    ax = fig.add_subplot(gs[2, 1]); g_ = maf.groupby("CHR").MAF
    q1, md, q3 = g_.quantile(0.25), g_.median(), g_.quantile(0.75); ch = q1.index.values
    ax.vlines(ch, q1, q3, color="#9E9E9E", lw=2.2); ax.scatter(ch, md, s=8, color="#333333", zorder=3, marker="_")
    ax.set_xticks([1, 5, 10, 15, 20, 22]); ax.set_xlabel("Chromosome"); ax.set_ylabel("MAF (median, IQR)"); ax.set_ylim(0, None); ax.margins(x=0.03)
    ax.set_title("Similar across chromosomes", loc="left")
    ax = fig.add_subplot(gs[2, 2]); cnt_ = maf.CHR.value_counts().sort_index()
    ax.bar(cnt_.index, cnt_.values, color="#9E9E9E", width=0.75); ax.set_xticks([1, 5, 10, 15, 20, 22])
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v / 1e3:.0f}k" if v else "0"))
    ax.set_xlabel("Chromosome"); ax.set_ylabel("Variants"); ax.set_title("Variants per chromosome", loc="left")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.07)
    return fig

fS2 = build_S2(qq=(Nq, qqt["rank"].values, qqt["p"].values))
a1, a2 = fS2.axes[0], fS2.axes[1]
h_, l_ = a1.get_legend_handles_labels(); a1.get_legend().remove()
a2.legend(h_, l_, loc="upper left", frameon=False, fontsize=6, markerscale=2.5,
          handletextpad=0.1, borderaxespad=0.1, labelspacing=0.2)

for t in fS2.findobj(mpl.text.Text):
    if t.get_text().startswith("λ = "):
        t.set_position((0.97, 0.05)); t.set_ha("right"); t.set_va("bottom")

lets = [t for t in fS2.findobj(mpl.text.Text) if t.get_text() in ("d", "e") and t.get_fontweight() in ("bold", 700)]
for t in lets:
    t.set_text({"d": "e", "e": "d"}[t.get_text()])

os.makedirs("suppfig", exist_ok=True)
fS2.savefig("figures/png/figS1_genomic_qc.png", dpi=300)
fS2.savefig("figures/pdf/figS1_genomic_qc.pdf", dpi=300)
