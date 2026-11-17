"""Targeted lookup of Bacteroides vulgatus and BMD (Lin et al., Nat Commun 2023;14:6853).

Guangzhou cohort: lumbar spine L1-L4, total hip and femoral neck BMD, with the covariates of the original
report (age, BMI, years since menopause, FSH, exercise, family income) and with age + weight.
LOS women (BGI subset = modeling discovery + internal holdout): lumbar spine L1-L4 (TOT_BMD), total hip and
femoral neck BMD, with age + BMI + exercise and with age + weight. Lumbar spine BMD is used in LOS only here.

Exposure: CLR-transformed relative abundance of B. vulgatus (zeros = half the global minimum non-zero value).
Models: (1) B. vulgatus alone, OLS with covariates; (2) joint zero-sum log-ratio model over the species whose
mean relative abundance is at or above a threshold (approximates the constrained regression of the original).

Inputs : data/phenotypes/{los_dxa_bmd.xlsx, los_clinical_imputed.csv, los_questionnaire.xlsx,
         guangzhou_clinical.csv, guangzhou_dxa.xlsx, guangzhou_questionnaire.csv};
         data/microbiome/{los_bgi_species_abundance.csv, metagenome_id_map.xlsx, guangzhou_species_abundance.csv,
         metagenome_id_corrections.csv (optional)};
         derived/analysis_sets/{ids_modeling_discovery.txt, ids_internal_holdout.txt};
         results/species_association.csv (species universe, from R/14_diversity_and_species.R)
Output : results/bvulgatus_lookup.csv
"""
import os, re
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

BASE = os.environ.get("HIPBMD_ROOT", os.path.abspath("."))
P = lambda *a: os.path.join(BASE, *a)
BV = "s__Bacteroides_vulgatus"


def yes(v):
    s = v.astype(str).str.strip()
    o = s.isin(["Yes", "yes", "1", "1.0", "True", "TRUE"]).astype(float)
    o[v.isna().values] = np.nan
    return o


def clr_mat(ra):
    x = ra.astype(float).copy()
    x = x.div(x.sum(axis=1), axis=0)
    mn = x[x > 0].min().min()
    lx = np.log(x.where(x > 0, mn / 2))
    return lx.sub(lx.mean(axis=1), axis=0)


res = []


def fit_single(d, y, cov, cohort, sub, covlab):
    dd = d[[y, "bv"] + cov].dropna()
    if len(dd) < 30:
        return
    m = smf.ols(f"{y} ~ bv + " + " + ".join(cov), data=dd).fit()
    ci = m.conf_int().loc["bv"]
    res.append(dict(cohort=cohort, subgroup=sub, outcome=y, model="single species (CLR)", covariates=covlab,
                    threshold=np.nan, n_species=np.nan, n=len(dd), beta=m.params["bv"], lo=ci[0], hi=ci[1],
                    p=m.pvalues["bv"]))


def fit_joint(counts, ids, d, thr, y, cov, cohort, sub, covlab):
    ra = counts.loc[ids].div(counts.loc[ids].sum(axis=1), axis=0)
    sp = ra.columns[ra.mean() >= thr]
    if BV not in sp:
        return
    x = counts.loc[ids, sp].astype(float)
    x = x.div(x.sum(axis=1), axis=0)
    mn = x[x > 0].min().min()
    x = np.log(x.where(x > 0, mn / 2))
    ref = sp[sp != BV][0]                      # with the zero-sum constraint the reference is arbitrary
    lr = x.drop(columns=[ref]).sub(x[ref], axis=0)
    names = list(lr.columns)
    lr.columns = [f"x{i}" for i in range(lr.shape[1])]
    bvc = lr.columns[names.index(BV)]
    dd = d.set_index("IID").loc[ids, [y] + cov].join(lr).dropna()
    m = sm.OLS(dd[y], sm.add_constant(dd[cov + list(lr.columns)])).fit()
    res.append(dict(cohort=cohort, subgroup=sub, outcome=y, model="joint zero-sum log-ratio model", covariates=covlab,
                    threshold=thr, n_species=len(sp), n=len(dd), beta=m.params[bvc], lo=m.conf_int().loc[bvc, 0],
                    hi=m.conf_int().loc[bvc, 1], p=m.pvalues[bvc]))


# ---------------------------------------------------------------- Guangzhou cohort
tab = pd.read_csv(P("data/microbiome/guangzhou_species_abundance.csv"))
MG = tab.set_index(tab.columns[0]).T
MG.columns = MG.columns.astype(str)
MG = MG.loc[:, MG.columns.str.startswith("s__")]
MG = MG.T.groupby(level=0).sum().T
keep = (MG > 0).mean()
keep = keep[keep >= 0.01].index                               # prevalence >= 1%
clin_g = pd.read_csv(P("data/phenotypes/guangzhou_clinical.csv"))
clin_g["IID"] = clin_g.ID.astype(str).str.strip()
spine_g = pd.read_excel(P("data/phenotypes/guangzhou_dxa.xlsx"), "Spine")
hip_g = pd.read_excel(P("data/phenotypes/guangzhou_dxa.xlsx"), "Hip")
quest_g = pd.read_csv(P("data/phenotypes/guangzhou_questionnaire.csv"), encoding="latin-1")
income = {"Less than 6,000": 1, "Between 6,000 and 12,000": 2, "Between 12,000 and 24,000": 3,
          "Between 24,000 and 36,000": 4, "Between 36,000 and 60,000": 5, "Between 60,000 and 120,000": 6,
          "Above 120,000": 7}
bv_g = clr_mat(MG[keep])[BV]
dG = (pd.DataFrame({"IID": bv_g.index, "bv": bv_g.values})
      .merge(spine_g[["Subject_ID", "L1_L4_BMD"]].rename(columns={"Subject_ID": "IID"}), on="IID")
      .merge(hip_g[["Subject_ID", "HTOT_BMD", "NECK_BMD"]].rename(columns={"Subject_ID": "IID"}), on="IID")
      .merge(quest_g[["ID", "LengthSinceMenopause", "FSH", "FamilyAnnualIncome"]].rename(columns={"ID": "IID"}), on="IID")
      .merge(clin_g[["IID", "Age", "BMI", "Weight_kg", "Exercise_cat"]], on="IID"))
dG["income"] = dG.FamilyAnnualIncome.map(income)
dG["YSM"] = dG.LengthSinceMenopause
ra_g = MG[keep].div(MG[keep].sum(axis=1), axis=0)
nonrare = ra_g.columns[ra_g.mean() >= 1e-4]                    # sensitivity: CLR over non-rare species
dG["bv_nonrare"] = dG.IID.map(clr_mat(MG[nonrare])[BV])
COV_LIN = ["Age", "BMI", "YSM", "FSH", "Exercise_cat", "income"]
for y in ["L1_L4_BMD", "HTOT_BMD", "NECK_BMD"]:
    fit_single(dG, y, COV_LIN, "Guangzhou", "all", "Lin 2023 set")
    fit_single(dG, y, ["Age", "Weight_kg"], "Guangzhou", "all", "age + weight")
    fit_single(dG.assign(bv=dG.bv_nonrare), y, COV_LIN, "Guangzhou", "all (CLR over non-rare species)", "Lin 2023 set")
ids_g = list(dG.IID)
for thr in [5e-3, 2e-3, 1e-3, 5e-4]:
    for y in ["L1_L4_BMD", "HTOT_BMD", "NECK_BMD"]:
        fit_joint(MG[keep], ids_g, dG, thr, y, COV_LIN, "Guangzhou", "all", "Lin 2023 set")

# ---------------------------------------------------------------- LOS women, BGI subset
ids = [l.strip() for f in ["ids_modeling_discovery.txt", "ids_internal_holdout.txt"]
       for l in open(P("derived/analysis_sets", f)) if l.strip()]
bmd = pd.read_excel(P("data/phenotypes/los_dxa_bmd.xlsx"))
bmd = pd.DataFrame({"IID": bmd.Sample_ID.astype(str).str.strip(),
                    **{c: pd.to_numeric(bmd[c], errors="coerce") for c in ["HTOT_BMD", "NECK_BMD", "TOT_BMD"]}}
                   ).drop_duplicates("IID")
clin = pd.read_csv(P("data/phenotypes/los_clinical_imputed.csv"))
clin["IID"] = clin.Subject_ID.astype(str).str.strip()
d = pd.DataFrame({"IID": ids}).merge(bmd[["IID", "HTOT_BMD"]].dropna(), on="IID").merge(clin, on="IID", how="left")
d = d[d.gender_dummy.astype(str).str.lower() == "female"].copy()
d["Caucasian"] = yes(d.Race_Caucasian)
quest = pd.read_excel(P("data/phenotypes/los_questionnaire.xlsx"), sheet_name=0)
idc = [c for c in quest.columns if re.search("subject|sample|^id", str(c), re.I)][0]
quest["IID"] = quest[idc].astype(str).str.strip()
d = d[["IID", "Age", "Weight_kg", "BMI", "Caucasian", "Exercise_Regular_cat"]].merge(
    quest[["IID", "Postmenopausal"]].drop_duplicates("IID"), on="IID", how="left")

mgs = pd.read_csv(P("data/microbiome/los_bgi_species_abundance.csv"))
mgs = mgs.rename(columns={mgs.columns[0]: "Metagenome_ID"})
mgs["Metagenome_ID"] = mgs.Metagenome_ID.astype(str).str.strip()
fx = P("data/microbiome/metagenome_id_corrections.csv")       # columns: old, new
if os.path.exists(fx):
    c_ = pd.read_csv(fx, dtype=str)
    mgs["Metagenome_ID"] = mgs.Metagenome_ID.replace(dict(zip(c_.old.str.strip(), c_.new.str.strip())))
idmap = pd.read_excel(P("data/microbiome/metagenome_id_map.xlsx"), sheet_name="bgi").rename(columns={"BGI_1794": "IID"})
idmap = idmap[["Metagenome_ID", "IID"]].astype(str).apply(lambda s: s.str.strip())
ML = idmap[idmap.IID.isin(set(ids))].merge(mgs, on="Metagenome_ID").drop(columns=["Metagenome_ID"]).set_index("IID")
ML.columns = [re.sub(r".*(s__[^.]+).*", r"\1", c) if "s__" in c else None for c in ML.columns]
ML = ML.loc[:, [c is not None for c in ML.columns]]
ML = ML.T.groupby(level=0).sum().T
universe = pd.read_csv(P("results/species_association.csv")).especie.tolist()
bv_l = clr_mat(ML[universe])[BV]
dL = d.merge(pd.DataFrame({"IID": bv_l.index, "bv": bv_l.values}), on="IID").merge(
    bmd[["IID", "TOT_BMD", "NECK_BMD"]], on="IID", how="left").merge(bmd[["IID", "HTOT_BMD"]], on="IID", how="left")
dL["Exercise_cat"] = yes(dL.Exercise_Regular_cat)
SUBS = {"women 45–60 y": (dL.Age >= 45) & (dL.Age <= 60),
        "postmenopausal, 45–60 y": (dL.Age >= 45) & (dL.Age <= 60) & (dL.Postmenopausal == "Yes"),
        "white postmenopausal, ≥60 y (Lin 2023 LOS definition)": (dL.Age >= 60) & (dL.Postmenopausal == "Yes") & (dL.Caucasian == 1),
        "all women": dL.Age > 0}
for y in ["TOT_BMD", "HTOT_BMD", "NECK_BMD"]:
    for sub, mask in SUBS.items():
        fit_single(dL[mask], y, ["Age", "BMI", "Exercise_cat"], "LOS (BGI)", sub, "age + BMI + exercise")
        fit_single(dL[mask], y, ["Age", "Weight_kg"], "LOS (BGI)", sub, "age + weight")
        for thr in [5e-3, 2e-3]:
            fit_joint(ML, list(dL.loc[mask, "IID"]), dL, thr, y, ["Age", "BMI", "Exercise_cat"],
                      "LOS (BGI)", sub, "age + BMI + exercise")

os.makedirs(P("results"), exist_ok=True)
pd.DataFrame(res).to_csv(P("results/bvulgatus_lookup.csv"), index=False)
