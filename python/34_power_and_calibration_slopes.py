# ======================================================================
# 34_power_and_calibration_slopes.py — Minimum detectable increments and calibration-slope stability
#
# Panel A: minimum detectable increment in R2 for one added predictor (non-central F, 80% power,
#          alpha = 0.05) and empirically as 2.8 x the bootstrap SE of Delta R2 (Elastic Net).
# Panel G: calibration slope (observed on predicted BMD) across 2,000 bootstrap resamples, Elastic Net.
#
# Inputs : results/holdout_ci.csv, lc_ci.csv, guangzhou_ci.csv, calibration_holdout_lc.csv,
#          guangzhou_calibration.csv, holdout_predictions.csv, lc_predictions.csv, guangzhou_predictions.csv
# Outputs: results/power_mdes.csv, results/calibration_slope_bootstrap.csv
# Paper  : Table S7 (panels A and G)
# Run from the repository root.
# ======================================================================
import os
import numpy as np, pandas as pd
from scipy import stats
from scipy.optimize import brentq

RES = os.path.join(os.environ.get("HIPBMD_ROOT", "."), "results")

def mdes(n, k, r2_base, power=0.8, alpha=0.05):
    df2 = n - k - 2
    fcrit = stats.f.ppf(1 - alpha, 1, df2)
    lam = brentq(lambda L: 1 - stats.ncf.cdf(fcrit, 1, df2, L) - power, 0.01, 100)
    f2 = lam / n
    return f2 * (1 - r2_base) / (1 + f2)

cal = pd.read_csv(os.path.join(RES, "calibration_holdout_lc.csv"))
gz = pd.read_csv(os.path.join(RES, "guangzhou_calibration.csv"))
r2_base = {"holdout": cal.query("cohorte == 'holdout' and modelo == 'S2_base__ElasticNet'").R2.item(),
           "LC": cal.query("cohorte == 'LC' and modelo == 'S2_base__ElasticNet'").R2.item(),
           "Guangzhou": gz.query("modelo == 'S2_base__ElasticNet'").R2_recal.item()}
n = {"holdout": 351, "LC": 495, "Guangzhou": 492}
k = {"holdout": 16, "LC": 13, "Guangzhou": 5}     # estimable base predictors in each cohort
ci = {"holdout": "holdout_ci.csv", "LC": "lc_ci.csv", "Guangzhou": "guangzhou_ci.csv"}
rows = []
for coh in ("holdout", "LC", "Guangzhou"):
    e = pd.read_csv(os.path.join(RES, ci[coh])).query("algoritmo == 'ElasticNet'").set_index("especificacion")
    se = {s: (e.loc[s, "IC95_sup"] - e.loc[s, "IC95_inf"]) / 3.92 for s in ("S3_base_PRS", "S3q_base_PRSqian", "S4_base_MRS")}
    rows.append(dict(cohort=coh, n=n[coh], base_predictors=k[coh], r2_base=round(r2_base[coh], 2),
                     mdes_f_test=round(mdes(n[coh], k[coh], r2_base[coh]), 4),
                     mdes_boot_los_prs=round(2.8 * se["S3_base_PRS"], 4),
                     mdes_boot_external_prs=round(2.8 * se["S3q_base_PRSqian"], 4),
                     mdes_boot_mrs=round(2.8 * se["S4_base_MRS"], 4)))
pd.DataFrame(rows).to_csv(os.path.join(RES, "power_mdes.csv"), index=False)

out = []
for coh, f in (("Guangzhou", "guangzhou_predictions.csv"), ("LC", "lc_predictions.csv"), ("holdout", "holdout_predictions.csv")):
    p = pd.read_csv(os.path.join(RES, f)); y = p["y"].values
    rng = np.random.default_rng(2026); idx = [rng.integers(0, len(y), len(y)) for _ in range(2000)]
    for c in [c for c in p.columns if c.endswith("ElasticNet")]:
        sl = np.array([np.polyfit(p[c].values[i], y[i], 1)[0] for i in idx])
        out.append(dict(cohort=coh, specification=c, min_slope=sl.min(), resamples_slope_le_0=int((sl <= 0).sum())))
s = pd.DataFrame(out)
s.groupby("cohort").agg(specifications=("specification", "count"), min_slope=("min_slope", "min"),
                        resamples_slope_le_0=("resamples_slope_le_0", "sum")).round(2).to_csv(
    os.path.join(RES, "calibration_slope_bootstrap.csv"))
