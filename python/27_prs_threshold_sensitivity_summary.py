# ======================================================================
# 27_prs_threshold_sensitivity_summary.py — Summary of the alternative-threshold runs
# #
# Collects association and Delta R2 of the external PRS at the three thresholds.
# #
# Inputs : results/sensitivity_threshold/<run>/
# Outputs: results/sensitivity_threshold/summary_3_thresholds_holdout_lc.csv
# Paper  : Table S7 (panel H)
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
import csv
from pathlib import Path
import sys

run = Path(sys.argv[1] if len(sys.argv) > 1 else "results").resolve()
sensitivity = run / "sensitivity_threshold"
thresholds = (
    ("0.000001", "PT1e-6", "PT_1e-6"),
    ("0.0001", "PT1e-4", "PT_1e-4"),
    ("0.1", "ctrl_PT0.1", "PT_0.1"),
)
cohorts = (
    ("Internal holdout", "holdout_score_association.csv", "holdout_ci.csv", "dR2_vs_S2", 351),
    ("LC subset", "lc_score_association.csv", "lc_ci.csv", "dR2", 495),
)


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def summary(folder, assoc_name, ci_name, delta, expected_n):
    assoc, = (r for r in rows(folder / assoc_name)
              if r["score"] == "PRS_qian" and int(r["n"]) == expected_n)
    ci, = (r for r in rows(folder / ci_name)
           if r["algoritmo"] == "ElasticNet" and r["especificacion"] == "S3q_base_PRSqian")
    return [int(assoc["n"]), *(assoc[k] for k in ("beta", "IC95_inf", "IC95_sup", "P")),
            *(ci[k] for k in (delta, "IC95_inf", "IC95_sup"))]


counts = {r["PT"]: r["variantes"] for r in rows(run / "prs_external_threshold_selection.csv")}
output = []
for threshold, directory, pt in thresholds:
    folder = sensitivity / directory
    for label, assoc_name, ci_name, delta, n in cohorts:
        values = summary(folder, assoc_name, ci_name, delta, n)
        if directory == "ctrl_PT0.1":
            assert values == summary(run, assoc_name, ci_name, delta, n)
        output.append([threshold, counts[pt], label, *values])

with (sensitivity / "summary_3_thresholds_holdout_lc.csv").open("w") as f:
    writer = csv.writer(f)
    writer.writerow(["PT", "variants", "cohort", "n", "beta", "lo", "hi", "P", "dR2", "dlo", "dhi"])
    writer.writerows(output)
