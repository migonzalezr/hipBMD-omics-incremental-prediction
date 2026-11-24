# ======================================================================
# figS4_prs_deciles.R — Figure S4
# #
# Total hip BMD across deciles of each polygenic score in modeling discovery and the internal holdout.
# #
# Inputs : results/*_predictions.csv, PRS score files
# Outputs: figures/png, figures/pdf
# Paper  : Figure S4
# Run from the repository root; paths are defined in config/paths.R.
# ======================================================================
library(data.table)
library(ggplot2)

source("config/paths.R")
root <- BASE
analysis <- file.path(BASE, "derived")
results <- file.path(BASE, "results")
cohorts <- c("Modeling discovery", "Internal holdout")
scores <- c("UK Biobank PRS", "LOS PRS")
file_labels <- c("discovery_predictions_oof.csv", "holdout_predictions.csv")
expected_n <- c(1430L, 351L)

samples <- rbindlist(lapply(seq_along(cohorts), function(i) {
  p <- fread(file.path(results, file_labels[i]), select = c("IID", "y"),
             colClasses = c(IID = "character"))
  p[, cohort := cohorts[i]]
  stopifnot(nrow(p) == expected_n[i], all(is.finite(p$y)))
  p
}))
stopifnot(!anyDuplicated(samples$IID))
samples[, cohort := factor(cohort, levels = cohorts)]
read_score <- function(path, name) {
  score <- fread(file.path(analysis, path), select = c("IID", "SCORE1_AVG"),
                 colClasses = c(IID = "character"))[IID %in% samples$IID]
  stopifnot(nrow(score) == nrow(samples), !anyDuplicated(score$IID),
            all(is.finite(score$SCORE1_AVG)))
  setnames(score, "SCORE1_AVG", name)
}
samples <- merge(samples, read_score("prs_external/external_prs.PT_0.1.sscore", "PRS_qian"), by = "IID")
samples <- merge(samples, read_score("prs_los/los_prs.PT_0.01.sscore", "PRS"), by = "IID")
d <- melt(samples, id.vars = c("IID", "cohort", "y"), measure.vars = c("PRS_qian", "PRS"),
          variable.name = "score_key", value.name = "score_value")
d[, score := factor(score_key, labels = scores)]
setorder(d, cohort, score, score_value, IID)
d[, decile := 1L + ((seq_len(.N) - 1L) * 10L) %/% .N, by = .(cohort, score)]
print(d[, .(tied_cutpoints = sum(score_value[-1L] == score_value[-.N] &
       decile[-1L] != decile[-.N])), by = .(cohort, score)])

means <- d[, .(n = .N, mean_BMD = mean(y), sd_BMD = sd(y)), by = .(score, cohort, decile)]
means[, margin := qt(0.975, n - 1L) * sd_BMD / sqrt(n)]
means[, `:=`(CI95_low = mean_BMD - margin, CI95_high = mean_BMD + margin, margin = NULL)]
stopifnot(nrow(means) == 40L, all(means$n %in% c(35L, 36L, 143L)))

clinical <- c("Age", "gender_dummy", "Weight_kg", "Height_cm", "Race_African",
              "Race_Caucasian", "Race_Asian", "Exercise_Regular_cat", "Smoke_cat",
              "Sun_Exposure_cat", "Alcohol_Drinking_cat", "Cheese_Eat_cat",
              "Milk_Drink_cat", "Yogurt_Eat_cat")
cl <- fread(file.path(root, "data/phenotypes/los_clinical_imputed.csv"),
            select = c("Subject_ID", clinical), colClasses = c(Subject_ID = "character"))
setnames(cl, "Subject_ID", "IID")
cl[, IID := trimws(IID)]
cl <- unique(cl[IID %in% samples$IID], by = "IID")
pcs <- fread(file.path(root, "data/genotypes/los_qc_pca.eigenvec"),
             header = FALSE, select = 2:4, colClasses = list(character = 2L))
setnames(pcs, c("IID", "gPC1", "gPC2"))
pcs[, IID := trimws(IID)]
pcs <- pcs[IID %in% samples$IID]
covariates <- c(clinical, "gPC1", "gPC2")
x <- merge(cl, pcs, by = "IID")
stopifnot(nrow(x) == nrow(samples), !anyDuplicated(x$IID), all(complete.cases(x)))
analysis_data <- merge(d, x, by = "IID")
discovery_sd <- analysis_data[cohort == cohorts[1], .(score_sd = sd(score_value)), by = score_key]
analysis_data <- merge(analysis_data, discovery_sd, by = "score_key")
analysis_data[, score_per_discovery_sd := score_value / score_sd]

fit_term <- function(data, term, adjust) {
  fit <- lm(reformulate(c(term, adjust), response = "y"), data = data)
  coef <- summary(fit)$coefficients[term, ]
  ci <- confint(fit, term)
  c(beta = unname(coef["Estimate"]), lower = ci[1], upper = ci[2],
    p = unname(coef["Pr(>|t|)"]))
}
continuous <- analysis_data[, as.list(fit_term(.SD, "score_per_discovery_sd", covariates)),
                            by = .(cohort, score_key)]
ref <- fread(file.path(results, "holdout_score_association.csv"))[score %chin% c("PRS", "PRS_qian")]
ref[, cohort := fifelse(conjunto == "discovery", cohorts[1], cohorts[2])]
check <- merge(continuous, ref, by.x = c("cohort", "score_key"), by.y = c("cohort", "score"))
stopifnot(nrow(check) == 4L, all(round(check$beta.x, 5) == check$beta.y))

trends <- analysis_data[, {
  a <- fit_term(.SD, "decile", character())
  b <- fit_term(.SD, "decile", covariates)
  c(list(n = .N), setNames(as.list(a), paste0("unadjusted_", names(a))),
    setNames(as.list(b), paste0("adjusted_", names(b))))
}, by = .(score, cohort)]
print(trends)
setorder(means, score, cohort, decile)
fwrite(means, file.path(results, "PRS_deciles_summary.csv"))
fwrite(trends, file.path(results, "PRS_deciles_trends.csv"))

k <- position_dodge(width = 0.26)
cohort_labels <- c("Modeling discovery (n = 1,430)", "Internal holdout (n = 351)")
plot <- ggplot(means, aes(decile, mean_BMD, colour = cohort, shape = cohort, group = cohort)) +
  geom_errorbar(aes(ymin = CI95_low, ymax = CI95_high), width = 0.14,
                linewidth = 0.62, position = k, alpha = 0.70) +
  geom_line(linewidth = 0.7, position = k) +
  geom_point(size = 2.8, position = k) +
  facet_wrap(~score, nrow = 1) +
  scale_colour_manual(values = c("#202020", "#2783B5"), labels = cohort_labels) +
  scale_shape_manual(values = c(16, 17), labels = cohort_labels) +
  scale_x_continuous(breaks = 1:10, minor_breaks = NULL) +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
  labs(x = "PRS decile (lowest to highest)", y = "Mean total hip BMD (g/cm²)",
       colour = NULL, shape = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        panel.grid.major.y = element_line(colour = "#e6ebef", linewidth = 0.4),
        panel.spacing = grid::unit(2.4, "lines"))

name <- "figS4_prs_deciles"
ggsave(file.path(BASE, "figures", "png", paste0(name, ".png")), plot,
       width = 10.5, height = 4.6, dpi = 320, bg = "white")
pdf_device <- if (Sys.info()[["sysname"]] == "Darwin") {
  function(filename, ...) quartz(type = "pdf", file = filename, ...)
} else {
  cairo_pdf
}
ggsave(file.path(BASE, "figures", "pdf", paste0(name, ".pdf")), plot,
       width = 10.5, height = 4.6, device = pdf_device, bg = "white")
