
install.packages("tidyverse")
install.packages("caret")
suppressPackageStartupMessages({
  library(tidyverse)   # data handling, ggplot2
  library(pheatmap)    # heatmap
  library(caret)       # stratified split, confusion matrix
  library(e1071)       # SVM with tuning
})

set.seed(2026)
FIG <- "C:/Users/User/Desktop/Gawesha/amr-resistance-r/amr-resistance-r/results/figures"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 11))

# ---- 1. Import ------------------------------------------------------------------
load("C:/Users/User/Desktop/Gawesha/amr-resistance-r/amr-resistance-r/data/example_isolates.rda")     
load("C:/Users/User/Desktop/Gawesha/amr-resistance-r/amr-resistance-r/data/antimicrobials.rda")       
raw <- example_isolates
drug_cols <- names(raw)[7:46]         

# Drug class: first listed class of each drug code (some drugs have several)
am <- antimicrobials[!duplicated(antimicrobials$ab), ]
drug_class <- setNames(
  vapply(am$group[match(drug_cols, am$ab)], function(x) x[1], character(1)),
  drug_cols)
drug_name <- setNames(am$name[match(drug_cols, am$ab)], drug_cols)

# Species codes (AMR package microorganism IDs) for the six most common species
species_lab <- c(B_ESCHR_COLI = "Escherichia coli",
                 B_STPHY_CONS = "Coag.-negative staphylococci",
                 B_STPHY_AURS = "Staphylococcus aureus",
                 B_STPHY_EPDR = "Staphylococcus epidermidis",
                 B_STRPT_PNMN = "Streptococcus pneumoniae",
                 B_STPHY_HMNS = "Staphylococcus haemolyticus")

# ---- 2. Quality control ----------------------------------------------------------
stopifnot(all(drug_cols %in% am$ab))
stopifnot(all(unlist(raw[drug_cols]) %in% c("S", "I", "R", NA)))

qc <- tibble(
  n_isolates          = nrow(raw),
  n_patients          = n_distinct(raw$patient),
  n_drug_columns      = length(drug_cols),
  date_min            = as.character(min(as.Date(raw$date))),
  date_max            = as.character(max(as.Date(raw$date))),
  n_missing_age       = sum(is.na(raw$age)),
  n_missing_gender    = sum(is.na(raw$gender)),
  n_dup_patient_species_date = sum(duplicated(raw[, c("patient", "mo", "date")])),
  n_species           = n_distinct(raw$mo)
)
write_csv(qc, "results/qc_summary.csv")
print(qc)

# Missingness per drug (not tested)
missing_drug <- tibble(drug = drug_cols,
                       n_tested = colSums(!is.na(raw[drug_cols])),
                       pct_tested = round(100 * n_tested / nrow(raw), 1))
write_csv(missing_drug, "results/drug_test_coverage.csv")

# ---- 3. Cleaning: first isolate per patient per species ----------------------------
iso <- raw |>
  mutate(date = as.Date(date),
         species = recode(mo, !!!species_lab, .default = "Other species")) |>
  arrange(date) |>
  group_by(patient, mo) |>
  slice(1) |>
  ungroup()
cat(sprintf("Isolates before first-isolate filter: %d; after: %d\n", nrow(raw), nrow(iso)))

# ---- 4. Prevalence of resistance ------------------------------------------------------
prev <- map_dfr(drug_cols, function(d) {
  x <- iso[[d]]; n <- sum(!is.na(x)); r <- sum(x == "R", na.rm = TRUE)
  ci <- if (n > 0) binom.test(r, n)$conf.int else c(NA, NA)
  tibble(drug = d, name = drug_name[[d]], class = drug_class[[d]],
         n_tested = n, n_resistant = r,
         pct_resistant = round(100 * r / n, 1),
         ci_low = round(100 * ci[1], 1), ci_high = round(100 * ci[2], 1))
}) |> arrange(desc(pct_resistant))
write_csv(prev, "results/drug_resistance_prevalence.csv")

# Fig 1: resistance prevalence by drug (drugs with >= 100 tests)
p1 <- prev |>
  filter(n_tested >= 100) |>
  mutate(drug = fct_reorder(drug, pct_resistant)) |>
  ggplot(aes(pct_resistant, drug)) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, colour = "grey50") +
  geom_point(colour = "#0072B2") +
  labs(title = "Fig 1. Resistance prevalence by antimicrobial",
       subtitle = "Among first isolates tested; bars are 95% exact confidence intervals",
       x = "% resistant (R)", y = NULL)
ggsave(file.path(FIG, "fig1_prevalence.png"), p1, width = 7, height = 8, dpi = 150)

# Fig 2: heatmap, % resistant by species (rows) and drug (columns)
top_drugs <- prev |> filter(n_tested >= 800) |> pull(drug)
top_species <- iso |> count(species, sort = TRUE) |>
  filter(species != "Other species") |> slice_head(n = 4) |> pull(species)
hm <- map_dfr(top_species, function(sp) {
  sub <- iso |> filter(species == sp)
  map_dfr(top_drugs, function(d) {
    x <- sub[[d]]; n <- sum(!is.na(x))
    tibble(species = sp, drug = d,
           pct = if (n >= 20) 100 * mean(x == "R", na.rm = TRUE) else NA_real_)
  })
})
hm_mat <- hm |> pivot_wider(names_from = drug, values_from = pct) |>
  column_to_rownames("species") |> as.matrix()
pheatmap(hm_mat, cluster_rows = FALSE, cluster_cols = TRUE,
         display_numbers = TRUE, number_format = "%.0f",
         color = colorRampPalette(c("white", "#D55E00"))(100),
         main = "Fig 2. % resistant by species and antimicrobial (species with n >= 20 per drug)",
         fontsize = 8, na_col = "grey90",
         filename = file.path(FIG, "fig2_species_heatmap.png"),
         width = 10, height = 4.5)
write_csv(as_tibble(hm_mat, rownames = "species"), "results/species_resistance_table.csv")

# ---- 5. Multidrug resistance (MDR) -------------------------------------------------------
res_class <- function(row) {
  classes <- unique(drug_class[drug_cols])
  tested <- sapply(classes, function(cl) any(!is.na(row[names(drug_class)[drug_class == cl]])))
  resist <- sapply(classes, function(cl) any(row[names(drug_class)[drug_class == cl]] == "R", na.rm = TRUE))
  c(n_tested = sum(tested), n_resist = sum(resist & tested))
}
cls_counts <- t(apply(iso[, drug_cols], 1, res_class))
iso <- iso |>
  mutate(n_classes_tested = cls_counts[, "n_tested"],
         n_classes_resist = cls_counts[, "n_resist"],
         MDR = case_when(n_classes_tested >= 3 & n_classes_resist >= 3 ~ 1L,
                         n_classes_tested >= 3 ~ 0L,
                         TRUE ~ NA_integer_))
cat(sprintf("MDR defined for %d isolates; MDR prevalence %.1f%%\n",
            sum(!is.na(iso$MDR)), 100 * mean(iso$MDR, na.rm = TRUE)))

mdr_ward <- iso |> filter(!is.na(MDR)) |>
  group_by(ward) |>
  summarise(n = n(), n_MDR = sum(MDR), pct_MDR = round(100 * mean(MDR), 1), .groups = "drop")
write_csv(mdr_ward, "results/mdr_by_ward.csv")

# Fig 3: MDR by ward
p3 <- iso |> filter(!is.na(MDR)) |>
  group_by(ward) |>
  summarise(pct = 100 * mean(MDR), n = n(), se = 100 * sqrt(mean(MDR) * (1 - mean(MDR)) / n()),
            .groups = "drop") |>
  ggplot(aes(ward, pct)) +
  geom_col(fill = "#0072B2", width = 0.6) +
  geom_errorbar(aes(ymin = pct - 1.96 * se, ymax = pct + 1.96 * se), width = 0.2) +
  geom_text(aes(label = paste0("n=", n)), vjust = -2.5, size = 3) +
  labs(title = "Fig 3. MDR prevalence by ward type",
       subtitle = "Error bars: approximate 95% CI; MDR = resistant to >= 3 classes",
       x = "Ward", y = "% MDR")
ggsave(file.path(FIG, "fig3_mdr_by_ward.png"), p3, width = 6, height = 5, dpi = 150)

# Association: logistic regression for MDR, adjusted for ward, age, gender, species
mdr_dat <- iso |>
  filter(!is.na(MDR), !is.na(age), !is.na(gender)) |>
  mutate(ward = factor(ward, levels = c("Outpatient", "Clinical", "ICU")),
         gender = factor(gender),
         # Species with almost no non-MDR isolates (E. coli, S. pneumoniae) cause
         # complete separation, so species enters as S. aureus versus all others.
         s_aureus = as.integer(species == "Staphylococcus aureus"))
fit_mdr <- glm(MDR ~ ward + age + gender + s_aureus, data = mdr_dat, family = binomial)
or_tbl <- tibble(term = names(coef(fit_mdr)),
                 odds_ratio = exp(coef(fit_mdr)),
                 ci_low = exp(confint.default(fit_mdr))[, 1],
                 ci_high = exp(confint.default(fit_mdr))[, 2],
                 p_value = summary(fit_mdr)$coefficients[, 4]) |>
  mutate(across(where(is.numeric), ~ signif(.x, 3)))
write_csv(or_tbl, "results/mdr_logistic_odds_ratios.csv")
print(or_tbl)

# ---- 6. Prediction of ciprofloxacin (CIP) resistance ---------------------------------
feat_drugs <- prev |>
  filter(drug != "CIP", class != "Fluoroquinolones", n_tested >= 1000) |>
  pull(drug)
pred_dat <- iso |> filter(!is.na(CIP))
pred_dat <- pred_dat[complete.cases(pred_dat[, feat_drugs]), ]
X <- as.data.frame(lapply(pred_dat[, feat_drugs], function(x) as.integer(x == "R")))
y <- factor(ifelse(pred_dat$CIP == "R", "R", "nonR"), levels = c("nonR", "R"))
cat(sprintf("CIP model: %d isolates, %d features, %d CIP-R (%.1f%%)\n",
            nrow(X), ncol(X), sum(y == "R"), 100 * mean(y == "R")))

idx <- createDataPartition(y, p = 0.7, list = FALSE)
Xtr <- X[idx, ]; Xte <- X[-idx, ]; ytr <- y[idx]; yte <- y[-idx]

# Baseline: always predict the majority class
baseline_acc <- max(prop.table(table(yte)))

fit_log <- glm(ytr ~ ., data = cbind(Xtr, ytr = as.integer(ytr == "R")), family = binomial)
p_log <- predict(fit_log, newdata = Xte, type = "response")

tuned <- tune.svm(x = Xtr, y = ytr, kernel = "radial",
                  cost = 2^(-2:4), gamma = 2^(-6:0) / ncol(Xtr))
svm_fit <- svm(x = Xtr, y = ytr, kernel = "radial",
               cost = tuned$best.parameters$cost,
               gamma = tuned$best.parameters$gamma, probability = TRUE)
svm_pred <- predict(svm_fit, Xte, probability = TRUE)
p_svm <- attr(svm_pred, "probabilities")[, "R"]

auc_fn <- function(score, y01) {
  r <- rank(score); n1 <- sum(y01 == 1); n0 <- sum(y01 == 0)
  (sum(r[y01 == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
eval_model <- function(name, p) {
  cls <- factor(ifelse(p >= 0.5, "R", "nonR"), levels = c("nonR", "R"))
  cm <- confusionMatrix(cls, yte, positive = "R")
  tibble(model = name,
         accuracy = unname(cm$overall["Accuracy"]),
         sensitivity = unname(cm$byClass["Sensitivity"]),
         specificity = unname(cm$byClass["Specificity"]),
         auc = auc_fn(p, as.integer(yte == "R")),
         n_test = length(yte))
}
perf <- bind_rows(eval_model("Logistic regression", p_log),
                  eval_model("RBF SVM (tuned)", p_svm)) |>
  mutate(baseline_accuracy = baseline_acc,
         across(where(is.numeric), ~ round(.x, 3)))
write_csv(perf, "results/cip_model_performance.csv")
print(perf)

# ---- 7. Session information --------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
cat("Analysis complete. Outputs written to results/.\n")
