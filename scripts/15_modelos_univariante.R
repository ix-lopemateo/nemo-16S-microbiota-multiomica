#!/usr/bin/env Rscript
# Cribado univariante de candidatas (clinicas y omicas)

library(readxl)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(corrplot)
library(knitr)
library(kableExtra)

set.seed(2026)

# Wide clínico AMPLIADO con las variables multiómicas (microbiota 16S, función
# PICRUSt2 y epigenética) generadas por 14_exportar_features_a_wide.R.
df <- read_excel("../data/nemo_wide_clean_con_leyenda_multiomics.xlsx",
                 sheet = "BASE DE DATOS")
cat("Dataset cargado")

# Excluir díadas marcadas
n_excl <- sum(df$mat_excluded == 1, na.rm = TRUE)
if (n_excl > 0) {
  excl_ids <- df$dyad_id[df$mat_excluded == 1]
  cat("Díadas excluidas:", paste(excl_ids, collapse = ", "), "\n")
  df <- df %>% filter(mat_excluded == 0 | is.na(mat_excluded))
}
cat("n analítica:", nrow(df), "\n")

# Codificación de factores
df <- df %>%
  mutate(
    inf_sex          = factor(inf_sex, levels = c("Female", "Male")),
    delivery_mode    = factor(delivery_mode, levels = c("Vaginal", "Cesarean")),
    feeding_6m = factor(feeding_6m, levels = c("Exclusive_BF",
                                            "BF_complementary",
                                            "Exclusive_Formula")),
    mat_bmi_cat      = factor(mat_bmi_cat, levels = c("Underweight", "Normal",
                                                       "Overweight", "Obese")),
    mat_gwg_adequacy = factor(mat_gwg_adequacy, levels = c("Insufficient",
                                                            "Adequate",
                                                            "Excessive"))
  )

# Transformaciones logarítmicas
df <- df %>%
  mutate(
    log_homa_ir = log(mat_homa_ir),
    log_pcr     = log(mat_pcr)
  )

n_total   <- nrow(df)
n_outcome <- sum(!is.na(df$dWLZ_0_6))
n_wlz6    <- sum(!is.na(df$WLZ_6))
max_pred  <- floor(n_outcome / 10)

continuous_candidates <- c(
  "mat_pregest_bmi",       # D1: Antropometría materna
  "mat_weight_gain",       # D1: Antropometría materna
  "log_homa_ir",           # D2: Metabolismo materno
  "mat_hba1c_pct",         # D2: Metabolismo materno
  "mat_tag",               # D2: Metabolismo materno (lipídico)
  "mat_apob_apoa1_ratio",  # D2: Metabolismo materno (ratio aterogénico)
  "log_pcr",               # D3: Inflamación
  "mat_age"                # D3: Clínica
)

categorical_candidates <- c(
  "delivery_mode",     # D4: Factores perinatales
  "inf_sex",           # D4: Factores perinatales (confusor DAG)
  "feeding_6m"         # D5: Alimentación infantil
)

multiomic_candidates <- grep("^(mb|div|met)_", names(df), value = TRUE)

# met_b1m_NRF1 se mide a 1 mes (b1m), un punto temporal ajeno al marco del estudio
# (madre basal, meconio y 6 meses); se excluye del cribado univariante.
multiomic_candidates <- setdiff(multiomic_candidates, "met_b1m_NRF1")

continuous_candidates <- c(continuous_candidates, multiomic_candidates)

all_candidates <- c(continuous_candidates, categorical_candidates)

cat("Candidatas continuas:", length(continuous_candidates),
    "(de ellas", length(multiomic_candidates), "multiómicas) | categóricas:",
    length(categorical_candidates), "\n")

summary_rows <- list()

for (var in continuous_candidates) {
  x <- df[[var]]
  summary_rows[[var]] <- tibble(
    Variable    = var,
    Tipo        = "Continua",
    n           = sum(!is.na(x)),
    `% ausente` = round(mean(is.na(x)) * 100, 1),
    Media       = round(mean(x, na.rm = TRUE), 2),
    DE          = round(sd(x, na.rm = TRUE), 2),
    Mediana     = round(median(x, na.rm = TRUE), 2),
    Mín         = round(min(x, na.rm = TRUE), 2),
    Máx         = round(max(x, na.rm = TRUE), 2)
  )
}

for (var in categorical_candidates) {
  x <- df[[var]]
  tbl <- table(x, useNA = "no")
  summary_rows[[var]] <- tibble(
    Variable    = var,
    Tipo        = "Categórica",
    n           = sum(!is.na(x)),
    `% ausente` = round(mean(is.na(x)) * 100, 1),
    Media       = NA_real_,
    DE          = NA_real_,
    Mediana     = NA_real_,
    Mín         = NA_real_,
    Máx         = NA_real_
  )
}

candidate_summary <- bind_rows(summary_rows)

kable(candidate_summary,
      caption = paste0("Resumen descriptivo de las ", nrow(candidate_summary),
                       " variables candidatas (clínicas + multiómicas).")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE)

for (var in categorical_candidates) {
  cat(var, ":\n")
  print(table(df[[var]], useNA = "ifany"))
  cat("\n")
}

cor_data <- df %>% select(all_of(continuous_candidates))

# Con el bloque multiómico añadido hay variables con cobertura disjunta (p. ej.
# madre vs bebé), por lo que `complete.cases` global vaciaría la matriz: se usa
# correlación por pares. Los pares sin solape suficiente quedan como 0 para poder
# dibujar la matriz.
cor_mat <- cor(cor_data, method = "spearman", use = "pairwise.complete.obs")
cor_mat[is.na(cor_mat)] <- 0

cat("Predictores continuos en la matriz:", ncol(cor_mat), "\n")

# Con muchas variables se omiten los coeficientes numéricos para legibilidad
mostrar_coef <- if (ncol(cor_mat) <= 12) "black" else NA
corrplot(cor_mat, method = "color", type = "lower",
         tl.cex = 0.5, tl.col = "black",
         addCoef.col = mostrar_coef, number.cex = 0.45,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         title = "Correlaciones de Spearman — predictores candidatos",
         mar = c(0, 0, 2, 0))

cor_pairs <- data.frame(Var1 = character(), Var2 = character(),
                        r = numeric(), stringsAsFactors = FALSE)
for (i in 1:(ncol(cor_mat) - 1)) {
  for (j in (i + 1):ncol(cor_mat)) {
    r <- cor_mat[i, j]
    if (!is.na(r) && abs(r) > 0.30) {
      cor_pairs <- rbind(cor_pairs,
        data.frame(Var1 = rownames(cor_mat)[i],
                   Var2 = colnames(cor_mat)[j],
                   r = round(r, 3)))
    }
  }
}
cor_pairs <- cor_pairs %>% arrange(desc(abs(r)))

kable(cor_pairs, caption = "Pares de predictores con |r| > 0,30.") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)

run_univariate <- function(outcome_var, predictor, data) {
  formula <- as.formula(paste(outcome_var, "~", predictor))
  mod <- lm(formula, data = data)
  tidy_out <- tidy(mod, conf.int = TRUE) %>%
    filter(term != "(Intercept)")
  glance_out <- glance(mod)
  tidy_out %>%
    mutate(
      variable = predictor,
      n        = nobs(mod),
      adj_r2   = round(glance_out$adj.r.squared, 4),
      .before  = 1
    )
}

run_conditional <- function(predictor, data) {
  formula <- as.formula(paste("WLZ_6 ~ WLZ_0 +", predictor))
  mod <- lm(formula, data = data)
  tidy_out <- tidy(mod, conf.int = TRUE) %>%
    filter(term != "(Intercept)" & term != "WLZ_0")
  glance_out <- glance(mod)
  tidy_out %>%
    mutate(
      variable = predictor,
      n        = nobs(mod),
      adj_r2   = round(glance_out$adj.r.squared, 4),
      .before  = 1
    )
}

# Estandarizar las variables continuas (media=0, DE=1)
df_std <- df
for (var in continuous_candidates) {
  s <- sd(df_std[[var]], na.rm = TRUE)
  m <- mean(df_std[[var]], na.rm = TRUE)
  if (!is.na(s) && s > 0) {
    df_std[[var]] <- (df_std[[var]] - m) / s
  }
}

univariate_cond <- lapply(all_candidates, function(var) {
  tryCatch(run_conditional(var, df_std), error = function(e) NULL)
}) %>%
  bind_rows() %>%
  mutate(sort_p = ifelse(variable == "feeding_6m", 
                         min(p.value[variable == "feeding_6m"]), 
                         p.value)) %>%
  arrange(sort_p, variable, p.value) %>%
  select(-sort_p) %>%
  mutate(
    selected = ifelse(p.value < 0.20, "Sí", "No"),
    dag      = ifelse(variable == "inf_sex", "Sí", "No")
  )

univariate_cond %>%
  select(variable, term, estimate, conf.low, conf.high, p.value,
         n, adj_r2, selected, dag) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~round(., 3)),
         p.value = round(p.value, 4)) %>%
  kable(caption = "Modelo 1 — WLZ_6 ~ WLZ_0 + predictor. Umbral: p < 0,20.",
        col.names = c("Variable", "Término", "β est", "IC inf", "IC sup",
                       "p-valor", "n", "R² aj.", "p<0.20", "DAG")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE) %>%
  row_spec(which(univariate_cond$selected == "Sí" |
                 univariate_cond$dag == "Sí"),
           bold = TRUE)

selected_cond <- univariate_cond %>%
  filter(selected == "Sí" | dag == "Sí") %>%
  pull(variable) %>%
  unique()

univariate_cont <- lapply(continuous_candidates, function(var) {
  tryCatch(run_univariate("dWLZ_0_6", var, df_std), error = function(e) NULL)
}) %>% bind_rows()

univariate_cat <- lapply(categorical_candidates, function(var) {
  tryCatch(run_univariate("dWLZ_0_6", var, df), error = function(e) NULL)
}) %>% bind_rows()

univariate_all <- bind_rows(univariate_cont, univariate_cat) %>%
  mutate(sort_p = ifelse(variable == "feeding_6m", 
                         min(p.value[variable == "feeding_6m"]), 
                         p.value)) %>%
  arrange(sort_p, variable, p.value) %>%
  select(-sort_p) %>%
  mutate(
    selected = ifelse(p.value < 0.20, "Sí", "No"),
    dag      = ifelse(variable == "inf_sex", "Sí", "No")
  )

univariate_all %>%
  select(variable, term, estimate, conf.low, conf.high, p.value,
         n, adj_r2, selected, dag) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~round(., 3)),
         p.value = round(p.value, 4)) %>%
  kable(caption = "Modelo 2 — dWLZ_0_6 ~ predictor. Umbral: p < 0,20.",
        col.names = c("Variable", "Término", "β est", "IC inf", "IC sup",
                       "p-valor", "n", "R² aj.", "p<0.20", "DAG")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE) %>%
  row_spec(which(univariate_all$selected == "Sí" |
                 univariate_all$dag == "Sí"),
           bold = TRUE)

selected_vars <- univariate_all %>%
  filter(selected == "Sí" | dag == "Sí") %>%
  pull(variable) %>%
  unique()

# Estandarizar continuas
df_std2 <- df
for (var in continuous_candidates) {
  s <- sd(df_std2[[var]], na.rm = TRUE)
  m <- mean(df_std2[[var]], na.rm = TRUE)
  if (!is.na(s) && s > 0) {
    df_std2[[paste0(var, "_z")]] <- (df_std2[[var]] - m) / s
  }
}

std_vars2 <- paste0(continuous_candidates, "_z")
std_vars2 <- std_vars2[std_vars2 %in% names(df_std2)]

forest_cond_cont <- lapply(std_vars2, function(var) {
  tryCatch({
    mod <- lm(as.formula(paste("WLZ_6 ~ WLZ_0 +", var)), data = df_std2)
    tidy(mod, conf.int = TRUE) %>%
      filter(term == var) %>%
      mutate(variable = gsub("_z$", "", var),
             type = "Continua (por 1 DE)")
  }, error = function(e) NULL)
}) %>% bind_rows()

forest_cond_cat <- lapply(categorical_candidates, function(var) {
  tryCatch({
    mod <- lm(as.formula(paste("WLZ_6 ~ WLZ_0 +", var)), data = df)
    tidy(mod, conf.int = TRUE) %>%
      filter(term != "(Intercept)" & term != "WLZ_0") %>%
      mutate(variable = var, type = "Categórica")
  }, error = function(e) NULL)
}) %>% bind_rows()

forest_cond_all <- bind_rows(forest_cond_cont, forest_cond_cat) %>%
  mutate(
    sig   = case_when(
      p.value < 0.05 ~ "p < 0.05",
      p.value < 0.20 ~ "p < 0.20",
      TRUE           ~ "p >= 0.20"
    ),
    label = ifelse(type == "Categórica", term, variable)
  )

ggplot(forest_cond_all, aes(x = estimate, y = reorder(label, estimate),
                             color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  scale_color_manual(
    values = c("p < 0.05" = "#B2182B",
               "p < 0.20" = "#F4A582",
               "p >= 0.20" = "#999999"),
    name = "Significación"
  ) +
  labs(
    x = expression(beta ~ "(WLZ 6 meses | WLZ 0)"),
    y = NULL,
    title = "Asociaciones con WLZ_6 ajustado por WLZ_0 (Modelo 1)",
    subtitle = paste0("n = ", n_outcome,
                      " · Predictores continuos estandarizados (β por 1 DE)")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom"
  )

df_std <- df
for (var in continuous_candidates) {
  s <- sd(df_std[[var]], na.rm = TRUE)
  m <- mean(df_std[[var]], na.rm = TRUE)
  if (!is.na(s) && s > 0) {
    df_std[[paste0(var, "_z")]] <- (df_std[[var]] - m) / s
  }
}

std_cont_vars <- paste0(continuous_candidates, "_z")
std_cont_vars <- std_cont_vars[std_cont_vars %in% names(df_std)]

forest_data <- lapply(std_cont_vars, function(var) {
  tryCatch({
    mod <- lm(as.formula(paste("dWLZ_0_6 ~", var)), data = df_std)
    tidy(mod, conf.int = TRUE) %>%
      filter(term != "(Intercept)") %>%
      mutate(variable = gsub("_z$", "", var),
             type = "Continua (por 1 DE)")
  }, error = function(e) NULL)
}) %>% bind_rows()

forest_cat <- lapply(categorical_candidates, function(var) {
  tryCatch({
    mod <- lm(as.formula(paste("dWLZ_0_6 ~", var)), data = df)
    tidy(mod, conf.int = TRUE) %>%
      filter(term != "(Intercept)") %>%
      mutate(variable = var, type = "Categórica")
  }, error = function(e) NULL)
}) %>% bind_rows()

forest_all <- bind_rows(forest_data, forest_cat) %>%
  mutate(
    sig   = case_when(
      p.value < 0.05 ~ "p < 0.05",
      p.value < 0.20 ~ "p < 0.20",
      TRUE           ~ "p >= 0.20"
    ),
    label = ifelse(type == "Categórica", term, variable)
  )

ggplot(forest_all, aes(x = estimate, y = reorder(label, estimate),
                        color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  scale_color_manual(
    values = c("p < 0.05" = "#B2182B",
               "p < 0.20" = "#F4A582",
               "p >= 0.20" = "#999999"),
    name = "Significación"
  ) +
  labs(
    x = expression(beta ~ "(ΔZ WLZ 0-6 meses)"),
    y = NULL,
    title = "Asociaciones con ΔZ(WLZ) nacimiento → 6 meses (Modelo 2)",
    subtitle = paste0("n = ", n_outcome,
                      " · Predictores continuos estandarizados (β por 1 DE)")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom"
  )

max_pred <- floor(sum(!is.na(df$dWLZ_0_6)) / 10)
n_events <- sum(df$overweight_risk_6m == 1, na.rm = TRUE)
n_ow     <- sum(!is.na(df$overweight_risk_6m))

from_direct <- univariate_all %>%
  filter(selected == "Sí") %>%
  pull(variable) %>%
  unique()

from_cond <- univariate_cond %>%
  filter(selected == "Sí") %>%
  pull(variable) %>%
  unique()

resumen_rows <- list()

for (v in union(union(from_direct, from_cond), "inf_sex")) {
  origenes <- c()
  if (v %in% from_cond)   origenes <- c(origenes, "Modelo 1: WLZ_6 ~ WLZ_0 + predictor (p < 0.20)")
  if (v %in% from_direct) origenes <- c(origenes, "Modelo 2: dWLZ_0_6 ~ predictor (p < 0.20)")
  if (v == "inf_sex")     origenes <- c(origenes, "Confusor DAG (forzado)")
  resumen_rows[[v]] <- tibble(Variable = v, Origen = paste(origenes, collapse = " + "))
}

resumen <- bind_rows(resumen_rows)

kable(resumen,
      caption = "Variables que entran en la fase multivariante.") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)

dir.create("../R_output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("../R_output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("../R_output/models",  recursive = TRUE, showWarnings = FALSE)

write.csv(candidate_summary,
          "../R_output/tables/candidate_variable_summary_omico.csv", row.names = FALSE)
write.csv(round(cor_mat, 3),
          "../R_output/tables/correlation_matrix_candidates_omico.csv")
write.csv(univariate_cond,
          "../R_output/tables/univariate_screening_WLZ6_omico.csv", row.names = FALSE)
write.csv(univariate_all,
          "../R_output/tables/univariate_screening_dWLZ06_omico.csv", row.names = FALSE)

pdf("../R_output/figures/fig_corr_candidates_omico.pdf", width = 10, height = 9)
corrplot(cor_mat, method = "color", type = "lower",
         tl.cex = 0.65, tl.col = "black",
         addCoef.col = "black", number.cex = 0.55,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         title = "Correlaciones de Spearman — predictores candidatos",
         mar = c(0, 0, 2, 0))
dev.off()

# --- Objetos para el script 16_modelos_multivariante.Rmd ---
# Dataset analítico limpio
saveRDS(df, "../R_output/models/df_analytic_omico.rds")

# Variables seleccionadas por cada estrategia
selected_cond   <- univariate_cond %>%
  filter(selected == "Sí" | dag == "Sí") %>% pull(variable) %>% unique()
selected_direct <- univariate_all %>%
  filter(selected == "Sí" | dag == "Sí") %>% pull(variable) %>% unique()

saveRDS(list(
  selected_cond   = selected_cond,
  selected_direct = selected_direct,
  forced          = c("WLZ_0", "inf_sex"),
  univariate_cond = univariate_cond,
  univariate_all  = univariate_all,
  resumen         = resumen,
  multiomic_candidates = multiomic_candidates
), "../R_output/models/variable_selection_omico.rds")

cat("Objetos guardados en output/models/ para el script 16.\n")
