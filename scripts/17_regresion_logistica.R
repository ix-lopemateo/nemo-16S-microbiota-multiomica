#!/usr/bin/env Rscript
# Regresion logistica penalizada de Firth (riesgo de sobrepeso)

library(tidyverse)
library(readxl)
library(logistf)
library(broom)
library(pROC)
library(boot)
library(knitr)
library(kableExtra)
library(patchwork)
library(ggplot2)

df  <- readRDS("../R_output/models/df_analytic_clinico.rds")
sel <- readRDS("../R_output/models/variable_selection_clinico.rds")
cat("Dataset cargado: n =", nrow(df), "\n")

df <- df %>%
  mutate(
    overweight_risk_6m = ifelse(WLZ_6 > 1, 1L, 0L)
  )

# Subconjuntos de trabajo
df_lr <- df %>% filter(!is.na(overweight_risk_6m))
df_cc <- df_lr %>%
  filter(complete.cases(overweight_risk_6m, mat_weight_gain, inf_sex))

dir.create("../R_output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("../R_output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("../R_output/models",  recursive = TRUE, showWarnings = FALSE)

n_events    <- sum(df$overweight_risk_6m == 1, na.rm = TRUE)
n_nonevents <- sum(df$overweight_risk_6m == 0, na.rm = TRUE)
n_total_lr  <- n_events + n_nonevents
prevalence  <- round(n_events / n_total_lr * 100, 1)
max_pred    <- floor(n_events / 10)

freq_table <- tibble(
  Categoría = c("Sin riesgo (WLZ ≤ +1)",
                "Riesgo de sobrepeso (WLZ > +1)",
                "Total"),
  n         = c(n_nonevents, n_events, n_total_lr),
  `%`       = c(round(n_nonevents / n_total_lr * 100, 1), prevalence, 100)
)

kable(freq_table,
      caption = "Distribución del desenlace binario.") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = FALSE)

continuous_vars <- c("mat_weight_gain", "mat_pregest_bmi", "log_homa_ir",
                     "mat_hba1c_pct", "mat_tag", "mat_apob_apoa1_ratio",
                     "log_pcr", "mat_age")

desc_cont <- df_lr %>%
  select(overweight_risk_6m, all_of(continuous_vars)) %>%
  pivot_longer(-overweight_risk_6m, names_to = "variable",
               values_to = "valor") %>%
  filter(!is.na(valor)) %>%
  group_by(variable, overweight_risk_6m) %>%
  summarise(
    n       = n(),
    media   = round(mean(valor), 2),
    DE      = round(sd(valor), 2),
    mediana = round(median(valor), 2),
    IQR     = round(IQR(valor), 2),
    .groups = "drop"
  ) %>%
  mutate(Grupo = ifelse(overweight_risk_6m == 1, "Riesgo", "Sin riesgo"))

kable(desc_cont %>% select(variable, Grupo, n, media, DE, mediana, IQR),
      caption = "Descriptiva de predictores continuos según desenlace binario.") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = FALSE)

p_box <- df_lr %>%
  select(overweight_risk_6m, all_of(continuous_vars)) %>%
  mutate(overweight_risk_6m = factor(overweight_risk_6m,
                                     labels = c("Sin riesgo", "Riesgo"))) %>%
  pivot_longer(-overweight_risk_6m, names_to = "variable",
               values_to = "valor") %>%
  filter(!is.na(valor)) %>%
  ggplot(aes(x = overweight_risk_6m, y = valor,
             fill = overweight_risk_6m)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Sin riesgo" = "#2166AC",
                               "Riesgo"     = "#B2182B")) +
  labs(x = NULL, y = "Valor", fill = "Desenlace",
       title = "Distribución de predictores continuos según riesgo de sobrepeso") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p_box

cat_vars <- c("inf_sex", "delivery_mode", "feeding_6m")

for (v in cat_vars) {
  tab <- table(df_lr[[v]], df_lr$overweight_risk_6m, useNA = "no")
  cat("\n---", v, "vs overweight_risk_6m ---\n")
  print(tab)
  ft <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
  cat("Fisher exact test p =", round(ft$p.value, 4), "\n")
}

all_candidates <- c(continuous_vars, cat_vars)

univariate_logit <- lapply(all_candidates, function(var) {
  
  df_temp <- df_lr %>%
    select(overweight_risk_6m, all_of(var)) %>%
    filter(complete.cases(.))
  
  if (length(unique(df_temp[[var]])) < 2) return(NULL)
  
  fml <- as.formula(paste("overweight_risk_6m ~", var))
  
  tryCatch({
    m <- logistf(fml, data = df_temp)
    
    terms <- names(coef(m))
    idx   <- which(terms != "(Intercept)")
    
    tibble(
      variable = var,
      term     = terms[idx],
      log_OR   = round(coef(m)[idx], 3),
      OR       = round(exp(coef(m)[idx]), 3),
      CI_inf   = round(exp(confint(m)[idx, 1]), 3),
      CI_sup   = round(exp(confint(m)[idx, 2]), 3),
      p_valor  = round(m$prob[idx], 4),
      n        = nrow(df_temp),
      n_events = sum(df_temp$overweight_risk_6m == 1),
      selected = ifelse(m$prob[idx] < 0.20, "Sí", "No")
    )
  }, error = function(e) {
    cat("  Error en", var, ":", conditionMessage(e), "\n")
    return(NULL)
  })
  
}) %>%
  bind_rows() %>%
  arrange(p_valor)

kable(univariate_logit,
      caption = "Cribado univariante logístico (Firth). Umbral p < 0,20.",
      col.names = c("Variable", "Término", "log(OR)", "OR",
                     "IC inf", "IC sup", "p-valor", "n",
                     "Eventos", "Seleccionada")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE)

selected_logit <- univariate_logit %>%
  filter(selected == "Sí") %>%
  pull(variable) %>%
  unique()

cat("Variables seleccionadas (p < 0,20):",
    paste(selected_logit, collapse = ", "), "\n")

m_firth_A <- logistf(overweight_risk_6m ~ mat_weight_gain, data = df_cc)

coef_A <- tibble(
  Término = names(coef(m_firth_A)),
  log_OR  = round(coef(m_firth_A), 3),
  OR      = round(exp(coef(m_firth_A)), 3),
  CI_inf  = round(exp(confint(m_firth_A)[, 1]), 3),
  CI_sup  = round(exp(confint(m_firth_A)[, 2]), 3),
  p_valor = round(m_firth_A$prob, 4)
)

kable(coef_A,
      caption = "Modelo A (Firth): overweight_risk_6m ~ GWG.",
      col.names = c("Término", "log(OR)", "OR",
                     "IC inf", "IC sup", "p-valor")) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = FALSE)

m_firth_B <- logistf(overweight_risk_6m ~ mat_weight_gain + inf_sex,
                     data = df_cc)

coef_B <- tibble(
  Término = names(coef(m_firth_B)),
  log_OR  = round(coef(m_firth_B), 3),
  OR      = round(exp(coef(m_firth_B)), 3),
  CI_inf  = round(exp(confint(m_firth_B)[, 1]), 3),
  CI_sup  = round(exp(confint(m_firth_B)[, 2]), 3),
  p_valor = round(m_firth_B$prob, 4)
)

kable(coef_B,
      caption = "Modelo B (Firth): overweight_risk_6m ~ GWG + sexo.",
      col.names = c("Término", "log(OR)", "OR",
                     "IC inf", "IC sup", "p-valor")) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = FALSE)

lr_AB <- anova(m_firth_A, m_firth_B, method = "PLR")
cat("LRT Modelo A vs B (efecto de añadir sexo):\n")
print(lr_AB)

df_cc$pred_prob_B <- predict(m_firth_B, type = "response")

p_pred <- ggplot(df_cc,
                 aes(x = pred_prob_B,
                     fill = factor(overweight_risk_6m,
                                   labels = c("Sin riesgo", "Riesgo")))) +
  geom_histogram(bins = 12, alpha = 0.7, position = "identity",
                 color = "white") +
  scale_fill_manual(values = c("Sin riesgo" = "#2166AC",
                               "Riesgo"     = "#B2182B")) +
  labs(x = "Probabilidad predicha", y = "Frecuencia",
       fill = "Desenlace observado",
       title = "Distribución de probabilidades predichas (Modelo B Firth)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p_pred

m_glm_A <- glm(overweight_risk_6m ~ mat_weight_gain,
               family = binomial(link = "logit"), data = df_cc)

m_glm_B <- glm(overweight_risk_6m ~ mat_weight_gain + inf_sex,
               family = binomial(link = "logit"), data = df_cc)
pred_glm_B <- predict(m_glm_B, type = "response")

roc_obj <- roc(df_cc$overweight_risk_6m,
               pred_glm_B,
               quiet = TRUE)
auc_val <- round(auc(roc_obj), 3)
auc_ci  <- round(ci.auc(roc_obj), 3)

cat("AUC (Modelo B):", auc_val, "\n")
cat("IC 95% (DeLong):", auc_ci[1], "-", auc_ci[3], "\n")

roc_df <- data.frame(
  sensibilidad    = roc_obj$sensitivities,
  especificidad_1 = 1 - roc_obj$specificities
)

p_roc <- ggplot(roc_df, aes(x = especificidad_1, y = sensibilidad)) +
  geom_line(color = "#B2182B", linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50") +
  annotate("text", x = 0.6, y = 0.3,
           label = paste0("AUC = ", auc_val,
                          "\n[", auc_ci[1], " – ", auc_ci[3], "]"),
           size = 4, color = "#B2182B") +
  labs(x = "1 − Especificidad", y = "Sensibilidad",
       title = "Curva ROC — Modelo B (GWG + sexo)") +
  coord_equal() +
  theme_minimal(base_size = 11)

p_roc

set.seed(42)

boot_auc <- function(data, indices) {
  d <- data[indices, ]
  if (sum(d$overweight_risk_6m == 1) < 1 |
      sum(d$overweight_risk_6m == 0) < 1) {
    return(NA_real_)
  }
  m <- glm(overweight_risk_6m ~ mat_weight_gain + inf_sex,
           family = binomial, data = d)
  pred <- predict(m, type = "response")
  tryCatch(
    as.numeric(auc(roc(d$overweight_risk_6m, pred, quiet = TRUE))),
    error = function(e) NA_real_
  )
}

boot_result <- boot(df_cc, boot_auc, R = 2000)

auc_boot_mean <- round(mean(boot_result$t, na.rm = TRUE), 3)
auc_boot_ci   <- round(quantile(boot_result$t, c(0.025, 0.975),
                                na.rm = TRUE), 3)
optimismo     <- round(auc_boot_mean - auc_val, 3)
auc_corregido <- round(auc_val - abs(optimismo), 3)

cat("AUC aparente (modelo):  ", auc_val, "\n")
cat("AUC bootstrap (media):  ", auc_boot_mean, "\n")
cat("IC 95% bootstrap:       ", auc_boot_ci[1], "-", auc_boot_ci[2], "\n")
cat("Optimismo estimado:     ", optimismo, "\n")
cat("AUC corregido:          ", auc_corregido, "\n")

cooksd_B  <- cooks.distance(m_glm_B)
umbral_lr <- 4 / nobs(m_glm_B)
n_infl_lr <- sum(cooksd_B > umbral_lr, na.rm = TRUE)

cat("Umbral Cook's D (4/n):", round(umbral_lr, 4), "\n")
cat("Observaciones influyentes:", n_infl_lr, "\n")



forest_data <- tibble(
  term    = names(coef(m_firth_B)),
  log_OR  = coef(m_firth_B),
  OR      = exp(coef(m_firth_B)),
  CI_inf  = exp(confint(m_firth_B)[, 1]),
  CI_sup  = exp(confint(m_firth_B)[, 2]),
  p_valor = m_firth_B$prob
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    sig = case_when(
      p_valor < 0.05 ~ "p < 0,05",
      p_valor < 0.10 ~ "p < 0,10",
      p_valor < 0.20 ~ "p < 0,20",
      TRUE           ~ "p ≥ 0,20"
    ),
    label = case_when(
      term == "mat_weight_gain" ~ "Ganancia peso (kg)",
      term == "inf_sexMale"     ~ "Sexo: varón",
      TRUE ~ term
    )
  )

p_forest <- ggplot(forest_data,
                   aes(x = OR, y = reorder(label, OR), color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_inf, xmax = CI_sup), height = 0.2) +
  scale_color_manual(
    values = c("p < 0,05" = "#B2182B",
               "p < 0,10" = "#F4A582",
               "p < 0,20" = "#FDB863",
               "p ≥ 0,20" = "#999999"),
    name = "Significación"
  ) +
  scale_x_log10() +
  labs(x = "OR (IC 95%, escala log)", y = NULL,
       title = "Modelo B (Firth): overweight_risk_6m ~ GWG + sexo",
       subtitle = paste0("n = ", n_total_lr, " | Eventos = ", n_events,
                         " | Regresión penalizada de Firth")) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom"
  )

p_forest

m_cond <- lm(WLZ_6 ~ WLZ_0 + mat_weight_gain + inf_sex, data = df)

coherencia <- tibble(
  Predictor      = c("mat_weight_gain", "inf_sexMale"),
  `β_lineal`     = round(coef(m_cond)[c("mat_weight_gain",
                                         "inf_sexMale")], 3),
  Dirección_lin  = ifelse(coef(m_cond)[c("mat_weight_gain",
                                          "inf_sexMale")] > 0,
                          "↑ WLZ₆", "↓ WLZ₆"),
  OR_logístico   = round(exp(coef(m_firth_B)[c("mat_weight_gain",
                                                "inf_sexMale")]), 3),
  Dirección_log  = ifelse(exp(coef(m_firth_B)[c("mat_weight_gain",
                                                 "inf_sexMale")]) > 1,
                          "↑ Riesgo", "↓ Riesgo"),
  Coherente      = ifelse(
    (coef(m_cond)[c("mat_weight_gain", "inf_sexMale")] > 0 &
       exp(coef(m_firth_B)[c("mat_weight_gain", "inf_sexMale")]) > 1) |
      (coef(m_cond)[c("mat_weight_gain", "inf_sexMale")] < 0 &
         exp(coef(m_firth_B)[c("mat_weight_gain", "inf_sexMale")]) < 1),
    "Sí", "No"
  )
)

kable(coherencia,
      caption = "Coherencia entre modelos lineales y logístico.") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = FALSE)

univariate_glm <- lapply(all_candidates, function(var) {
  
  df_temp <- df_lr %>%
    select(overweight_risk_6m, all_of(var)) %>%
    filter(complete.cases(.))
  
  if (length(unique(df_temp[[var]])) < 2) return(NULL)
  
  fml <- as.formula(paste("overweight_risk_6m ~", var))
  
  tryCatch({
    m <- glm(fml, family = binomial(link = "logit"), data = df_temp)
    s <- summary(m)
    
    terms <- names(coef(m))
    idx   <- which(terms != "(Intercept)")
    
    # IC de Wald (por defecto en glm)
    ci <- confint.default(m)
    
    tibble(
      variable = var,
      term     = terms[idx],
      log_OR   = round(coef(m)[idx], 3),
      OR       = round(exp(coef(m)[idx]), 3),
      CI_inf   = round(exp(ci[idx, 1]), 3),
      CI_sup   = round(exp(ci[idx, 2]), 3),
      p_valor  = round(s$coefficients[idx, "Pr(>|z|)"], 4),
      n        = nrow(df_temp),
      n_events = sum(df_temp$overweight_risk_6m == 1),
      selected = ifelse(s$coefficients[idx, "Pr(>|z|)"] < 0.20, "Sí", "No")
    )
  }, error = function(e) {
    cat("  Error en", var, ":", conditionMessage(e), "\n")
    return(NULL)
  })
  
}) %>%
  bind_rows() %>%
  arrange(p_valor)

kable(univariate_glm,
      caption = "Cribado univariante logístico (glm clásico). Umbral p < 0,20.",
      col.names = c("Variable", "Término", "log(OR)", "OR",
                     "IC inf", "IC sup", "p-valor", "n",
                     "Eventos", "Seleccionada")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE)

selected_glm <- univariate_glm %>%
  filter(selected == "Sí") %>%
  pull(variable) %>%
  unique()

cat("Variables seleccionadas (p < 0,20):",
    paste(selected_glm, collapse = ", "), "\n")

m_glm_A <- glm(overweight_risk_6m ~ mat_weight_gain,
               family = binomial(link = "logit"), data = df_cc)

ci_A <- confint.default(m_glm_A)

coef_A_glm <- tibble(
  Término = names(coef(m_glm_A)),
  log_OR  = round(coef(m_glm_A), 3),
  OR      = round(exp(coef(m_glm_A)), 3),
  CI_inf  = round(exp(ci_A[, 1]), 3),
  CI_sup  = round(exp(ci_A[, 2]), 3),
  p_valor = round(summary(m_glm_A)$coefficients[, "Pr(>|z|)"], 4)
)

kable(coef_A_glm,
      caption = "Modelo A (glm clásico): overweight_risk_6m ~ GWG.",
      col.names = c("Término", "log(OR)", "OR",
                     "IC inf", "IC sup", "p-valor")) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = FALSE)

m_glm_B <- glm(overweight_risk_6m ~ mat_weight_gain + inf_sex,
               family = binomial(link = "logit"), data = df_cc)

ci_B <- confint.default(m_glm_B)

coef_B_glm <- tibble(
  Término = names(coef(m_glm_B)),
  log_OR  = round(coef(m_glm_B), 3),
  OR      = round(exp(coef(m_glm_B)), 3),
  CI_inf  = round(exp(ci_B[, 1]), 3),
  CI_sup  = round(exp(ci_B[, 2]), 3),
  p_valor = round(summary(m_glm_B)$coefficients[, "Pr(>|z|)"], 4)
)

kable(coef_B_glm,
      caption = "Modelo B (glm clásico): overweight_risk_6m ~ GWG + sexo.",
      col.names = c("Término", "log(OR)", "OR",
                     "IC inf", "IC sup", "p-valor")) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = FALSE)

# Detectar separación cuasi-completa
if (any(abs(coef(m_glm_B)) > 10)) {
  cat("\n\n> **ADVERTENCIA:** Al menos un coeficiente tiene un valor absoluto",
      "superior a 10, lo que indica separación cuasi-completa.",
      "Las estimaciones del glm clásico son inestables y no deben",
      "interpretarse. Esto justifica el uso de la penalización de Firth.\n\n")
}

lr_AB_glm <- anova(m_glm_A, m_glm_B, test = "LRT")
cat("LRT clásico Modelo A vs B (efecto de añadir sexo):\n")
print(lr_AB_glm)

df_cc$pred_prob_glm_B <- predict(m_glm_B, type = "response")

p_pred_glm <- ggplot(df_cc,
                     aes(x = pred_prob_glm_B,
                         fill = factor(overweight_risk_6m,
                                       labels = c("Sin riesgo", "Riesgo")))) +
  geom_histogram(bins = 12, alpha = 0.7, position = "identity",
                 color = "white") +
  scale_fill_manual(values = c("Sin riesgo" = "#2166AC",
                               "Riesgo"     = "#B2182B")) +
  labs(x = "Probabilidad predicha", y = "Frecuencia",
       fill = "Desenlace observado",
       title = "Distribución de probabilidades predichas (Modelo B glm clásico)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p_pred_glm

cooksd_glm  <- cooks.distance(m_glm_B)
umbral_glm  <- 4 / nobs(m_glm_B)
n_infl_glm  <- sum(cooksd_glm > umbral_glm, na.rm = TRUE)

cat("Umbral Cook's D (4/n):", round(umbral_glm, 4), "\n")
cat("Observaciones influyentes:", n_infl_glm, "\n")



ci_forest_glm <- confint.default(m_glm_B)

forest_glm <- tibble(
  term    = names(coef(m_glm_B)),
  log_OR  = coef(m_glm_B),
  OR      = exp(coef(m_glm_B)),
  CI_inf  = exp(ci_forest_glm[, 1]),
  CI_sup  = exp(ci_forest_glm[, 2]),
  p_valor = summary(m_glm_B)$coefficients[, "Pr(>|z|)"]
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    sig = case_when(
      p_valor < 0.05 ~ "p < 0,05",
      p_valor < 0.10 ~ "p < 0,10",
      p_valor < 0.20 ~ "p < 0,20",
      TRUE           ~ "p ≥ 0,20"
    ),
    label = case_when(
      term == "mat_weight_gain" ~ "Ganancia peso (kg)",
      term == "inf_sexMale"     ~ "Sexo: varón",
      TRUE ~ term
    )
  )

p_forest_glm <- ggplot(forest_glm,
                       aes(x = OR, y = reorder(label, OR), color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_inf, xmax = CI_sup), height = 0.2) +
  scale_color_manual(
    values = c("p < 0,05" = "#B2182B",
               "p < 0,10" = "#F4A582",
               "p < 0,20" = "#FDB863",
               "p ≥ 0,20" = "#999999"),
    name = "Significación"
  ) +
  scale_x_log10() +
  labs(x = "OR (IC 95% Wald, escala log)", y = NULL,
       title = "Modelo B (glm clásico): overweight_risk_6m ~ GWG + sexo",
       subtitle = paste0("n = ", nobs(m_glm_B), " | Eventos = ",
                         sum(df_cc$overweight_risk_6m == 1),
                         " | Regresión logística clásica (sin penalización)")) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom"
  )

p_forest_glm

coherencia_glm <- tibble(
  Predictor      = c("mat_weight_gain", "inf_sexMale"),
  `β_lineal`     = round(coef(m_cond)[c("mat_weight_gain",
                                         "inf_sexMale")], 3),
  Dirección_lin  = ifelse(coef(m_cond)[c("mat_weight_gain",
                                          "inf_sexMale")] > 0,
                          "↑ WLZ₆", "↓ WLZ₆"),
  OR_logístico   = round(exp(coef(m_glm_B)[c("mat_weight_gain",
                                              "inf_sexMale")]), 3),
  Dirección_log  = ifelse(exp(coef(m_glm_B)[c("mat_weight_gain",
                                               "inf_sexMale")]) > 1,
                          "↑ Riesgo", "↓ Riesgo"),
  Coherente      = ifelse(
    (coef(m_cond)[c("mat_weight_gain", "inf_sexMale")] > 0 &
       exp(coef(m_glm_B)[c("mat_weight_gain", "inf_sexMale")]) > 1) |
      (coef(m_cond)[c("mat_weight_gain", "inf_sexMale")] < 0 &
         exp(coef(m_glm_B)[c("mat_weight_gain", "inf_sexMale")]) < 1),
    "Sí", "No"
  )
)

kable(coherencia_glm,
      caption = "Coherencia entre modelos lineales y logístico (glm clásico).") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = FALSE)

# Unir resultados de ambos cribados
comp_cribado <- univariate_logit %>%
  select(variable, term, OR_firth = OR, p_firth = p_valor) %>%
  left_join(
    univariate_glm %>%
      select(variable, term, OR_glm = OR, p_glm = p_valor),
    by = c("variable", "term")
  ) %>%
  mutate(
    delta_OR_pct = round((OR_glm - OR_firth) / OR_firth * 100, 1),
    delta_p      = round(p_glm - p_firth, 4)
  )

kable(comp_cribado,
      caption = "Cribado univariante: comparación Firth vs glm clásico.",
      col.names = c("Variable", "Término",
                     "OR Firth", "p Firth",
                     "OR glm", "p glm",
                     "Δ OR (%)", "Δ p")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE)

ci_firth_B <- confint(m_firth_B)
ci_glm_B   <- confint.default(m_glm_B)

comp_modelo_B <- tibble(
  Término = names(coef(m_firth_B)),
  
  OR_Firth     = round(exp(coef(m_firth_B)), 3),
  CI_Firth     = paste0("[", round(exp(ci_firth_B[, 1]), 3), " – ",
                        round(exp(ci_firth_B[, 2]), 3), "]"),
  p_Firth      = round(m_firth_B$prob, 4),
  
  OR_glm       = round(exp(coef(m_glm_B)), 3),
  CI_glm       = paste0("[", round(exp(ci_glm_B[, 1]), 3), " – ",
                        round(exp(ci_glm_B[, 2]), 3), "]"),
  p_glm        = round(summary(m_glm_B)$coefficients[, "Pr(>|z|)"], 4),
  
  `Δ_OR_%`     = round((exp(coef(m_glm_B)) - exp(coef(m_firth_B))) /
                          exp(coef(m_firth_B)) * 100, 1)
)

kable(comp_modelo_B,
      caption = "Modelo B: comparación Firth vs glm clásico.",
      col.names = c("Término",
                     "OR Firth", "IC 95% Firth", "p Firth",
                     "OR glm", "IC 95% glm", "p glm",
                     "Δ OR (%)")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE)

# Preparar datos de Firth
forest_firth <- tibble(
  term    = names(coef(m_firth_B)),
  OR      = exp(coef(m_firth_B)),
  CI_inf  = exp(confint(m_firth_B)[, 1]),
  CI_sup  = exp(confint(m_firth_B)[, 2]),
  Método  = "Firth"
) %>% filter(term != "(Intercept)")

# Preparar datos de glm
forest_glm_comp <- tibble(
  term    = names(coef(m_glm_B)),
  OR      = exp(coef(m_glm_B)),
  CI_inf  = exp(ci_glm_B[, 1]),
  CI_sup  = exp(ci_glm_B[, 2]),
  Método  = "glm clásico"
) %>% filter(term != "(Intercept)")

# Unir y formatear
forest_comp <- bind_rows(forest_firth, forest_glm_comp) %>%
  mutate(
    label = case_when(
      term == "mat_weight_gain" ~ "Ganancia peso (kg)",
      term == "inf_sexMale"     ~ "Sexo: varón",
      TRUE ~ term
    ),
    # Truncar IC extremos para visualización
    CI_inf = pmax(CI_inf, 0.0001),
    CI_sup = pmin(CI_sup, 100)
  )

p_comp <- ggplot(forest_comp,
                 aes(x = OR, y = label, color = Método,
                     shape = Método)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = CI_inf, xmax = CI_sup),
                 height = 0.2,
                 position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Firth"       = "#B2182B",
                                "glm clásico" = "#2166AC")) +
  scale_shape_manual(values = c("Firth" = 16, "glm clásico" = 17)) +
  scale_x_log10() +
  labs(x = "OR (IC 95%, escala log)", y = NULL,
       title = "Comparación de OR: Firth vs glm clásico (Modelo B)",
       subtitle = "La penalización de Firth estabiliza las estimaciones hacia valores finitos") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom"
  )

p_comp

# Tablas
write.csv(univariate_logit,
          "../R_output/tables/univariate_screening_logistic.csv",
          row.names = FALSE)
write.csv(coef_A,
          "../R_output/tables/logistic_model_A_firth.csv",
          row.names = FALSE)
write.csv(coef_B,
          "../R_output/tables/logistic_model_B_firth.csv",
          row.names = FALSE)

if (exists("coef_C")) {
  write.csv(coef_C,
            "../R_output/tables/logistic_model_C_firth.csv",
            row.names = FALSE)
}

# Figuras
ggsave("../R_output/figures/fig_boxplots_logistica.pdf", p_box,
       width = 14, height = 8)
ggsave("../R_output/figures/fig_pred_prob_logistica.pdf", p_pred,
       width = 8, height = 5)
ggsave("../R_output/figures/fig_roc_logistica.pdf", p_roc,
       width = 6, height = 6)
ggsave("../R_output/figures/fig_forest_logistica.pdf", p_forest,
       width = 8, height = 4)

if (exists("p_forest_C")) {
  ggsave("../R_output/figures/fig_forest_logistica_C.pdf", p_forest_C,
         width = 8, height = 5)
}

# Modelos
saveRDS(m_firth_A, "../R_output/models/logistic_firth_A.rds")
saveRDS(m_firth_B, "../R_output/models/logistic_firth_B.rds")
if (exists("m_firth_C")) {
  saveRDS(m_firth_C, "../R_output/models/logistic_firth_C.rds")
}

# Resumen ejecutivo
resumen_logistic <- list(
  n_total     = n_total_lr,
  n_events    = n_events,
  prevalence  = prevalence,
  auc_B       = auc_val,
  auc_ci      = auc_ci,
  auc_boot_ci = auc_boot_ci,
  auc_corr    = auc_corregido,
  OR_gwg_B    = round(exp(coef(m_firth_B)["mat_weight_gain"]), 3),
  OR_sex_B    = round(exp(coef(m_firth_B)["inf_sexMale"]), 3)
)
saveRDS(resumen_logistic, "../R_output/models/logistic_summary.rds")

cat("Tablas guardadas en output/tables/\n")
cat("Figuras guardadas en output/figures/\n")
cat("Modelos guardados en output/models/\n")
