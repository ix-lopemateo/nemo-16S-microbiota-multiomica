#!/usr/bin/env Rscript
# Fase 1 clinica. Modelo multivariante; genera los residuos de adiposidad

library(tidyverse)
library(broom)
library(car)
library(knitr)
library(kableExtra)
library(patchwork)

df  <- readRDS("../R_output/models/df_analytic_clinico.rds")
sel <- readRDS("../R_output/models/variable_selection_clinico.rds")

cat("Dataset cargado")

cat("Variables seleccionadas (Modelo WLZ6 ajustado por WLZ0):", paste(sel$selected_cond, collapse = ", "), "\n")
cat("Variables seleccionadas (Modelo dWLZ_6_0):",     paste(sel$selected_direct, collapse = ", "), "\n")
cat("Variables forzadas:",                     paste(sel$forced, collapse = ", "), "\n")

m_cond <- lm(WLZ_6 ~ WLZ_0 + mat_weight_gain + inf_sex + feeding_6m, data = df)

tidy(m_cond, conf.int = TRUE) %>%
  mutate(across(c(estimate, conf.low, conf.high, statistic), ~round(., 3)),
         p.value = round(p.value, 4)) %>%
  kable(caption = "Modelo condicional: WLZ_6 ~ WLZ_0 + GWG + sexo + alimentación.",
        col.names = c("Término", "β", "SE", "t", "p-valor",
                       "IC inf", "IC sup")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE)

g_cond <- glance(m_cond)
cat("n =", nobs(m_cond),
    "| R² aj. =", round(g_cond$adj.r.squared, 4),
    "| F =", round(g_cond$statistic, 2),
    "| p(F) =", format.pval(g_cond$p.value, digits = 4), "\n")

vif_vals <- vif(m_cond)
vif_df <- as.data.frame(vif_vals)

if (ncol(vif_df) > 1) {
  vif_df$Variable <- rownames(vif_df)
  vif_df <- vif_df %>% select(Variable, everything())
}

kable(vif_df, digits = 3,
      caption = "VIF/GVIF del modelo condicional.") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)

resid_df <- data.frame(
  fitted    = fitted(m_cond),
  residuals = residuals(m_cond),
  std_resid = rstandard(m_cond)
)

p_qq <- ggplot(resid_df, aes(sample = std_resid)) +
  stat_qq() + stat_qq_line(color = "#B2182B") +
  labs(title = "QQ-plot", x = "Cuantiles teóricos", y = "Cuantiles observados") +
  theme_minimal(base_size = 11)

p_hist <- ggplot(resid_df, aes(x = residuals)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 20, fill = "#2166AC", color = "white", alpha = 0.6) +
  geom_density(color = "#B2182B", linewidth = 1) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(resid_df$residuals), 
                            sd = sd(resid_df$residuals)),
                color = "black", linetype = "dashed", linewidth = 0.8) +
  labs(title = "Histograma de residuos", 
       x = "Residuos", y = "Densidad",
       caption = "Línea roja: densidad empírica · Línea negra: normal teórica") +
  theme_minimal(base_size = 11)

p_qq + p_hist

shapiro_p <- shapiro.test(residuals(m_cond))$p.value
cat("Shapiro-Wilk p =", round(shapiro_p, 4), "\n")

ggplot(resid_df, aes(x = fitted, y = std_resid)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "loess", se = FALSE, color = "#B2182B") +
  labs(title = "Residuos vs predichos",
       x = "Valores predichos", y = "Residuos estandarizados") +
  theme_minimal(base_size = 11)

bp_p <- car::ncvTest(m_cond)$p
cat("Breusch-Pagan p =", round(bp_p, 4), "\n")

cooks_d <- cooks.distance(m_cond)
umbral  <- 4 / nobs(m_cond)

infl_df <- data.frame(obs = seq_along(cooks_d), cooks_d = cooks_d) %>%
  mutate(influyente = cooks_d > umbral)

ggplot(infl_df, aes(x = obs, y = cooks_d, color = influyente)) +
  geom_point(size = 2) +
  geom_hline(yintercept = umbral, linetype = "dashed", color = "#B2182B") +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#B2182B"),
                     labels = c("Normal", "Influyente")) +
  labs(title = "Distancia de Cook",
       subtitle = paste0("Umbral: 4/n = ", round(umbral, 3)),
       x = "Observación", y = "Cook's D", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

n_infl <- sum(infl_df$influyente)
cat("Observaciones influyentes:", n_infl, "\n")

df_cc <- df[complete.cases(df[, all.vars(formula(m_cond))]), ]

# Inspección de las observaciones influyentes
obs_infl <- which(infl_df$influyente)

infl_inspect <- data.frame(
  obs        = obs_infl,
  dyad_id    = df_cc$dyad_id[obs_infl],
  cooks_d    = round(cooks_d[obs_infl], 3),
  residuo    = round(residuals(m_cond)[obs_infl], 3),
  std_resid  = round(rstandard(m_cond)[obs_infl], 3),
  leverage   = round(hatvalues(m_cond)[obs_infl], 3),
  WLZ_6      = df_cc$WLZ_6[obs_infl],
  WLZ_0      = df_cc$WLZ_0[obs_infl],
  GWG        = df_cc$mat_weight_gain[obs_infl],
  sexo       = df_cc$inf_sex[obs_infl],
  feeding    = df_cc$feeding_6m[obs_infl]
)

kable(infl_inspect,
      caption = "Inspección de las observaciones potencialmente influyentes.") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)

crPlots(m_cond, terms = "mat_weight_gain",
        main = "Residuos parciales: ganancia de peso gestacional",
        col = "#2166AC", col.lines = "#B2182B")

# Filtrar solo observaciones influyentes CON residuo grande (|std_resid| > 2)
obs_infl_real <- which(infl_df$influyente & abs(rstandard(m_cond)) > 2)
n_infl_real <- length(obs_infl_real)

cat("Observaciones influyentes con residuo grande:", n_infl_real, "\n")

if (n_infl_real > 0 & n_infl_real <= 5) {
  df_cc <- df[complete.cases(df[, all.vars(formula(m_cond))]), ]
  m_sens <- update(m_cond, data = df_cc[-obs_infl_real, ])

  coef_comp <- left_join(
    tidy(m_cond) %>% select(term, est_orig = estimate),
    tidy(m_sens) %>% select(term, est_sens = estimate),
    by = "term"
  ) %>%
    mutate(`Δ%` = round((est_sens - est_orig) / abs(est_orig) * 100, 1))

  kable(coef_comp, digits = 3,
        caption = "Coeficientes con y sin outliers reales (|std_resid| > 2).") %>%
    kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)
} else {
  cat("No se detectaron outliers reales con residuo grande.\n")
}

forest_df <- tidy(m_cond, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    sig = case_when(p.value < 0.05 ~ "p < 0.05",
                    p.value < 0.10 ~ "p < 0.10",
                    TRUE           ~ "p ≥ 0.10"),
    label = case_when(
      term == "WLZ_0"                         ~ "WLZ nacimiento",
      term == "mat_weight_gain"                ~ "Ganancia peso (kg)",
      term == "inf_sexMale"                    ~ "Sexo: varón",
      term == "feeding_6mBF_complementary"     ~ "Alim: LM + complementaria",
      term == "feeding_6mExclusive_Formula"    ~ "Alim: artificial exclusiva",
      TRUE ~ term
    )
  )

ggplot(forest_df, aes(x = estimate, y = reorder(label, estimate), color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  scale_color_manual(values = c("p < 0.05" = "#B2182B",
                                "p < 0.10" = "#F4A582",
                                "p ≥ 0.10" = "#999999"),
                     name = "Significación") +
  labs(x = expression(beta ~ "(IC 95%)"), y = NULL,
       title = "Modelo condicional: WLZ_6 ~ WLZ_0 + GWG + sexo + alimentación",
       subtitle = paste0("n = ", nobs(m_cond),
                         " | R² aj. = ", round(g_cond$adj.r.squared, 4))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_line(color = "grey92"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

m_dir <- lm(dWLZ_0_6 ~ mat_weight_gain + inf_sex + feeding_6m, data = df)

tidy(m_dir, conf.int = TRUE) %>%
  mutate(across(c(estimate, conf.low, conf.high, statistic), ~round(., 3)),
         p.value = round(p.value, 4)) %>%
  kable(caption = "Modelo directo: dWLZ_0_6 ~ GWG + sexo + alimentación.",
        col.names = c("Término", "β", "SE", "t", "p-valor",
                       "IC inf", "IC sup")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE)

g_dir <- glance(m_dir)
cat("n =", nobs(m_dir),
    "| R² aj. =", round(g_dir$adj.r.squared, 4),
    "| F =", round(g_dir$statistic, 2),
    "| p(F) =", format.pval(g_dir$p.value, digits = 4), "\n")

vif_vals_dir <- vif(m_dir)
vif_df_dir <- as.data.frame(vif_vals_dir)

if (ncol(vif_df_dir) > 1) {
  vif_df_dir$Variable <- rownames(vif_df_dir)
  vif_df_dir <- vif_df_dir %>% select(Variable, everything())
}

kable(vif_df_dir, digits = 3,
      caption = "VIF/GVIF del modelo directo.") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)

resid_dir <- data.frame(
  fitted    = fitted(m_dir),
  residuals = residuals(m_dir),
  std_resid = rstandard(m_dir)
)

p1 <- ggplot(resid_dir, aes(sample = std_resid)) +
  stat_qq() + stat_qq_line(color = "#B2182B") +
  labs(title = "QQ-plot (modelo directo)",
       x = "Cuantiles teóricos", y = "Cuantiles observados") +
  theme_minimal(base_size = 11)

p2 <- ggplot(resid_dir, aes(x = residuals)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 20, fill = "#2166AC", color = "white", alpha = 0.6) +
  geom_density(color = "#B2182B", linewidth = 1) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(resid_dir$residuals), 
                            sd = sd(resid_dir$residuals)),
                color = "black", linetype = "dashed", linewidth = 0.8) +
  labs(title = "Histograma de residuos (modelo directo)", 
       x = "Residuos", y = "Densidad",
       caption = "Línea roja: densidad empírica · Línea negra: normal teórica") +
  theme_minimal(base_size = 11)

p1 + p2

shapiro_p_dir <- shapiro.test(residuals(m_dir))$p.value
cat("Shapiro-Wilk p =", round(shapiro_p_dir, 4), "\n")

ggplot(resid_dir, aes(x = fitted, y = std_resid)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "loess", se = FALSE, color = "#B2182B") +
  labs(title = "Residuos vs predichos (modelo directo)",
       x = "Valores predichos", y = "Residuos estandarizados") +
  theme_minimal(base_size = 11)

bp_p_dir <- car::ncvTest(m_dir)$p
cat("Breusch-Pagan p =", round(bp_p_dir, 4), "\n")

cooks_d_dir <- cooks.distance(m_dir)
umbral_dir  <- 4 / nobs(m_dir)

infl_df_dir <- data.frame(obs = seq_along(cooks_d_dir), cooks_d = cooks_d_dir) %>%
  mutate(influyente = cooks_d > umbral_dir)

ggplot(infl_df_dir, aes(x = obs, y = cooks_d, color = influyente)) +
  geom_point(size = 2) +
  geom_hline(yintercept = umbral_dir, linetype = "dashed", color = "#B2182B") +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#B2182B"),
                     labels = c("Normal", "Influyente")) +
  labs(title = "Distancia de Cook (modelo directo)",
       subtitle = paste0("Umbral: 4/n = ", round(umbral_dir, 3)),
       x = "Observación", y = "Cook's D", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

n_infl_dir <- sum(infl_df_dir$influyente)
cat("Observaciones influyentes:", n_infl_dir, "\n")

crPlots(m_dir, terms = "mat_weight_gain",
        main = "Residuos parciales: GWG (modelo directo)",
        col = "#2166AC", col.lines = "#B2182B")

if (n_infl_dir > 0 & n_infl_dir <= 5) {
  obs_infl_dir <- which(infl_df_dir$influyente)
  df_cc_dir <- df[complete.cases(df[, all.vars(formula(m_dir))]), ]
  m_sens_dir <- update(m_dir, data = df_cc_dir[-obs_infl_dir, ])

  coef_comp_dir <- left_join(
    tidy(m_dir)      %>% select(term, est_orig = estimate),
    tidy(m_sens_dir) %>% select(term, est_sens = estimate),
    by = "term"
  ) %>%
    mutate(`Δ%` = round((est_sens - est_orig) / abs(est_orig) * 100, 1))

  kable(coef_comp_dir, digits = 3,
        caption = "Coeficientes con y sin observaciones influyentes (modelo directo).") %>%
    kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)
} else {
  cat("No se detectaron observaciones influyentes (o más de 5).\n")
}

forest_dir <- tidy(m_dir, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    sig = case_when(p.value < 0.05 ~ "p < 0.05",
                    p.value < 0.10 ~ "p < 0.10",
                    TRUE           ~ "p ≥ 0.10"),
    label = case_when(
      term == "mat_weight_gain"                ~ "Ganancia peso (kg)",
      term == "inf_sexMale"                    ~ "Sexo: varón",
      term == "feeding_6mBF_complementary"     ~ "Alim: LM + complementaria",
      term == "feeding_6mExclusive_Formula"    ~ "Alim: artificial exclusiva",
      TRUE ~ term
    )
  )

ggplot(forest_dir, aes(x = estimate, y = reorder(label, estimate), color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  scale_color_manual(values = c("p < 0.05" = "#B2182B",
                                "p < 0.10" = "#F4A582",
                                "p ≥ 0.10" = "#999999"),
                     name = "Significación") +
  labs(x = expression(beta ~ "(IC 95%)"), y = NULL,
       title = "Modelo directo: dWLZ_0_6 ~ GWG + sexo + alimentación",
       subtitle = paste0("n = ", nobs(m_dir),
                         " | R² aj. = ", round(g_dir$adj.r.squared, 4))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_line(color = "grey92"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

comp_df <- bind_rows(
  tidy(m_cond, conf.int = TRUE) %>%
    filter(term != "(Intercept)" & term != "WLZ_0") %>%
    mutate(Modelo = "Condicional (WLZ_6 | WLZ_0)"),
  tidy(m_dir, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(Modelo = "Directo (dWLZ_0_6)")
) %>%
  mutate(
    label = case_when(
      term == "mat_weight_gain"                ~ "GWG",
      term == "inf_sexMale"                    ~ "Sexo: varón",
      term == "feeding_6mBF_complementary"     ~ "LM + complementaria",
      term == "feeding_6mExclusive_Formula"    ~ "Artificial exclusiva",
      TRUE ~ term
    )
  )

ggplot(comp_df, aes(x = estimate, y = label, color = Modelo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  labs(x = expression(beta ~ "(IC 95%)"), y = NULL,
       title = "Comparación de coeficientes: condicional vs directo") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

dir.create("../R_output/models", recursive = TRUE, showWarnings = FALSE)

# Modelos
saveRDS(m_cond,  "../R_output/models/m_condicional_clinico.rds")
saveRDS(m_dir,   "../R_output/models/m_directo_clinico.rds")

# Tablas
write.csv(tidy(m_cond, conf.int = TRUE),
          "../R_output/tables/coef_modelo_condicional_clinico.csv", row.names = FALSE)
write.csv(tidy(m_dir, conf.int = TRUE),
          "../R_output/tables/coef_modelo_directo_clinico.csv", row.names = FALSE)

# Residuos del modelo CLÍNICO de la Fase 1 (n=29, sin ómicas). Son los que
# consume el bloque ómico (scripts 09/09b/11/13b/13c y el reporte): se exportan
# con el MISMO nombre y ruta desde los que esos scripts los leen
# (R_output/models/residuos_clinicos_fase1.rds), de modo que A2 es su
# generador local y el proyecto queda autocontenido.
df_resid <- df_cc %>%
  mutate(residuos  = residuals(m_cond),
         predichos = fitted(m_cond)) %>%
  select(dyad_id, residuos, predichos)
write.csv(df_resid, "../R_output/tables/residuos_modelo_condicional_clinico.csv",
          row.names = FALSE)
dir.create(file.path("..", "R_output", "models"), showWarnings = FALSE, recursive = TRUE)
saveRDS(df_resid, file.path("..", "R_output", "models", "residuos_clinicos_fase1.rds"))

# Figura
ggsave("../R_output/figures/fig_forest_condicional_clinico.pdf",
       last_plot(), width = 8, height = 5, dpi = 300)

cat("Exportación completada.\n")
