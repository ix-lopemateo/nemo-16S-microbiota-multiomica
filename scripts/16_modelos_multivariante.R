#!/usr/bin/env Rscript
# Modelo multivariante de la adiposidad infantil (panel completo y nucleo clinico)

library(tidyverse)
library(broom)
library(car)
library(knitr)
library(kableExtra)
library(patchwork)
library(sandwich)   # vcovHC: matriz de covarianzas robusta a heterocedasticidad
library(lmtest)     # coeftest: contraste de coeficientes con SE robustos

df  <- readRDS("../R_output/models/df_analytic_omico.rds")

# Panel completo de variables significativas en el univariante (sin reducir),
# excluyendo CRH (suelo del ensayo), el módulo PICRUSt2 (función inferida) y la
# diversidad alfa (colineal).
panel_16b <- c("mat_weight_gain", "feeding_6m", "met_b6m_BCAT1", "inf_sex")

# Condicional: WLZ_0 se fuerza como línea base (el desenlace es WLZ_6 ajustado por
# WLZ_0). Directo: el desenlace dWLZ_0_6 ya incorpora WLZ_0, por lo que NO se añade.
form_cond <- reformulate(c("WLZ_0", panel_16b), response = "WLZ_6")
form_dir  <- reformulate(panel_16b,             response = "dWLZ_0_6")

cat("Modelo condicional:", paste(deparse(form_cond), collapse = " "), "\n")
cat("Modelo directo:    ", paste(deparse(form_dir),  collapse = " "), "\n")

m_cond <- lm(form_cond, data = df)

tidy(m_cond, conf.int = TRUE) %>%
  mutate(across(c(estimate, conf.low, conf.high, statistic), ~round(., 3)),
         p.value = round(p.value, 4)) %>%
  kable(caption = "Modelo condicional (panel univariante completo).",
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
      caption = "VIF/GVIF del modelo condicional (panel completo).") %>%
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

resid_vals <- resid_df$residuals
n_resid    <- length(resid_vals)
nbin       <- 20
brks       <- seq(min(resid_vals), max(resid_vals), length.out = nbin + 1)
bw         <- brks[2] - brks[1]
hist_c     <- hist(resid_vals, breaks = brks, include.lowest = TRUE, plot = FALSE)
hist_df    <- data.frame(mid = hist_c$mids, count = hist_c$counts)

p_hist <- ggplot() +
  geom_col(data = hist_df, aes(x = mid, y = count),
           width = bw, fill = "#2166AC", color = "white", alpha = 0.6) +
  geom_text(data = subset(hist_df, count > 0),
            aes(x = mid, y = count, label = count),
            vjust = -0.4, size = 3, color = "grey20") +
  geom_rug(data = resid_df, aes(x = residuals),
           sides = "b", alpha = 0.8, color = "#B2182B") +
  stat_function(fun = function(x)
                  dnorm(x, mean(resid_vals), sd(resid_vals)) * n_resid * bw,
                color = "black", linetype = "dashed", linewidth = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Histograma de residuos", x = "Residuos", y = "Recuento de residuos") +
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

robust_cond <- coeftest(m_cond, vcov = vcovHC(m_cond, type = "HC3"))

robust_tab_cond <- data.frame(
  term   = rownames(robust_cond),
  se_HC3 = robust_cond[, "Std. Error"],
  p_HC3  = robust_cond[, "Pr(>|t|)"],
  row.names = NULL
)

hc3_comp_cond <- tidy(m_cond) %>%
  select(term, estimate, se_OLS = std.error, p_OLS = p.value) %>%
  left_join(robust_tab_cond, by = "term") %>%
  filter(term != "(Intercept)") %>%
  mutate(across(c(estimate, se_OLS, se_HC3), ~round(., 3)),
         across(c(p_OLS, p_HC3), ~round(., 4)))

kable(hc3_comp_cond,
      caption = "Modelo condicional: errores clásicos (OLS) frente a robustos (HC3).",
      col.names = c("Término", "β", "SE OLS", "p OLS", "SE HC3", "p HC3")) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)

p_hc3_cond <- setNames(robust_cond[, "Pr(>|t|)"], rownames(robust_cond))

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

crPlots(m_cond, terms = "mat_weight_gain",
        main = "Residuos parciales: ganancia de peso gestacional",
        col = "#2166AC", col.lines = "#B2182B")

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
        caption = "Coeficientes con y sin el outlier real (|std_resid| > 2).") %>%
    kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)
} else {
  cat("No se detectaron outliers reales con residuo grande.\n")
}

etiquetas <- function(term) dplyr::case_when(
  term == "WLZ_0"                       ~ "WLZ nacimiento",
  term == "mat_weight_gain"             ~ "Ganancia peso (kg)",
  term == "feeding_6mBF_complementary"  ~ "Alim: LM + complementaria",
  term == "feeding_6mExclusive_Formula" ~ "Alim: artificial exclusiva",
  term == "met_b6m_BCAT1"               ~ "Metilación BCAT1 (bebé 6m)",
  term == "met_b1m_NRF1"                ~ "Metilación NRF1 (bebé 1m)",
  term == "inf_sexMale"                 ~ "Sexo: varón",
  TRUE ~ term
)

forest_df <- tidy(m_cond, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(p.value < 0.05 ~ "p < 0.05",
                         p.value < 0.10 ~ "p < 0.10",
                         TRUE           ~ "p ≥ 0.10"),
         label = etiquetas(term))

p_forest_cond <- ggplot(forest_df, aes(x = estimate, y = reorder(label, estimate), color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  scale_color_manual(values = c("p < 0.05" = "#B2182B",
                                "p < 0.10" = "#F4A582",
                                "p ≥ 0.10" = "#999999"),
                     name = "Significación") +
  labs(x = expression(beta ~ "(IC 95%)"), y = NULL,
       title = "Modelo condicional (panel univariante completo)",
       subtitle = paste0("n = ", nobs(m_cond),
                         " | R² aj. = ", round(g_cond$adj.r.squared, 4))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_line(color = "grey92"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
p_forest_cond

m_dir <- lm(form_dir, data = df)

tidy(m_dir, conf.int = TRUE) %>%
  mutate(across(c(estimate, conf.low, conf.high, statistic), ~round(., 3)),
         p.value = round(p.value, 4)) %>%
  kable(caption = "Modelo directo (panel univariante completo).",
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
      caption = "VIF/GVIF del modelo directo (panel completo).") %>%
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

resid_vals_d <- resid_dir$residuals
n_resid_d    <- length(resid_vals_d)
brks_d       <- seq(min(resid_vals_d), max(resid_vals_d), length.out = 21)
bw_d         <- brks_d[2] - brks_d[1]
hist_d       <- hist(resid_vals_d, breaks = brks_d, include.lowest = TRUE, plot = FALSE)
hist_df_d    <- data.frame(mid = hist_d$mids, count = hist_d$counts)

p2 <- ggplot() +
  geom_col(data = hist_df_d, aes(x = mid, y = count),
           width = bw_d, fill = "#2166AC", color = "white", alpha = 0.6) +
  geom_text(data = subset(hist_df_d, count > 0),
            aes(x = mid, y = count, label = count),
            vjust = -0.4, size = 3, color = "grey20") +
  geom_rug(data = resid_dir, aes(x = residuals),
           sides = "b", alpha = 0.8, color = "#B2182B") +
  stat_function(fun = function(x)
                  dnorm(x, mean(resid_vals_d), sd(resid_vals_d)) * n_resid_d * bw_d,
                color = "black", linetype = "dashed", linewidth = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Histograma de residuos (modelo directo)",
       x = "Residuos", y = "Recuento de residuos") +
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

robust_dir <- coeftest(m_dir, vcov = vcovHC(m_dir, type = "HC3"))

robust_tab_dir <- data.frame(
  term   = rownames(robust_dir),
  se_HC3 = robust_dir[, "Std. Error"],
  p_HC3  = robust_dir[, "Pr(>|t|)"],
  row.names = NULL
)

hc3_comp_dir <- tidy(m_dir) %>%
  select(term, estimate, se_OLS = std.error, p_OLS = p.value) %>%
  left_join(robust_tab_dir, by = "term") %>%
  filter(term != "(Intercept)") %>%
  mutate(across(c(estimate, se_OLS, se_HC3), ~round(., 3)),
         across(c(p_OLS, p_HC3), ~round(., 4)))

kable(hc3_comp_dir,
      caption = "Modelo directo: errores clásicos (OLS) frente a robustos (HC3).",
      col.names = c("Término", "β", "SE OLS", "p OLS", "SE HC3", "p HC3")) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)

p_hc3_dir <- setNames(robust_dir[, "Pr(>|t|)"], rownames(robust_dir))

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

obs_infl_real_dir <- which(infl_df_dir$influyente & abs(rstandard(m_dir)) > 2)
n_infl_real_dir <- length(obs_infl_real_dir)

cat("Observaciones influyentes con residuo grande:", n_infl_real_dir, "\n")

if (n_infl_real_dir > 0 & n_infl_real_dir <= 5) {
  df_cc_dir <- df[complete.cases(df[, all.vars(formula(m_dir))]), ]
  m_sens_dir <- update(m_dir, data = df_cc_dir[-obs_infl_real_dir, ])

  coef_comp_dir <- left_join(
    tidy(m_dir)      %>% select(term, est_orig = estimate),
    tidy(m_sens_dir) %>% select(term, est_sens = estimate),
    by = "term"
  ) %>%
    mutate(`Δ%` = round((est_sens - est_orig) / abs(est_orig) * 100, 1))

  kable(coef_comp_dir, digits = 3,
        caption = "Coeficientes con y sin el outlier real (modelo directo).") %>%
    kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)
} else {
  cat("No se detectaron outliers reales con residuo grande.\n")
}

forest_dir <- tidy(m_dir, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(p.value < 0.05 ~ "p < 0.05",
                         p.value < 0.10 ~ "p < 0.10",
                         TRUE           ~ "p ≥ 0.10"),
         label = etiquetas(term))

ggplot(forest_dir, aes(x = estimate, y = reorder(label, estimate), color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  scale_color_manual(values = c("p < 0.05" = "#B2182B",
                                "p < 0.10" = "#F4A582",
                                "p ≥ 0.10" = "#999999"),
                     name = "Significación") +
  labs(x = expression(beta ~ "(IC 95%)"), y = NULL,
       title = "Modelo directo (panel univariante completo)",
       subtitle = paste0("n = ", nobs(m_dir),
                         " | R² aj. = ", round(g_dir$adj.r.squared, 4))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_line(color = "grey92"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

panel_red <- c("mat_weight_gain", "inf_sex")
form_cond_red <- reformulate(c("WLZ_0", panel_red), response = "WLZ_6")
form_dir_red  <- reformulate(panel_red,             response = "dWLZ_0_6")
m_cond_red <- lm(form_cond_red, data = df)
m_dir_red  <- lm(form_dir_red,  data = df)
cat("Reducido condicional:", paste(deparse(form_cond_red), collapse = " "), "\n")
cat("Reducido directo:    ", paste(deparse(form_dir_red),  collapse = " "), "\n")

tidy(m_cond_red, conf.int = TRUE) %>%
  mutate(across(c(estimate, conf.low, conf.high, statistic), ~round(., 3)),
         p.value = round(p.value, 4)) %>%
  kable(caption = "Modelo condicional REDUCIDO (núcleo clínico).",
        col.names = c("Término", "β", "SE", "t", "p-valor", "IC inf", "IC sup")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)

tidy(m_dir_red, conf.int = TRUE) %>%
  mutate(across(c(estimate, conf.low, conf.high, statistic), ~round(., 3)),
         p.value = round(p.value, 4)) %>%
  kable(caption = "Modelo directo REDUCIDO (núcleo clínico).",
        col.names = c("Término", "β", "SE", "t", "p-valor", "IC inf", "IC sup")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)

gc_red <- glance(m_cond_red); gd_red <- glance(m_dir_red)
cat("Condicional reducido: n =", nobs(m_cond_red), "| R² aj. =", round(gc_red$adj.r.squared, 3),
    "| p(F) =", format.pval(gc_red$p.value, digits = 3), "\n")
cat("Directo reducido:     n =", nobs(m_dir_red), "| R² aj. =", round(gd_red$adj.r.squared, 3),
    "| p(F) =", format.pval(gd_red$p.value, digits = 3), "\n")

data.frame(
  Predictor = names(vif(m_cond_red)),
  `VIF condicional` = round(as.numeric(vif(m_cond_red)), 3),
  check.names = FALSE
) %>%
  kable(caption = "VIF del modelo condicional reducido.") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)

hc3_tab <- function(m) {
  rob <- coeftest(m, vcov = vcovHC(m, type = "HC3"))
  tidy(m) %>% select(term, estimate, p_OLS = p.value) %>%
    left_join(data.frame(term = rownames(rob), p_HC3 = rob[, "Pr(>|t|)"]), by = "term") %>%
    filter(term != "(Intercept)") %>%
    mutate(across(c(estimate), ~round(., 3)), across(c(p_OLS, p_HC3), ~round(., 4)))
}
kable(hc3_tab(m_cond_red),
      caption = "Condicional reducido: p clásico (OLS) vs robusto (HC3).",
      col.names = c("Término", "β", "p OLS", "p HC3")) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)
kable(hc3_tab(m_dir_red),
      caption = "Directo reducido: p clásico (OLS) vs robusto (HC3).",
      col.names = c("Término", "β", "p OLS", "p HC3")) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)

p_hc3_cond_red <- setNames(coeftest(m_cond_red, vcov = vcovHC(m_cond_red, "HC3"))[, "Pr(>|t|)"],
                           rownames(coeftest(m_cond_red, vcov = vcovHC(m_cond_red, "HC3"))))
p_hc3_dir_red  <- setNames(coeftest(m_dir_red,  vcov = vcovHC(m_dir_red,  "HC3"))[, "Pr(>|t|)"],
                           rownames(coeftest(m_dir_red,  vcov = vcovHC(m_dir_red,  "HC3"))))

cooks_red <- cooks.distance(m_cond_red)
u_red <- 4 / nobs(m_cond_red)
infl_red <- data.frame(obs = seq_along(cooks_red), cooks_d = cooks_red) %>%
  mutate(influyente = cooks_d > u_red)

ggplot(infl_red, aes(obs, cooks_d, color = influyente)) +
  geom_point(size = 2) +
  geom_hline(yintercept = u_red, linetype = "dashed", color = "#B2182B") +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#B2182B"),
                     labels = c("Normal", "Influyente")) +
  labs(title = "Distancia de Cook (condicional reducido)",
       subtitle = paste0("Umbral: 4/n = ", round(u_red, 3)),
       x = "Observación", y = "Cook's D", color = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

# Sensibilidad: reajuste sin el outlier real (|std_resid| > 2)
obs_real_red <- which(infl_red$influyente & abs(rstandard(m_cond_red)) > 2)
if (length(obs_real_red) > 0) {
  dcc_red <- df[complete.cases(df[, all.vars(formula(m_cond_red))]), ]
  m_sens_red <- update(m_cond_red, data = dcc_red[-obs_real_red, ])
  left_join(tidy(m_cond_red) %>% select(term, est_orig = estimate),
            tidy(m_sens_red)  %>% select(term, est_sens = estimate), by = "term") %>%
    mutate(`Δ%` = round((est_sens - est_orig) / abs(est_orig) * 100, 1)) %>%
    kable(digits = 3, caption = "Condicional reducido: con y sin el outlier real.") %>%
    kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)
}

forest_red <- tidy(m_cond_red, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(p.value < 0.05 ~ "p < 0.05",
                         p.value < 0.10 ~ "p < 0.10",
                         TRUE           ~ "p ≥ 0.10"),
         label = etiquetas(term))

p_forest_red <- ggplot(forest_red, aes(x = estimate, y = reorder(label, estimate), color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  scale_color_manual(values = c("p < 0.05" = "#B2182B",
                                "p < 0.10" = "#F4A582",
                                "p ≥ 0.10" = "#999999"),
                     name = "Significación") +
  labs(x = expression(beta ~ "(IC 95%)"), y = NULL,
       title = "Modelo condicional reducido (núcleo clínico)",
       subtitle = paste0("n = ", nobs(m_cond_red),
                         " | R² aj. = ", round(gc_red$adj.r.squared, 3))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_line(color = "grey92"),
        panel.grid.minor = element_blank(), legend.position = "bottom")
p_forest_red

dir.create("../R_output/models", recursive = TRUE, showWarnings = FALSE)

# Modelos 
saveRDS(m_cond,  "../R_output/models/m_condicional_16b.rds")
saveRDS(m_dir,   "../R_output/models/m_directo_16b.rds")

# Tablas
write.csv(tidy(m_cond, conf.int = TRUE),
          "../R_output/tables/coef_modelo_condicional_16b.csv", row.names = FALSE)
write.csv(tidy(m_dir, conf.int = TRUE),
          "../R_output/tables/coef_modelo_directo_16b.csv", row.names = FALSE)
write.csv(hc3_comp_cond,
          "../R_output/tables/hc3_modelo_condicional_16b.csv", row.names = FALSE)
write.csv(hc3_comp_dir,
          "../R_output/tables/hc3_modelo_directo_16b.csv", row.names = FALSE)

# Figura
ggsave("../R_output/figures/fig_forest_condicional_16b.pdf",
       p_forest_cond, width = 8, height = 5, dpi = 300)

# Modelo reducido 
saveRDS(m_cond_red, "../R_output/models/m_condicional_16b_reducido.rds")
saveRDS(m_dir_red,  "../R_output/models/m_directo_16b_reducido.rds")
write.csv(tidy(m_cond_red, conf.int = TRUE),
          "../R_output/tables/coef_modelo_condicional_16b_reducido.csv", row.names = FALSE)
write.csv(tidy(m_dir_red, conf.int = TRUE),
          "../R_output/tables/coef_modelo_directo_16b_reducido.csv", row.names = FALSE)
ggsave("../R_output/figures/fig_forest_condicional_16b_reducido.pdf",
       p_forest_red, width = 8, height = 5, dpi = 300)

cat("Exportación 16b completada.\n")
