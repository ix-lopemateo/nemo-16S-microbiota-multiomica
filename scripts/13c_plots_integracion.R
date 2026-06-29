#!/usr/bin/env Rscript
# Visualización de los 4 patrones sub-umbral más coherentes de la integración

suppressPackageStartupMessages({
  library(readxl)
  library(phyloseq)
  library(tidyverse)
})

# Rutas
.file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
PROJECT_DIR <- if (length(.file_arg))
  normalizePath(file.path(dirname(sub("^--file=", "", .file_arg)), "..")) else normalizePath("..")
MODELS_DIR  <- file.path(PROJECT_DIR, "R_output", "models")
MET_FILE    <- file.path(PROJECT_DIR, "data",
                         "pirosecuenciacion_metilacion_crudo_limpio.xlsx")
WIDE_FILE   <- file.path(PROJECT_DIR, "data",
                         "nemo_wide_clean_con_leyenda_multiomics.xlsx")
PIC_DIR     <- file.path(PROJECT_DIR, "qiime2_output", "picrust2_output")
OUT_DIR     <- file.path(PROJECT_DIR, "R_output", "integracion_metilacion",
                         "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Tema base 
theme_paper <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = NA),
        plot.title.position = "plot")

# Helper: scatter + Spearman 
scatter_spearman <- function(df, x, y, title = NULL, subtitle = NULL,
                             xlab = x, ylab = y, color = "#2E86AB") {
  ct <- suppressWarnings(cor.test(df[[x]], df[[y]],
                                  method = "spearman", exact = FALSE))
  ann <- sprintf("rho == %.2f * ',' ~ italic(p) == %.3f * ',' ~ italic(n) == %d",
                 unname(ct$estimate), ct$p.value, sum(complete.cases(df[, c(x, y)])))
  ggplot(df, aes(.data[[x]], .data[[y]])) +
    geom_smooth(method = "lm", se = TRUE, colour = color, fill = color,
                alpha = 0.15, linewidth = 0.6) +
    geom_point(size = 2.2, alpha = 0.85, colour = color) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
             label = ann, parse = TRUE, size = 3.2) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme_paper
}

# Carga de datos
meth <- read_excel(MET_FILE, sheet = "wide")
# La hoja "wide" nombra la primera columna como `x`; la renombramos a `Pareja`.
if (!"Pareja" %in% names(meth) && "x" %in% names(meth)) {
  meth <- dplyr::rename(meth, Pareja = x)
}

ps <- readRDS(file.path(PROJECT_DIR, "R_output", "models", "ps_clean.rds"))
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
tt <- as.data.frame(tax_table(ps_genus), stringsAsFactors = FALSE)
genus_labels <- ifelse(is.na(tt$Genus) | tt$Genus == "",
                       paste0("UnknownG_", seq_len(nrow(tt))),
                       tt$Genus)
taxa_names(ps_genus) <- make.unique(genus_labels, sep = "_")

counts <- as(otu_table(ps_genus), "matrix")
if (taxa_are_rows(ps_genus)) counts <- t(counts)
rel <- as.data.frame(sweep(counts, 1, rowSums(counts), "/"))
meta_micro <- as(sample_data(ps_genus), "data.frame")
meta_micro$sample_id <- rownames(meta_micro)

residuals_df <- readRDS(file.path(MODELS_DIR, "residuos_clinicos_fase1.rds")) %>%
  mutate(Pareja = sub("^NEMO_", "N", dyad_id)) %>%
  select(Pareja, residuos)

# Enfoque DIRECTO (sin residuos): adiposidad observada (WLZ_6) por Pareja
clin <- read_excel(WIDE_FILE, sheet = "BASE DE DATOS") %>%
  transmute(Pareja = sub("^NEMO_", "N", dyad_id),
            WLZ_6  = as.numeric(WLZ_6)) %>%
  filter(!is.na(WLZ_6))

# FIG 1: FTO_mean (4 tiempos) × residuos clínicos
build_fto_panel <- function(tiempo, label_tiempo) {
  meth %>%
    filter(Tiempo == tiempo) %>%
    select(Pareja, FTO_mean) %>%
    inner_join(residuals_df, by = "Pareja") %>%
    drop_na() %>%
    mutate(panel = label_tiempo)
}
fto_df <- bind_rows(
  build_fto_panel("m2", "Madre basal (m2)"),
  build_fto_panel("b2", "Bebé basal/meconio (b2)"),
  build_fto_panel("b3", "Bebé 1 mes (b3)"),
  build_fto_panel("b4", "Bebé 6 meses (b4)")
) %>%
  mutate(panel = factor(panel, levels = c("Madre basal (m2)",
                                          "Bebé basal/meconio (b2)",
                                          "Bebé 1 mes (b3)",
                                          "Bebé 6 meses (b4)")))

# Calcular rho/p por panel para anotar
stats_per_panel <- fto_df %>%
  group_by(panel) %>%
  summarise(rho = cor(residuos, FTO_mean, method = "spearman"),
            pval = suppressWarnings(cor.test(residuos, FTO_mean,
                                             method = "spearman",
                                             exact = FALSE)$p.value),
            n = n(),
            .groups = "drop") %>%
  mutate(label = sprintf("rho==%.2f*','~italic(p)==%.3f*','~italic(n)==%d",
                         rho, pval, n))

fig1 <- ggplot(fto_df, aes(residuos, FTO_mean)) +
  geom_smooth(method = "lm", se = TRUE, colour = "#A23B72",
              fill = "#A23B72", alpha = 0.15, linewidth = 0.6) +
  geom_point(size = 2.2, alpha = 0.85, colour = "#A23B72") +
  geom_text(data = stats_per_panel,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.4, size = 3, parse = TRUE,
            inherit.aes = FALSE) +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  labs(title = "Metilación de FTO en saliva vs.\\ residuos clínicos de adiposidad",
       subtitle = "Patrón consistente a lo largo de 4 puntos temporales",
       x = "Residuos clínicos de WLZ a 6 m (Fase 1)",
       y = "Metilación promedio de FTO (%)") +
  theme_paper
ggsave(file.path(OUT_DIR, "fig1_FTO_vs_residuos.pdf"),
       fig1, width = 8.5, height = 6.5)

# FIG 1b: FTO_mean (4 tiempos) × WLZ_6 observado (enfoque SIN residuos)
build_fto_panel_wlz <- function(tiempo, label_tiempo) {
  meth %>%
    filter(Tiempo == tiempo) %>%
    select(Pareja, FTO_mean) %>%
    inner_join(clin, by = "Pareja") %>%
    drop_na() %>%
    mutate(panel = label_tiempo)
}
fto_wlz <- bind_rows(
  build_fto_panel_wlz("m2", "Madre basal (m2)"),
  build_fto_panel_wlz("b2", "Bebé basal/meconio (b2)"),
  build_fto_panel_wlz("b3", "Bebé 1 mes (b3)"),
  build_fto_panel_wlz("b4", "Bebé 6 meses (b4)")
) %>%
  mutate(panel = factor(panel, levels = c("Madre basal (m2)",
                                          "Bebé basal/meconio (b2)",
                                          "Bebé 1 mes (b3)",
                                          "Bebé 6 meses (b4)")))

stats_wlz <- fto_wlz %>%
  group_by(panel) %>%
  summarise(rho = cor(WLZ_6, FTO_mean, method = "spearman"),
            pval = suppressWarnings(cor.test(WLZ_6, FTO_mean,
                                             method = "spearman",
                                             exact = FALSE)$p.value),
            n = n(),
            .groups = "drop") %>%
  mutate(label = sprintf("rho==%.2f*','~italic(p)==%.3f*','~italic(n)==%d",
                         rho, pval, n))

fig1b <- ggplot(fto_wlz, aes(WLZ_6, FTO_mean)) +
  geom_smooth(method = "lm", se = TRUE, colour = "#2E86AB",
              fill = "#2E86AB", alpha = 0.15, linewidth = 0.6) +
  geom_point(size = 2.2, alpha = 0.85, colour = "#2E86AB") +
  geom_text(data = stats_wlz,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.4, size = 3, parse = TRUE,
            inherit.aes = FALSE) +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  labs(title = "Metilación de FTO en saliva vs.\\ adiposidad observada (WLZ a 6 m)",
       subtitle = "Enfoque sin residuos — comparar con el panel basado en residuos",
       x = "WLZ a 6 meses del lactante",
       y = "Metilación promedio de FTO (%)") +
  theme_paper
ggsave(file.path(OUT_DIR, "fig1b_FTO_vs_WLZ6.pdf"),
       fig1b, width = 8.5, height = 6.5)

# FIG 5: NRF1_mean (4 tiempos) × WLZ_6 observado (enfoque DIRECTO)
# El hit de la integración (B4_adiposidad) es NRF1 a 1 mes (b3) vs WLZ_6:
# rho = -0,51, q = 0,15. A diferencia de FTO, el signo NO es consistente en los
# cuatro tiempos: solo alcanza significación en b3.
build_nrf1_panel_wlz <- function(tiempo, label_tiempo) {
  meth %>%
    filter(Tiempo == tiempo) %>%
    select(Pareja, NRF1_mean) %>%
    inner_join(clin, by = "Pareja") %>%
    drop_na() %>%
    mutate(panel = label_tiempo)
}
nrf1_wlz <- bind_rows(
  build_nrf1_panel_wlz("m2", "Madre basal (m2)"),
  build_nrf1_panel_wlz("b2", "Bebé basal/meconio (b2)"),
  build_nrf1_panel_wlz("b3", "Bebé 1 mes (b3)"),
  build_nrf1_panel_wlz("b4", "Bebé 6 meses (b4)")
) %>%
  mutate(panel = factor(panel, levels = c("Madre basal (m2)",
                                          "Bebé basal/meconio (b2)",
                                          "Bebé 1 mes (b3)",
                                          "Bebé 6 meses (b4)")))

stats_nrf1 <- nrf1_wlz %>%
  group_by(panel) %>%
  summarise(rho = cor(WLZ_6, NRF1_mean, method = "spearman"),
            pval = suppressWarnings(cor.test(WLZ_6, NRF1_mean,
                                             method = "spearman",
                                             exact = FALSE)$p.value),
            n = n(),
            .groups = "drop") %>%
  mutate(label = sprintf("rho==%.2f*','~italic(p)==%.3f*','~italic(n)==%d",
                         rho, pval, n))

fig5 <- ggplot(nrf1_wlz, aes(WLZ_6, NRF1_mean)) +
  geom_smooth(method = "lm", se = TRUE, colour = "#6A51A3",
              fill = "#6A51A3", alpha = 0.15, linewidth = 0.6) +
  geom_point(size = 2.2, alpha = 0.85, colour = "#6A51A3") +
  geom_text(data = stats_nrf1,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.4, size = 3, parse = TRUE,
            inherit.aes = FALSE) +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  labs(title = "Metilación de NRF1 en saliva vs.\\ adiposidad observada (WLZ a 6 m)",
       subtitle = "Enfoque directo — solo alcanza significación a 1 mes (b3); sin la consistencia temporal de FTO",
       x = "WLZ a 6 meses del lactante",
       y = "Metilación promedio de NRF1 (%)") +
  theme_paper
ggsave(file.path(OUT_DIR, "fig5_NRF1_vs_WLZ6.pdf"),
       fig5, width = 8.5, height = 6.5)
ggsave(file.path(OUT_DIR, "fig5_NRF1_vs_WLZ6.png"),
       fig5, width = 8.5, height = 6.5, dpi = 150)

# FIG 2: GCA-900066575 materna × MC4R bebé 6m (eje vertical)
gca_ids <- meta_micro$sample_id[meta_micro$Origen_tiempo == "Madre_Basal"]
gca_df <- tibble(Pareja = meta_micro[gca_ids, "Pareja"],
                 micro  = rel[gca_ids, "GCA-900066575"]) %>%
  inner_join(meth %>% filter(Tiempo == "b4") %>%
               select(Pareja, MC4R_mean), by = "Pareja") %>%
  drop_na()

fig2 <- scatter_spearman(
  gca_df, "micro", "MC4R_mean",
  title    = "Eje vertical madre → hijo (a 6 meses)",
  subtitle = "Abundancia de GCA-900066575 en la microbiota materna basal vs.\\ metilación de MC4R en saliva del bebé a 6 m",
  xlab     = expression("Abundancia relativa "*italic("GCA-900066575")*" (madre basal)"),
  ylab     = "Metilación promedio de MC4R en bebé a 6 m (%)",
  color    = "#2E86AB"
)
ggsave(file.path(OUT_DIR, "fig2_GCA_madre_MC4R_bebe6m.pdf"),
       fig2, width = 6.5, height = 5)

# FIG 3: Erysipelatoclostridium × MC4R en bebés 6m (intra-individuo) 
ery_ids <- meta_micro$sample_id[meta_micro$Origen_tiempo == "Bebe_mes_6"]
ery_df <- tibble(Pareja = meta_micro[ery_ids, "Pareja"],
                 micro  = rel[ery_ids, "Erysipelatoclostridium"]) %>%
  inner_join(meth %>% filter(Tiempo == "b4") %>%
               select(Pareja, MC4R_mean), by = "Pareja") %>%
  drop_na()

fig3 <- scatter_spearman(
  ery_df, "micro", "MC4R_mean",
  title    = "Asociación intra-individuo en bebés a 6 meses",
  subtitle = "Abundancia de Erysipelatoclostridium vs.\\ metilación de MC4R (misma muestra)",
  xlab     = expression("Abundancia relativa "*italic("Erysipelatoclostridium")),
  ylab     = "Metilación promedio de MC4R (%)",
  color    = "#F18F01"
)
ggsave(file.path(OUT_DIR, "fig3_Erysipelato_MC4R_bebe6m.pdf"),
       fig3, width = 6.5, height = 5)

# ── FIG 4: Rutas lipogénicas maternas × BCAT1 materno (coherencia funcional) ─
pw_file <- file.path(PIC_DIR, "pathways_out", "path_abun_unstrat.tsv.gz")
if (file.exists(pw_file)) {
  pw_raw <- read_tsv(pw_file, show_col_types = FALSE)
  pw_feat_names <- pw_raw[[1]]
  pw_mat <- as.data.frame(t(as.matrix(pw_raw[, -1])))
  colnames(pw_mat) <- pw_feat_names
  pw_mat <- pw_mat / rowSums(pw_mat)

  fa_pathways <- c("FASYN-INITIAL-PWY", "PWY-5971", "PWY-5989",
                   "PWY-6282", "PWY-7664", "PWY0-862", "PWYG-321")
  fa_labels   <- c("FA biosynth. initiation", "Palmitate biosynth.",
                   "Stearate biosynth.", "Palmitoleate biosynth.",
                   "Oleate biosynth.", "Dodec-5-enoate biosynth.",
                   "Mycolate biosynth.")
  names(fa_labels) <- fa_pathways

  mad_ids <- meta_micro$sample_id[meta_micro$Origen_tiempo == "Madre_Basal"]
  mad_ids <- intersect(mad_ids, rownames(pw_mat))
  fa_df <- pw_mat[mad_ids, fa_pathways, drop = FALSE] %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(-sample_id, names_to = "pathway", values_to = "abundancia") %>%
    mutate(Pareja = meta_micro[sample_id, "Pareja"]) %>%
    inner_join(meth %>% filter(Tiempo == "m2") %>%
                 select(Pareja, BCAT1_mean), by = "Pareja") %>%
    drop_na() %>%
    mutate(pathway_lbl = fa_labels[pathway],
           pathway_lbl = factor(pathway_lbl, levels = fa_labels))

  stats_fa <- fa_df %>%
    group_by(pathway_lbl) %>%
    summarise(rho = cor(abundancia, BCAT1_mean, method = "spearman"),
              pval = suppressWarnings(cor.test(abundancia, BCAT1_mean,
                                               method = "spearman",
                                               exact = FALSE)$p.value),
              n = n(),
              .groups = "drop") %>%
    mutate(label = sprintf("rho==%.2f*','~italic(p)==%.3f", rho, pval))

  fig4 <- ggplot(fa_df, aes(abundancia, BCAT1_mean)) +
    geom_smooth(method = "lm", se = TRUE, colour = "#3B8132",
                fill = "#3B8132", alpha = 0.15, linewidth = 0.6) +
    geom_point(size = 1.9, alpha = 0.8, colour = "#3B8132") +
    geom_text(data = stats_fa,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.05, vjust = 1.4, size = 2.8, parse = TRUE,
              inherit.aes = FALSE) +
    facet_wrap(~ pathway_lbl, scales = "free", ncol = 3) +
    labs(title = "Rutas lipogénicas maternas vs.\\ metilación BCAT1 (madre)",
         subtitle = "Coherencia funcional intra-individuo: lipogénesis bacteriana ↔ metabolismo de aminoácidos ramificados (BCAT1)",
         x = "Abundancia relativa predicha de la ruta MetaCyc",
         y = "Metilación promedio de BCAT1 (%)") +
    theme_paper +
    theme(strip.text = element_text(size = 8))

  ggsave(file.path(OUT_DIR, "fig4_FA_pathways_BCAT1_madres.pdf"),
         fig4, width = 9, height = 6.5)
} else {
  message("PICRUSt2 pathways no encontrado; se omite fig4.")
}

cat("\nFiguras generadas:\n")
for (f in list.files(OUT_DIR, pattern = "^fig.*\\.pdf$", full.names = TRUE)) {
  cat(sprintf("  %s\n", f))
}
cat("\n13c_plots_integracion.R completado.\n")
