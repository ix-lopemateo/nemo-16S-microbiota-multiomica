#!/usr/bin/env Rscript
# Re-ejecución de MaAsLin2 con prevalencia mínima reducida a 10%.
# Mismo diseño que 09_maaslin2.R pero permitiendo géneros menos abundantes

suppressPackageStartupMessages({
  library(Maaslin2)
  library(phyloseq)
  library(readxl)
  library(tidyverse)
})

# Raíz del proyecto (relativa a la ubicación del script). Para ejecución
# interactiva sin Rscript, asignar PROJECT_DIR manualmente.
.file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
PROJECT_DIR <- if (length(.file_arg))
  normalizePath(file.path(dirname(sub("^--file=", "", .file_arg)), "..")) else normalizePath("..")
MODELS_DIR  <- file.path(PROJECT_DIR, "R_output", "models")
OUT_BASE    <- file.path(PROJECT_DIR, "R_output", "maaslin2_prev10")
WIDE_FILE   <- file.path(PROJECT_DIR, "data",
                          "nemo_wide_clean_con_leyenda_multiomics.xlsx")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

ps <- readRDS(file.path(PROJECT_DIR, "R_output", "models", "ps_clean.rds"))
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)

tt <- as.data.frame(tax_table(ps_genus), stringsAsFactors = FALSE)
genus_labels <- ifelse(is.na(tt$Genus) | tt$Genus == "",
                       paste0("UnknownG_", seq_len(nrow(tt))),
                       tt$Genus)
genus_labels <- make.unique(genus_labels, sep = "_")
taxa_names(ps_genus) <- genus_labels

residuals_df <- readRDS(file.path(MODELS_DIR, "residuos_clinicos_fase1.rds"))
residuals_df <- residuals_df %>%
  mutate(Pareja = sub("^NEMO_", "N", dyad_id))

# Enfoque DIRECTO (sin residuos): adiposidad clínica unida por Pareja desde el wide.
clin <- read_excel(WIDE_FILE, sheet = "BASE DE DATOS") %>%
  transmute(
    Pareja = sub("^NEMO_", "N", dyad_id),
    WLZ_6  = as.numeric(WLZ_6),
    WLZ_0  = as.numeric(WLZ_0),
    overweight_risk_6m = factor(ifelse(as.numeric(WLZ_6) > 1, "riesgo", "no_riesgo"),
                                levels = c("no_riesgo", "riesgo"))
  )

run_maaslin <- function(ps_sub, meta_extra, fixed_effects, label,
                        out_subdir, random_effects = NULL) {
  out_path <- file.path(OUT_BASE, out_subdir, label)
  if (dir.exists(out_path)) {
    message(sprintf("  [%s/%s] ya existe — omitido", out_subdir, label))
    return(invisible(NULL))
  }
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

  feat_mat <- as(otu_table(ps_sub), "matrix")
  if (taxa_are_rows(ps_sub)) feat_mat <- t(feat_mat)
  feat_mat <- as.data.frame(feat_mat)

  meta_df <- as(sample_data(ps_sub), "data.frame")
  meta_df$sample_id <- rownames(meta_df)
  if (!is.null(meta_extra)) {
    meta_df <- left_join(meta_df, meta_extra, by = "Pareja")
  }
  rownames(meta_df) <- meta_df$sample_id
  meta_df$sample_id <- NULL

  vars_needed <- c(fixed_effects, random_effects)
  keep <- complete.cases(meta_df[, vars_needed, drop = FALSE])
  meta_df  <- meta_df[keep, ]
  feat_mat <- feat_mat[rownames(meta_df), ]

  if (nrow(meta_df) < 10) {
    message(sprintf("  [%s/%s] saltado (n=%d)", out_subdir, label, nrow(meta_df)))
    return(invisible(NULL))
  }

  message(sprintf("  [%s/%s] n=%d, géneros=%d",
                  out_subdir, label, nrow(meta_df), ncol(feat_mat)))

  Maaslin2(
    input_data       = feat_mat,
    input_metadata   = meta_df,
    output           = out_path,
    fixed_effects    = fixed_effects,
    random_effects   = random_effects,
    normalization    = "TSS",
    transform        = "LOG",
    analysis_method  = "LM",
    min_prevalence   = 0.10,
    min_abundance    = 0.0,
    max_significance = 0.25,
    correction       = "BH",
    plot_heatmap     = TRUE,
    plot_scatter     = TRUE,
    standardize      = TRUE
  )
}

ps_b2 <- subset_samples(ps_genus, Origen_tiempo == "Bebe_Basal")
ps_b4 <- subset_samples(ps_genus, Origen_tiempo == "Bebe_mes_6")
ps_m2 <- subset_samples(ps_genus,
                        Origen_tiempo == "Madre_Basal" |
                          grepl("M1$", sample_names(ps_genus)))

meta_res <- residuals_df %>% select(Pareja, residuos)

message("\n=== Bebé Basal (meconio) — prev 10% ===")
run_maaslin(ps_b2, meta_res, c("residuos", "Peso_Pre"),
            "residuos_y_pesoPre", "babies_basal")
run_maaslin(ps_b2, clin, c("WLZ_6", "WLZ_0", "Sexo_bebe"),
            "WLZ6_adj", "babies_basal")
run_maaslin(ps_b2, clin, c("overweight_risk_6m", "Sexo_bebe"),
            "overweight_risk", "babies_basal")
run_maaslin(ps_b2, NULL, "Peso_Pre",  "pesoPre_solo", "babies_basal")
run_maaslin(ps_b2, NULL, "Sexo_bebe", "sexo_solo",    "babies_basal")

message("\n=== Bebé 6 meses — prev 10% ===")
run_maaslin(ps_b4, meta_res, c("residuos", "Peso_Pre"),
            "residuos_y_pesoPre", "babies_6m")
run_maaslin(ps_b4, clin, c("WLZ_6", "WLZ_0", "Sexo_bebe"),
            "WLZ6_adj", "babies_6m")
run_maaslin(ps_b4, clin, c("overweight_risk_6m", "Sexo_bebe"),
            "overweight_risk", "babies_6m")
run_maaslin(ps_b4, NULL, "Peso_Pre",  "pesoPre_solo", "babies_6m")
run_maaslin(ps_b4, NULL, "Sexo_bebe", "sexo_solo",    "babies_6m")

message("\n=== Madre Basal — prev 10% ===")
run_maaslin(ps_m2, meta_res, c("residuos", "Peso_Pre"),
            "residuos_y_pesoPre", "mothers_basal")
run_maaslin(ps_m2, clin, c("WLZ_6", "WLZ_0", "Sexo_bebe"),
            "WLZ6_adj", "mothers_basal")
run_maaslin(ps_m2, clin, c("overweight_risk_6m", "Sexo_bebe"),
            "overweight_risk", "mothers_basal")
run_maaslin(ps_m2, NULL, "Peso_Pre", "pesoPre_solo", "mothers_basal")

message("\n09b_maaslin2_prev10.R completado.")
