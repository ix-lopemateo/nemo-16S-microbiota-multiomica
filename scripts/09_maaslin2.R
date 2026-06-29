#!/usr/bin/env Rscript
# Abundancia diferencial con MaAsLin2 

suppressPackageStartupMessages({
  library(Maaslin2)
  library(phyloseq)
  library(readxl)
  library(tidyverse)
})

# Rutas 
.file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
PROJECT_DIR <- if (length(.file_arg))
  normalizePath(file.path(dirname(sub("^--file=", "", .file_arg)), "..")) else normalizePath("..")
MODELS_DIR    <- file.path(PROJECT_DIR, "R_output", "models")
TABLES_DIR    <- file.path(PROJECT_DIR, "R_output", "models")
OUT_BASE      <- file.path(PROJECT_DIR, "R_output", "maaslin2")
WIDE_FILE     <- file.path(PROJECT_DIR, "data",
                            "nemo_wide_clean_con_leyenda_multiomics.xlsx")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

# Cargar phyloseq limpio (de 08_composicion.R) 
ps <- readRDS(file.path(PROJECT_DIR, "R_output", "models", "ps_clean.rds"))
message(sprintf("Phyloseq cargado: %d muestras, %d ASVs",
                nsamples(ps), ntaxa(ps)))

# Agregar a nivel de género
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
message(sprintf("Tras agregar a género: %d géneros", ntaxa(ps_genus)))

# Renombrar features a "Genus_<idx>" con etiqueta legible para MaAsLin2
tt <- as.data.frame(tax_table(ps_genus), stringsAsFactors = FALSE)
genus_labels <- ifelse(is.na(tt$Genus) | tt$Genus == "",
                       paste0("UnknownG_", seq_len(nrow(tt))),
                       tt$Genus)
# Hacer únicos
genus_labels <- make.unique(genus_labels, sep = "_")
taxa_names(ps_genus) <- genus_labels

# Cargar residuos de Fase 1 (n=30) 
residuals_df <- readRDS(file.path(MODELS_DIR, "residuos_clinicos_fase1.rds"))
message(sprintf("Residuos de Fase 1 cargados: %d díadas", nrow(residuals_df)))

# Mapa NEMO_XX → NXX para unir con IDs de microbioma
residuals_df <- residuals_df %>%
  mutate(Pareja = sub("^NEMO_", "N", dyad_id))

# Enfoque DIRECTO (sin residuos): adiposidad clínica unida por Pareja desde el wide.
# overweight_risk_6m = WLZ_6 > +1 (referencia "no_riesgo").
clin <- read_excel(WIDE_FILE, sheet = "BASE DE DATOS") %>%
  transmute(
    Pareja = sub("^NEMO_", "N", dyad_id),
    WLZ_6  = as.numeric(WLZ_6),
    WLZ_0  = as.numeric(WLZ_0),
    overweight_risk_6m = factor(ifelse(as.numeric(WLZ_6) > 1, "riesgo", "no_riesgo"),
                                levels = c("no_riesgo", "riesgo"))
  )
message(sprintf("Adiposidad directa cargada: %d díadas con WLZ_6", sum(!is.na(clin$WLZ_6))))

# Función genérica para correr MaAsLin2 sobre un subconjunto 
run_maaslin <- function(ps_sub, meta_extra, fixed_effects, label,
                        out_subdir, random_effects = NULL) {
  out_path <- file.path(OUT_BASE, out_subdir, label)
  if (dir.exists(out_path)) {
    message(sprintf("  [%s/%s] ya existe — omitido", out_subdir, label))
    return(invisible(NULL))
  }
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

  # Matriz de features: samples × genus
  feat_mat <- as(otu_table(ps_sub), "matrix")
  if (taxa_are_rows(ps_sub)) feat_mat <- t(feat_mat)
  feat_mat <- as.data.frame(feat_mat)

  # Metadata: combinar sample_data con variables extra
  meta_df <- as(sample_data(ps_sub), "data.frame")
  meta_df$sample_id <- rownames(meta_df)
  if (!is.null(meta_extra)) {
    meta_df <- left_join(meta_df, meta_extra, by = "Pareja")
  }
  rownames(meta_df) <- meta_df$sample_id
  meta_df$sample_id <- NULL

  # Filtrar a muestras con datos completos en todas las variables fijas
  vars_needed <- c(fixed_effects, random_effects)
  keep <- complete.cases(meta_df[, vars_needed, drop = FALSE])
  meta_df  <- meta_df[keep, ]
  feat_mat <- feat_mat[rownames(meta_df), ]

  if (nrow(meta_df) < 10) {
    message(sprintf("  [%s/%s] saltado (n=%d insuficiente)", out_subdir, label, nrow(meta_df)))
    return(invisible(NULL))
  }

  message(sprintf("  [%s/%s] n=%d, géneros=%d",
                  out_subdir, label, nrow(meta_df), ncol(feat_mat)))

  Maaslin2(
    input_data      = feat_mat,
    input_metadata  = meta_df,
    output          = out_path,
    fixed_effects   = fixed_effects,
    random_effects  = random_effects,
    normalization   = "TSS",
    transform       = "LOG",
    analysis_method = "LM",
    min_prevalence  = 0.20,
    min_abundance   = 0.0,
    max_significance = 0.25,
    correction      = "BH",
    plot_heatmap    = TRUE,
    plot_scatter    = TRUE,
    standardize     = TRUE
  )
}

# Subconjuntos por grupo 
ps_b2 <- subset_samples(ps_genus, Origen_tiempo == "Bebe_Basal")
ps_b4 <- subset_samples(ps_genus, Origen_tiempo == "Bebe_mes_6")
ps_m2 <- subset_samples(ps_genus,
                        Origen_tiempo == "Madre_Basal" |
                          grepl("M1$", sample_names(ps_genus)))

# Variables clínicas a unir desde residuals_df (sólo "residuos")
meta_res <- residuals_df %>% select(Pareja, residuos)

# BEBÉ BASAL (B2)
message("\n=== Bebé Basal (meconio) ===")
# (A) Enfoque CON residuos de Fase 1
run_maaslin(ps_b2, meta_res, fixed_effects = c("residuos", "Peso_Pre"),
            label = "residuos_y_pesoPre", out_subdir = "babies_basal")
# (B) Enfoque DIRECTO (adiposidad observada, ajustada)
run_maaslin(ps_b2, clin, fixed_effects = c("WLZ_6", "WLZ_0", "Sexo_bebe"),
            label = "WLZ6_adj", out_subdir = "babies_basal")
run_maaslin(ps_b2, clin, fixed_effects = c("overweight_risk_6m", "Sexo_bebe"),
            label = "overweight_risk", out_subdir = "babies_basal")
# (C) Modelos simples (referencia)
run_maaslin(ps_b2, NULL, fixed_effects = "Peso_Pre",
            label = "pesoPre_solo", out_subdir = "babies_basal")
run_maaslin(ps_b2, NULL, fixed_effects = "Sexo_bebe",
            label = "sexo_solo", out_subdir = "babies_basal")

# BEBÉ 6 MESES (B4)
message("\n=== Bebé 6 meses ===")
# (A) Enfoque CON residuos de Fase 1
run_maaslin(ps_b4, meta_res, fixed_effects = c("residuos", "Peso_Pre"),
            label = "residuos_y_pesoPre", out_subdir = "babies_6m")
# (B) Enfoque DIRECTO (adiposidad observada, ajustada)
run_maaslin(ps_b4, clin, fixed_effects = c("WLZ_6", "WLZ_0", "Sexo_bebe"),
            label = "WLZ6_adj", out_subdir = "babies_6m")
run_maaslin(ps_b4, clin, fixed_effects = c("overweight_risk_6m", "Sexo_bebe"),
            label = "overweight_risk", out_subdir = "babies_6m")
# (C) Modelos simples (referencia)
run_maaslin(ps_b4, NULL, fixed_effects = "Peso_Pre",
            label = "pesoPre_solo", out_subdir = "babies_6m")
run_maaslin(ps_b4, NULL, fixed_effects = "Sexo_bebe",
            label = "sexo_solo", out_subdir = "babies_6m")

# MADRE BASAL (M2/M1)
message("\n=== Madre Basal ===")
# (A) Enfoque CON residuos de Fase 1
run_maaslin(ps_m2, meta_res, fixed_effects = c("residuos", "Peso_Pre"),
            label = "residuos_y_pesoPre", out_subdir = "mothers_basal")
# (B) Enfoque DIRECTO (microbiota materna vs adiposidad del lactante)
run_maaslin(ps_m2, clin, fixed_effects = c("WLZ_6", "WLZ_0", "Sexo_bebe"),
            label = "WLZ6_adj", out_subdir = "mothers_basal")
run_maaslin(ps_m2, clin, fixed_effects = c("overweight_risk_6m", "Sexo_bebe"),
            label = "overweight_risk", out_subdir = "mothers_basal")
# (C) Modelos simples (referencia)
run_maaslin(ps_m2, NULL, fixed_effects = "Peso_Pre",
            label = "pesoPre_solo", out_subdir = "mothers_basal")

message("\n09_maaslin2.R completado.")
message("Resultados en: ", OUT_BASE)
