#!/usr/bin/env Rscript
# Análisis funcional con MaAsLin2 sobre las predicciones de PICRUSt2.

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
MODELS_DIR  <- file.path(PROJECT_DIR, "R_output", "models")
PICRUST_DIR <- file.path(PROJECT_DIR, "qiime2_output", "picrust2_output")
WIDE_FILE   <- file.path(PROJECT_DIR, "data",
                         "nemo_wide_clean_con_leyenda_multiomics.xlsx")
OUT_BASE    <- file.path(PROJECT_DIR, "R_output", "picrust2_maaslin2")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

# (1) Cargar metadata desde phyloseq 
ps <- readRDS(file.path(PROJECT_DIR, "R_output", "models", "ps_clean.rds"))
meta_all <- as(sample_data(ps), "data.frame")

residuals_df <- readRDS(file.path(MODELS_DIR, "residuos_clinicos_fase1.rds")) %>%
  mutate(Pareja = sub("^NEMO_", "N", dyad_id)) %>%
  select(Pareja, residuos)

# Enfoque DIRECTO (sin residuos): adiposidad clínica observada, unida por Pareja.
clin <- read_excel(WIDE_FILE, sheet = "BASE DE DATOS") %>%
  transmute(
    Pareja = sub("^NEMO_", "N", dyad_id),
    WLZ_6  = as.numeric(WLZ_6),
    WLZ_0  = as.numeric(WLZ_0),
    overweight_risk_6m = factor(ifelse(as.numeric(WLZ_6) > 1, "riesgo", "no_riesgo"),
                                levels = c("no_riesgo", "riesgo"))
  )

# (2) QC de NSTI 
nsti_file <- file.path(PICRUST_DIR, "marker_predicted_and_nsti.tsv.gz")
if (file.exists(nsti_file)) {
  nsti <- read_tsv(nsti_file, show_col_types = FALSE)
  message(sprintf("\nNSTI: media=%.3f, mediana=%.3f, max=%.3f",
                  mean(nsti$metadata_NSTI, na.rm = TRUE),
                  median(nsti$metadata_NSTI, na.rm = TRUE),
                  max(nsti$metadata_NSTI, na.rm = TRUE)))
  message(sprintf("ASVs con NSTI > 2 (candidatas a excluir): %d / %d",
                  sum(nsti$metadata_NSTI > 2, na.rm = TRUE),
                  nrow(nsti)))
  write_csv(nsti, file.path(OUT_BASE, "nsti_summary.csv"))
}

# (3) Función para cargar una tabla de PICRUSt2 
load_picrust_table <- function(path) {
  tab <- read_tsv(path, show_col_types = FALSE)
  feat_ids <- tab[[1]]
  mat <- as.data.frame(tab[, -1])
  rownames(mat) <- feat_ids
  # Filas = features; columnas = muestras. Para MaAsLin2 invertimos.
  as.data.frame(t(mat))
}

ko_path <- file.path(PICRUST_DIR, "KO_metagenome_out",
                     "pred_metagenome_unstrat.tsv.gz")
pw_path <- file.path(PICRUST_DIR, "pathways_out",
                     "path_abun_unstrat.tsv.gz")
stopifnot(file.exists(ko_path), file.exists(pw_path))

ko_tab <- load_picrust_table(ko_path)
pw_tab <- load_picrust_table(pw_path)
message(sprintf("\nKO: %d muestras × %d funciones",  nrow(ko_tab), ncol(ko_tab)))
message(sprintf("MetaCyc pathways: %d muestras × %d rutas",
                nrow(pw_tab), ncol(pw_tab)))

# (4) Wrapper MaAsLin2 
run_maaslin <- function(feat_mat, meta_df, fixed_effects, label,
                        out_subdir, level, random_effects = NULL) {
  out_path <- file.path(OUT_BASE, level, out_subdir, label)
  if (dir.exists(out_path)) {
    message(sprintf("  [%s/%s/%s] ya existe — omitido",
                    level, out_subdir, label))
    return(invisible(NULL))
  }
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

  vars_needed <- c(fixed_effects, random_effects)
  keep <- complete.cases(meta_df[, vars_needed, drop = FALSE])
  meta_df  <- meta_df[keep, , drop = FALSE]
  feat_mat <- feat_mat[rownames(meta_df), , drop = FALSE]

  message(sprintf("  [%s/%s/%s] n=%d, features=%d",
                  level, out_subdir, label,
                  nrow(meta_df), ncol(feat_mat)))

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
    plot_scatter     = FALSE,
    standardize      = TRUE
  )
}

# (5) Preparar metadata por estrato 
meta_all$sample_id <- rownames(meta_all)
meta_aug <- meta_all %>%
  left_join(residuals_df, by = "Pareja") %>%
  left_join(clin, by = "Pareja") %>%
  column_to_rownames("sample_id")

split_meta <- function(strato) {
  idx <- rownames(meta_aug)[meta_aug$Origen_tiempo == strato]
  meta_aug[idx, , drop = FALSE]
}
meta_b2 <- split_meta("Bebe_Basal")
meta_b4 <- split_meta("Bebe_mes_6")
meta_m2 <- split_meta("Madre_Basal")

subset_feat <- function(feat_mat, meta_df) {
  ids <- intersect(rownames(meta_df), rownames(feat_mat))
  feat_mat[ids, , drop = FALSE]
}

# (6) Ejecutar MaAsLin2 sobre pathways y KO 
run_block <- function(feat_mat, level) {
  message(sprintf("\n=== %s — Bebé Basal (meconio) ===", level))
  ft <- subset_feat(feat_mat, meta_b2)
  run_maaslin(ft, meta_b2[rownames(ft), ], c("residuos","Peso_Pre"),       "residuos_y_pesoPre", "babies_basal", level)
  run_maaslin(ft, meta_b2[rownames(ft), ], c("WLZ_6","WLZ_0","Sexo_bebe"),     "WLZ6_adj",        "babies_basal", level)
  run_maaslin(ft, meta_b2[rownames(ft), ], c("overweight_risk_6m","Sexo_bebe"), "overweight_risk", "babies_basal", level)
  run_maaslin(ft, meta_b2[rownames(ft), ], "Peso_Pre",            "pesoPre_solo",        "babies_basal", level)
  run_maaslin(ft, meta_b2[rownames(ft), ], "Sexo_bebe",           "sexo_solo",           "babies_basal", level)

  message(sprintf("\n=== %s — Bebé 6 meses ===", level))
  ft <- subset_feat(feat_mat, meta_b4)
  run_maaslin(ft, meta_b4[rownames(ft), ], c("residuos","Peso_Pre"),       "residuos_y_pesoPre", "babies_6m", level)
  run_maaslin(ft, meta_b4[rownames(ft), ], c("WLZ_6","WLZ_0","Sexo_bebe"),     "WLZ6_adj",        "babies_6m", level)
  run_maaslin(ft, meta_b4[rownames(ft), ], c("overweight_risk_6m","Sexo_bebe"), "overweight_risk", "babies_6m", level)
  run_maaslin(ft, meta_b4[rownames(ft), ], "Peso_Pre",            "pesoPre_solo",        "babies_6m", level)
  run_maaslin(ft, meta_b4[rownames(ft), ], "Sexo_bebe",           "sexo_solo",           "babies_6m", level)

  message(sprintf("\n=== %s — Madre Basal ===", level))
  ft <- subset_feat(feat_mat, meta_m2)
  run_maaslin(ft, meta_m2[rownames(ft), ], c("residuos","Peso_Pre"),       "residuos_y_pesoPre", "mothers_basal", level)
  run_maaslin(ft, meta_m2[rownames(ft), ], c("WLZ_6","WLZ_0","Sexo_bebe"),     "WLZ6_adj",        "mothers_basal", level)
  run_maaslin(ft, meta_m2[rownames(ft), ], c("overweight_risk_6m","Sexo_bebe"), "overweight_risk", "mothers_basal", level)
  run_maaslin(ft, meta_m2[rownames(ft), ], "Peso_Pre",            "pesoPre_solo",        "mothers_basal", level)
}

run_block(pw_tab, "pathways_metacyc")
run_block(ko_tab, "ko")

# (7) Resumen consolidado de hits 
collect_hits <- function(level) {
  base <- file.path(OUT_BASE, level)
  files <- list.files(base, pattern = "significant_results\\.tsv$",
                      recursive = TRUE, full.names = TRUE)
  purrr::map_df(files, function(f) {
    rel <- sub(paste0(base, "/"), "", f)
    parts <- strsplit(rel, "/")[[1]]
    df <- tryCatch(read_tsv(f, show_col_types = FALSE),
                   error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$grupo  <- parts[1]
    df$modelo <- parts[2]
    df$nivel  <- level
    df
  })
}

hits <- bind_rows(
  collect_hits("pathways_metacyc"),
  collect_hits("ko")
) %>% arrange(nivel, grupo, modelo, qval)

write_csv(hits, file.path(OUT_BASE, "hits_consolidados.csv"))

message(sprintf("\nTotal hits significativos (q<0.25): %d", nrow(hits)))
if (nrow(hits) > 0) print(hits, n = Inf)

cat("\n11_funcional.R completado.\n")
