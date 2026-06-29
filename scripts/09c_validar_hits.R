#!/usr/bin/env Rscript
# Validación de los hits de MaAsLin2 (prev10):

suppressPackageStartupMessages({
  library(Maaslin2)
  library(phyloseq)
  library(tidyverse)
})

# Rutas
.file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
PROJECT_DIR <- if (length(.file_arg))
  normalizePath(file.path(dirname(sub("^--file=", "", .file_arg)), "..")) else normalizePath("..")
OUT_BASE    <- file.path(PROJECT_DIR, "R_output", "maaslin2_prev_sanity")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

ps <- readRDS(file.path(PROJECT_DIR, "R_output", "models", "ps_clean.rds"))
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)

tt <- as.data.frame(tax_table(ps_genus), stringsAsFactors = FALSE)
genus_labels <- ifelse(is.na(tt$Genus) | tt$Genus == "",
                       paste0("UnknownG_", seq_len(nrow(tt))),
                       tt$Genus)
genus_labels <- make.unique(genus_labels, sep = "_")
taxa_names(ps_genus) <- genus_labels
tt$label <- genus_labels

# (a) Validación taxonómica 
hits <- c("UnknownG_442", "Gemella", "GCA-900066575", "GCA.900066575")
cat("\n=== Taxonomía de los hits ===\n")
hit_rows <- tt[tt$label %in% hits |
                 tt$Genus %in% gsub("\\.", "-", hits), ]
print(hit_rows[, intersect(c("Kingdom","Phylum","Class","Order","Family","Genus","label"),
                            colnames(hit_rows))])

# Guardamos la tabla completa para referencia (sólo géneros desconocidos)
unk <- tt[grepl("^UnknownG_", tt$label), ]
write.csv(unk, file.path(OUT_BASE, "unknown_genera_taxonomy.csv"), row.names = FALSE)
cat(sprintf("\nGéneros 'UnknownG_*' guardados en %s (n=%d)\n",
            file.path(OUT_BASE, "unknown_genera_taxonomy.csv"), nrow(unk)))

# (b) Sanity check con prev 15% y 20% 
run_maaslin_prev <- function(ps_sub, fixed_effects, label, out_subdir, prev) {
  out_path <- file.path(OUT_BASE, sprintf("prev%02d", round(prev*100)),
                        out_subdir, label)
  if (dir.exists(out_path)) {
    message(sprintf("  [prev%02d/%s/%s] ya existe — omitido",
                    round(prev*100), out_subdir, label))
    return(invisible(NULL))
  }
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

  feat_mat <- as(otu_table(ps_sub), "matrix")
  if (taxa_are_rows(ps_sub)) feat_mat <- t(feat_mat)
  feat_mat <- as.data.frame(feat_mat)

  meta_df <- as(sample_data(ps_sub), "data.frame")
  keep <- complete.cases(meta_df[, fixed_effects, drop = FALSE])
  meta_df  <- meta_df[keep, ]
  feat_mat <- feat_mat[rownames(meta_df), ]

  message(sprintf("  [prev%02d/%s/%s] n=%d, géneros=%d",
                  round(prev*100), out_subdir, label,
                  nrow(meta_df), ncol(feat_mat)))

  Maaslin2(
    input_data       = feat_mat,
    input_metadata   = meta_df,
    output           = out_path,
    fixed_effects    = fixed_effects,
    normalization    = "TSS",
    transform        = "LOG",
    analysis_method  = "LM",
    min_prevalence   = prev,
    min_abundance    = 0.0,
    max_significance = 0.25,
    correction       = "BH",
    plot_heatmap     = TRUE,
    plot_scatter     = TRUE,
    standardize      = TRUE
  )
}

ps_b4 <- subset_samples(ps_genus, Origen_tiempo == "Bebe_mes_6")
ps_m2 <- subset_samples(ps_genus,
                        Origen_tiempo == "Madre_Basal" |
                          grepl("M1$", sample_names(ps_genus)))

for (prev in c(0.15, 0.20)) {
  message(sprintf("\n=== Bebé 6 m — prev %d%% ===", round(prev*100)))
  run_maaslin_prev(ps_b4, "Peso_Pre", "pesoPre_solo", "babies_6m", prev)
  message(sprintf("\n=== Madre Basal — prev %d%% ===", round(prev*100)))
  run_maaslin_prev(ps_m2, "Peso_Pre", "pesoPre_solo", "mothers_basal", prev)
}

cat("\n09c_validar_hits.R completado.\n")
