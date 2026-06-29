#!/usr/bin/env Rscript
# Importa artefactos QIIME2 a phyloseq y prepara objetos por grupo
# Requiere: 01_table.qza, taxonomy.qza, rooted-tree.qza
# Produce:  ps_babies_basal.rds, ps_babies_6m.rds, ps_mothers.rds
#
# Dependencias: qiime2R, phyloseq, microbiome, tidyverse

library(qiime2R)
library(phyloseq)
library(microbiome)
library(tidyverse)

# ── Rutas ─────────────────────────────────────────────────────────────────────
# Raíz del proyecto (relativa a la ubicación del script). Para ejecución
# interactiva sin Rscript, asignar PROJECT_DIR manualmente.
.file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
PROJECT_DIR <- if (length(.file_arg))
  normalizePath(file.path(dirname(sub("^--file=", "", .file_arg)), "..")) else normalizePath("..")

QIIME_DIR    <- file.path(PROJECT_DIR, "qiime2_output")
META_FILE    <- file.path(PROJECT_DIR, "metadata", "sample-metadata-basal-6m.tsv")
MODELS_DIR   <- file.path(PROJECT_DIR, "R_output", "models")
OUT_DIR      <- file.path(PROJECT_DIR, "R_output", "models")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Importar desde QIIME2 ──────────────────────────────────────────────────────
message("Importando artefactos QIIME2...")

ps_raw <- qza_to_phyloseq(
  features = file.path(QIIME_DIR, "01_table.qza"),
  tree     = file.path(QIIME_DIR, "rooted-tree.qza"),
  taxonomy = file.path(QIIME_DIR, "taxonomy.qza"),
  metadata = META_FILE
)

message(sprintf("Phyloseq importado: %d muestras, %d taxa",
                nsamples(ps_raw), ntaxa(ps_raw)))

# ── Unir datos clínicos ────────────────────────────────────────────────────────
# Mapa: ID microbioma N01B2 → Pareja N01 → NEMO_01 en datos clínicos
clinical <- readRDS(file.path(MODELS_DIR, "df_analytic_fase1.rds"))

sdata <- sample_data(ps_raw) %>%
  as("data.frame") %>%
  rownames_to_column("sample_id") %>%
  mutate(dyad_id = paste0("NEMO_", sub("^N0*(\\d+)[BM]\\d$", "\\1", sample_id))) %>%
  left_join(
    clinical %>%
      select(dyad_id,
             WLZ_0, WLZ_6,
             GWG          = mat_weight_gain,
             alimentacion = feeding_6m,
             sexo         = inf_sex,
             Peso_Pre,
             overweight_risk),
    by = "dyad_id"
  ) %>%
  column_to_rownames("sample_id")

sample_data(ps_raw) <- sample_data(sdata)

# ── Filtrado de calidad ────────────────────────────────────────────────────────
# Eliminar controles y muestras con < 1000 reads
ps_filt <- ps_raw %>%
  prune_samples(sample_sums(.) >= 1000, .) %>%
  # Taxa presentes en al menos 10% de las muestras
  filter_taxa(function(x) sum(x > 0) >= 0.10 * nsamples(.), prune = TRUE)

message(sprintf("Tras filtrado: %d muestras, %d taxa",
                nsamples(ps_filt), ntaxa(ps_filt)))

# ── Agregar a nivel de género ─────────────────────────────────────────────────
ps_genus <- tax_glom(ps_filt, taxrank = "Genus", NArm = FALSE)
message(sprintf("Nivel género: %d géneros", ntaxa(ps_genus)))

# ── Dividir por grupo ─────────────────────────────────────────────────────────
ps_babies_basal <- subset_samples(ps_genus, Origen == "Bebe"  & Tiempo == "Basal")
ps_babies_6m    <- subset_samples(ps_genus, Origen == "Bebe"  & Tiempo == "mes_6")
ps_mothers      <- subset_samples(ps_genus, Origen == "Madre")

message(sprintf("Bebés basal: %d | Bebés 6m: %d | Madres: %d",
                nsamples(ps_babies_basal),
                nsamples(ps_babies_6m),
                nsamples(ps_mothers)))

# ── Guardar ───────────────────────────────────────────────────────────────────
saveRDS(ps_filt,         file.path(OUT_DIR, "ps_filtered.rds"))
saveRDS(ps_genus,        file.path(OUT_DIR, "ps_genus.rds"))
saveRDS(ps_babies_basal, file.path(OUT_DIR, "ps_babies_basal.rds"))
saveRDS(ps_babies_6m,    file.path(OUT_DIR, "ps_babies_6m.rds"))
saveRDS(ps_mothers,      file.path(OUT_DIR, "ps_mothers.rds"))

message("06_phyloseq_import.R completado. Objetos guardados en R_output/models/")
