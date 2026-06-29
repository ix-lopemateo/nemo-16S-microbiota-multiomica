#!/usr/bin/env Rscript
# Análisis descriptivo de composición taxonómica.
# Importa los artefactos QIIME2 a phyloseq, calcula abundancias relativas y genera figuras y tablas de composición.

suppressPackageStartupMessages({
  library(qiime2R)
  library(phyloseq)
  library(tidyverse)
  library(scales)
  library(RColorBrewer)
  library(viridis)
})

# Rutas
.file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
PROJECT_DIR <- if (length(.file_arg))
  normalizePath(file.path(dirname(sub("^--file=", "", .file_arg)), "..")) else normalizePath("..")
QIIME_DIR     <- file.path(PROJECT_DIR, "qiime2_output")
OUT_FIG       <- file.path(PROJECT_DIR, "R_output", "figures")
OUT_TAB       <- file.path(PROJECT_DIR, "R_output", "tables")
OUT_MOD       <- file.path(PROJECT_DIR, "R_output", "models")
for (d in c(OUT_FIG, OUT_TAB, OUT_MOD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Importar a phyloseq (sin metadata) y limpiar taxonomía
message("Importando artefactos QIIME2 a phyloseq...")
ps <- qza_to_phyloseq(
  features = file.path(QIIME_DIR, "04_table.qza"),
  tree     = file.path(QIIME_DIR, "05_rooted-tree.qza"),
  taxonomy = file.path(QIIME_DIR, "03_taxonomy.qza")
)

# Cargar y limpiar el metadata manualmente
meta_file <- file.path(PROJECT_DIR, "metadata", "sample-metadata-basal-6m.tsv")
meta <- read.delim(meta_file, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE,
                   comment.char = "")
# Descartar fila opcional "#q2:types" si está presente
meta <- meta[!grepl("^#", meta[[1]]) | grepl("^#SampleID", meta[[1]]), ]
meta <- meta[!grepl("^#q2:types", meta[[1]]), ]
# Renombrar primera columna a SampleID
names(meta)[1] <- "SampleID"
rownames(meta) <- meta$SampleID
# Inyectar en phyloseq, conservando solo las muestras comunes
common_samples <- intersect(rownames(meta), sample_names(ps))
ps <- prune_samples(common_samples, ps)
meta <- meta[common_samples, , drop = FALSE]
sample_data(ps) <- sample_data(meta)

message(sprintf("Phyloseq importado: %d muestras, %d ASVs, %d niveles taxonómicos",
                nsamples(ps), ntaxa(ps), length(rank_names(ps))))
message("Columnas de metadata: ", paste(colnames(sample_data(ps)), collapse = ", "))

# Limpiar prefijos taxonómicos (d__, p__, ...)
tax_clean <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
tax_clean[] <- lapply(tax_clean, function(x) gsub("^[a-z]__", "", x))
tax_table(ps) <- as.matrix(tax_clean)

# Guardar objeto phyloseq limpio para análisis posteriores
saveRDS(ps, file.path(OUT_MOD, "ps_clean.rds"))
message("Phyloseq guardado en R_output/models/ps_clean.rds")

# Normalizar a abundancia relativa 
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Orden de los niveles del factor Origen_tiempo para todas las figuras
sd <- as(sample_data(ps_rel), "data.frame")
sd$Origen_tiempo <- factor(sd$Origen_tiempo,
                            levels = c("Bebe_Basal", "Bebe_mes_6", "Madre_Basal"))
sample_data(ps_rel) <- sample_data(sd)

# Figura 1: stacked barplot por Filo (por muestra)
message("Generando barplot de Filos...")

ps_phylum <- tax_glom(ps_rel, taxrank = "Phylum", NArm = FALSE)
df_phylum <- psmelt(ps_phylum) %>%
  mutate(Phylum = ifelse(is.na(Phylum) | Phylum == "", "Sin asignación", Phylum))

top_phyla <- df_phylum %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 8) %>%
  pull(Phylum)

df_phylum <- df_phylum %>%
  mutate(Phylum_top = ifelse(Phylum %in% top_phyla, Phylum, "Otros"),
         Phylum_top = factor(Phylum_top, levels = c(top_phyla, "Otros")))

p_phylum <- ggplot(df_phylum, aes(x = Sample, y = Abundance, fill = Phylum_top)) +
  geom_col(width = 1) +
  facet_grid(~ Origen_tiempo, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set2", name = "Filo") +
  labs(x = NULL, y = "Abundancia relativa") +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        strip.background = element_rect(fill = "grey95"),
        panel.grid = element_blank())

ggsave(file.path(OUT_FIG, "composicion_phylum_barplot.pdf"),
       p_phylum, width = 13, height = 5, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "composicion_phylum_barplot.png"),
       p_phylum, width = 13, height = 5, dpi = 300)

# Figura 2: stacked barplot agregado por grupo (Género top 15)
message("Generando barplot de Géneros top 15...")

ps_genus <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)
df_genus <- psmelt(ps_genus) %>%
  mutate(Genus  = ifelse(is.na(Genus)  | Genus  == "", "Sin asignación", Genus),
         Family = ifelse(is.na(Family) | Family == "", "Sin asignación", Family))

top_genera <- df_genus %>%
  filter(Genus != "Sin asignación") %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

df_genus_grp <- df_genus %>%
  mutate(Genus_top = ifelse(Genus %in% top_genera, Genus, "Otros"),
         Genus_top = factor(Genus_top, levels = c(top_genera, "Otros"))) %>%
  group_by(Origen_tiempo, Genus_top) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop")

palette_15 <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set3")[1:7], "grey75")

p_genus <- ggplot(df_genus_grp,
                  aes(x = Origen_tiempo, y = mean_abund, fill = Genus_top)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = palette_15, name = "Género") +
  labs(x = NULL, y = "Abundancia relativa media") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank())

ggsave(file.path(OUT_FIG, "composicion_genus_barplot.pdf"),
       p_genus, width = 8, height = 6, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "composicion_genus_barplot.png"),
       p_genus, width = 8, height = 6, dpi = 300)

# Figura 3: boxplots top 15 géneros × Origen_tiempo
message("Generando boxplots de géneros top 15...")

df_top <- df_genus %>%
  filter(Genus %in% top_genera) %>%
  mutate(Genus = factor(Genus, levels = top_genera),
         pct   = Abundance * 100)

p_box <- ggplot(df_top, aes(x = Origen_tiempo, y = pct, fill = Origen_tiempo)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 5) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(x = NULL, y = "Abundancia relativa (%)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1),
        strip.text   = element_text(face = "italic", size = 8),
        legend.position = "none")

ggsave(file.path(OUT_FIG, "composicion_genus_boxplots.pdf"),
       p_box, width = 12, height = 8, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "composicion_genus_boxplots.png"),
       p_box, width = 12, height = 8, dpi = 300)

# Figura 4: heatmap top 25 géneros por muestra
message("Generando heatmap de géneros top 25...")

top25 <- df_genus %>%
  filter(Genus != "Sin asignación") %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 25) %>%
  pull(Genus)

df_hm <- df_genus %>%
  filter(Genus %in% top25) %>%
  mutate(abund_log = log10(Abundance + 1e-5),
         Genus = factor(Genus, levels = rev(top25)))

p_hm <- ggplot(df_hm, aes(x = Sample, y = Genus, fill = abund_log)) +
  geom_tile() +
  facet_grid(~ Origen_tiempo, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(option = "magma",
                       name = expression(log[10](abund))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        axis.text.y      = element_text(face = "italic", size = 8),
        strip.background = element_rect(fill = "grey95"),
        panel.grid       = element_blank())

ggsave(file.path(OUT_FIG, "composicion_genus_heatmap.pdf"),
       p_hm, width = 13, height = 7, device = cairo_pdf)
ggsave(file.path(OUT_FIG, "composicion_genus_heatmap.png"),
       p_hm, width = 13, height = 7, dpi = 300)

# Tablas de resumen
message("Generando tablas de abundancia y prevalencia por grupo...")

tab_abund <- df_genus %>%
  group_by(Origen_tiempo, Genus) %>%
  summarise(mean_pct = mean(Abundance) * 100, .groups = "drop") %>%
  pivot_wider(names_from = Origen_tiempo, values_from = mean_pct, values_fill = 0) %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  ungroup() %>%
  arrange(desc(total)) %>%
  select(-total)

write_csv(tab_abund, file.path(OUT_TAB, "abundancias_genero_por_grupo.csv"))

tab_prev <- df_genus %>%
  group_by(Origen_tiempo, Genus) %>%
  summarise(prev_pct = mean(Abundance > 0) * 100, .groups = "drop") %>%
  pivot_wider(names_from = Origen_tiempo, values_from = prev_pct, values_fill = 0)

write_csv(tab_prev, file.path(OUT_TAB, "prevalencia_genero_por_grupo.csv"))

message("\n08_composicion.R completado.")
message("Figuras: ", OUT_FIG)
message("Tablas:  ", OUT_TAB)
message("Phyloseq RDS: ", file.path(OUT_MOD, "ps_clean.rds"))
