#!/usr/bin/env Rscript
# Análisis de transmisión vertical madre → hijo.

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(tidyverse)
})

# Rutas
.file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
PROJECT_DIR <- if (length(.file_arg))
  normalizePath(file.path(dirname(sub("^--file=", "", .file_arg)), "..")) else normalizePath("..")
OUT_BASE    <- file.path(PROJECT_DIR, "R_output", "transmision_vertical")
FIG_DIR     <- file.path(OUT_BASE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

ps <- readRDS(file.path(PROJECT_DIR, "R_output", "models", "ps_clean.rds"))
message(sprintf("Phyloseq cargado: %d ASVs × %d muestras",
                ntaxa(ps), nsamples(ps)))

meta <- as(sample_data(ps), "data.frame")
meta$sample_id <- rownames(meta)

# Validación
required_cols <- c("Pareja", "Origen", "Origen_tiempo")
stopifnot(all(required_cols %in% colnames(meta)))

# (1) Función: distancias intra vs inter para un par de tiempos 
calc_intra_inter <- function(ps_obj, meta_obj, tiempo_madre, tiempo_bebe,
                              metric_name, dist_obj = NULL) {
  ids_m <- meta_obj$sample_id[meta_obj$Origen_tiempo == tiempo_madre]
  ids_b <- meta_obj$sample_id[meta_obj$Origen_tiempo == tiempo_bebe]

  parejas_m <- setNames(meta_obj$Pareja[match(ids_m, meta_obj$sample_id)], ids_m)
  parejas_b <- setNames(meta_obj$Pareja[match(ids_b, meta_obj$sample_id)], ids_b)

  if (is.null(dist_obj)) {
    if (metric_name == "bray") {
      dist_obj <- phyloseq::distance(ps_obj, method = "bray")
    } else if (metric_name == "wunifrac") {
      dist_obj <- phyloseq::distance(ps_obj, method = "wunifrac")
    } else if (metric_name == "unifrac") {
      dist_obj <- phyloseq::distance(ps_obj, method = "unifrac")
    } else stop("Métrica no reconocida: ", metric_name)
  }
  D <- as.matrix(dist_obj)

  records <- list()
  for (b_id in ids_b) {
    pareja <- parejas_b[b_id]
    own_m  <- names(parejas_m)[parejas_m == pareja]
    other_m <- names(parejas_m)[parejas_m != pareja]
    if (length(own_m) == 0) next
    own_d   <- mean(D[b_id, own_m])
    other_d <- mean(D[b_id, other_m])
    records[[length(records) + 1]] <- tibble(
      bebe_id      = b_id,
      Pareja       = pareja,
      tiempo_bebe  = tiempo_bebe,
      tiempo_madre = tiempo_madre,
      metric       = metric_name,
      dist_intra   = own_d,
      dist_inter   = other_d
    )
  }
  bind_rows(records)
}

# (2) Pre-cálculo de matrices de distancia (todas las muestras necesarias)
# Subconjunto: madres basal + bebés basal + bebés 6m
ps_sub <- subset_samples(ps, Origen_tiempo %in%
                           c("Madre_Basal", "Bebe_Basal", "Bebe_mes_6"))
meta_sub <- as(sample_data(ps_sub), "data.frame")
meta_sub$sample_id <- rownames(meta_sub)
message(sprintf("Subset trabajo: %d muestras", nsamples(ps_sub)))
print(table(meta_sub$Origen_tiempo))

message("\nCalculando matrices de distancia (puede tardar 1-2 min)...")
d_bray  <- phyloseq::distance(ps_sub, method = "bray")
d_wuf   <- phyloseq::distance(ps_sub, method = "wunifrac")
d_uuf   <- phyloseq::distance(ps_sub, method = "unifrac")

# (3) A1: intra vs inter para los dos tiempos del bebé × 3 métricas 
message("\n=== A1: distancias intra-díada vs inter-díada ===")
A1 <- bind_rows(
  calc_intra_inter(ps_sub, meta_sub, "Madre_Basal", "Bebe_Basal", "bray",     d_bray),
  calc_intra_inter(ps_sub, meta_sub, "Madre_Basal", "Bebe_Basal", "wunifrac", d_wuf),
  calc_intra_inter(ps_sub, meta_sub, "Madre_Basal", "Bebe_Basal", "unifrac",  d_uuf),
  calc_intra_inter(ps_sub, meta_sub, "Madre_Basal", "Bebe_mes_6", "bray",     d_bray),
  calc_intra_inter(ps_sub, meta_sub, "Madre_Basal", "Bebe_mes_6", "wunifrac", d_wuf),
  calc_intra_inter(ps_sub, meta_sub, "Madre_Basal", "Bebe_mes_6", "unifrac",  d_uuf)
)

# Wilcoxon pareado intra vs inter por (tiempo_bebe × metric)
tests_A1 <- A1 %>%
  group_by(tiempo_bebe, metric) %>%
  summarise(
    n            = n(),
    media_intra  = mean(dist_intra),
    media_inter  = mean(dist_inter),
    delta        = mean(dist_intra - dist_inter),  # negativo = madre propia más cercana
    pval_wilcox  = wilcox.test(dist_intra, dist_inter, paired = TRUE,
                                alternative = "less")$p.value,
    .groups = "drop"
  ) %>%
  mutate(qval = p.adjust(pval_wilcox, method = "BH"))

write_csv(A1,        file.path(OUT_BASE, "distancias_intra_vs_inter.csv"))
write_csv(tests_A1,  file.path(OUT_BASE, "tests_A1_wilcoxon.csv"))

cat("\nResultados Wilcoxon pareado (alternative = less; intra < inter):\n")
print(tests_A1)

# (4) A2: compartición de ASVs y de géneros 
message("\n=== A2: compartición de ASVs y géneros madre→hijo ===")

calc_sharing <- function(ps_obj, meta_obj, tiempo_bebe, level = "ASV") {
  if (level == "Genus") {
    ps_obj <- tax_glom(ps_obj, taxrank = "Genus", NArm = FALSE)
  }
  mat <- as(otu_table(ps_obj), "matrix")
  if (taxa_are_rows(ps_obj)) mat <- t(mat)
  pa <- mat > 0  # presencia/ausencia

  ids_m <- meta_obj$sample_id[meta_obj$Origen_tiempo == "Madre_Basal"]
  ids_b <- meta_obj$sample_id[meta_obj$Origen_tiempo == tiempo_bebe]
  parejas_m <- setNames(meta_obj$Pareja[match(ids_m, meta_obj$sample_id)], ids_m)
  parejas_b <- setNames(meta_obj$Pareja[match(ids_b, meta_obj$sample_id)], ids_b)

  records <- list()
  for (b_id in ids_b) {
    pareja  <- parejas_b[b_id]
    own_m   <- names(parejas_m)[parejas_m == pareja]
    other_m <- names(parejas_m)[parejas_m != pareja]
    if (length(own_m) == 0) next
    taxa_b      <- which(pa[b_id, ])
    if (length(taxa_b) == 0) next
    taxa_own_m  <- which(pa[own_m[1], ])
    pct_own     <- length(intersect(taxa_b, taxa_own_m)) / length(taxa_b)
    # Background: media del % compartido con cada una de las otras madres
    pct_other_v <- sapply(other_m, function(m) {
      length(intersect(taxa_b, which(pa[m, ]))) / length(taxa_b)
    })
    pct_other <- mean(pct_other_v, na.rm = TRUE)
    records[[length(records) + 1]] <- tibble(
      bebe_id      = b_id,
      Pareja       = pareja,
      tiempo_bebe  = tiempo_bebe,
      level        = level,
      n_taxa_bebe  = length(taxa_b),
      pct_compartido_propia = pct_own * 100,
      pct_compartido_otras  = pct_other * 100
    )
  }
  bind_rows(records)
}

A2 <- bind_rows(
  calc_sharing(ps_sub, meta_sub, "Bebe_Basal", "ASV"),
  calc_sharing(ps_sub, meta_sub, "Bebe_Basal", "Genus"),
  calc_sharing(ps_sub, meta_sub, "Bebe_mes_6", "ASV"),
  calc_sharing(ps_sub, meta_sub, "Bebe_mes_6", "Genus")
)

tests_A2 <- A2 %>%
  group_by(tiempo_bebe, level) %>%
  summarise(
    n              = n(),
    media_propia   = mean(pct_compartido_propia),
    media_otras    = mean(pct_compartido_otras),
    delta_pct      = mean(pct_compartido_propia - pct_compartido_otras),
    pval_wilcox    = wilcox.test(pct_compartido_propia, pct_compartido_otras,
                                  paired = TRUE, alternative = "greater")$p.value,
    .groups = "drop"
  ) %>%
  mutate(qval = p.adjust(pval_wilcox, method = "BH"))

write_csv(A2,       file.path(OUT_BASE, "asv_sharing.csv"))
write_csv(tests_A2, file.path(OUT_BASE, "tests_A2_sharing.csv"))

cat("\nCompartición de taxa (Wilcoxon pareado; propia > otras):\n")
print(tests_A2)

# (5) A3: ASVs candidatos a transmisión vertical 
message("\n=== A3: ASVs candidatos a transmisión vertical ===")

find_vertical_candidates <- function(ps_obj, meta_obj, tiempo_bebe,
                                      min_dyad_prevalence = 0.30) {
  mat <- as(otu_table(ps_obj), "matrix")
  if (taxa_are_rows(ps_obj)) mat <- t(mat)
  pa <- mat > 0

  ids_m <- meta_obj$sample_id[meta_obj$Origen_tiempo == "Madre_Basal"]
  ids_b <- meta_obj$sample_id[meta_obj$Origen_tiempo == tiempo_bebe]
  parejas_m <- setNames(meta_obj$Pareja[match(ids_m, meta_obj$sample_id)], ids_m)
  parejas_b <- setNames(meta_obj$Pareja[match(ids_b, meta_obj$sample_id)], ids_b)

  # díadas con AMBAS muestras disponibles
  parejas_completas <- intersect(parejas_m, parejas_b)
  n_dyads <- length(parejas_completas)

  results <- map_df(taxa_names(ps_obj), function(asv) {
    n_madre <- 0; n_bebe <- 0; n_dyad_both <- 0
    for (p in parejas_completas) {
      m_id <- names(parejas_m)[parejas_m == p][1]
      b_id <- names(parejas_b)[parejas_b == p][1]
      in_m <- pa[m_id, asv]
      in_b <- pa[b_id, asv]
      if (in_m) n_madre <- n_madre + 1
      if (in_b) n_bebe  <- n_bebe + 1
      if (in_m && in_b) n_dyad_both <- n_dyad_both + 1
    }
    tibble(
      ASV         = asv,
      n_dyads     = n_dyads,
      n_madre     = n_madre,
      n_bebe      = n_bebe,
      n_compartido = n_dyad_both,
      prev_compartida = n_dyad_both / n_dyads
    )
  })

  # Anotamos taxonomía
  tt <- as.data.frame(tax_table(ps_obj), stringsAsFactors = FALSE)
  results$Phylum <- tt[results$ASV, "Phylum"]
  results$Family <- tt[results$ASV, "Family"]
  results$Genus  <- tt[results$ASV, "Genus"]
  results$tiempo_bebe <- tiempo_bebe

  results %>%
    filter(prev_compartida >= min_dyad_prevalence) %>%
    arrange(desc(prev_compartida))
}

A3 <- bind_rows(
  find_vertical_candidates(ps_sub, meta_sub, "Bebe_Basal", min_dyad_prevalence = 0.30),
  find_vertical_candidates(ps_sub, meta_sub, "Bebe_mes_6", min_dyad_prevalence = 0.30)
)
write_csv(A3, file.path(OUT_BASE, "asvs_candidatos_transmision.csv"))

cat(sprintf("\nASVs candidatos a transmisión (compartidos en ≥30 %% de las díadas):\n"))
cat(sprintf("  Bebé basal: %d ASVs\n",
            sum(A3$tiempo_bebe == "Bebe_Basal")))
cat(sprintf("  Bebé 6m   : %d ASVs\n",
            sum(A3$tiempo_bebe == "Bebe_mes_6")))

# (6) Figuras 
theme_paper <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = NA),
        plot.title.position = "plot")

# Fig A: boxplots intra vs inter
plot_A1 <- A1 %>%
  pivot_longer(c(dist_intra, dist_inter),
               names_to = "tipo", values_to = "distancia") %>%
  mutate(tipo = recode(tipo,
                       dist_intra = "Madre propia",
                       dist_inter = "Otras madres"),
         metric = factor(metric, levels = c("bray", "wunifrac", "unifrac"),
                         labels = c("Bray-Curtis",
                                    "UniFrac ponderada",
                                    "UniFrac no ponderada")),
         tiempo_bebe = factor(tiempo_bebe,
                              levels = c("Bebe_Basal", "Bebe_mes_6"),
                              labels = c("Meconio (basal)", "6 meses")))

p_dist <- ggplot(plot_A1, aes(tipo, distancia, fill = tipo)) +
  geom_boxplot(width = 0.55, outlier.size = 0.8, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 0.7, alpha = 0.45) +
  facet_grid(metric ~ tiempo_bebe, scales = "free_y") +
  scale_fill_manual(values = c("Madre propia" = "#A23B72",
                                "Otras madres" = "#E0E0E0")) +
  labs(title = "Distancias microbianas: madre propia vs otras madres",
       subtitle = "Bebés más próximos a su madre → evidencia de transmisión vertical",
       x = NULL, y = "Distancia") +
  theme_paper + theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "fig_A1_distancias_intra_vs_inter.pdf"),
       p_dist, width = 8, height = 8)

# Fig B: % compartido propia vs otras
plot_A2 <- A2 %>%
  pivot_longer(c(pct_compartido_propia, pct_compartido_otras),
               names_to = "tipo", values_to = "pct") %>%
  mutate(tipo = recode(tipo,
                       pct_compartido_propia = "Madre propia",
                       pct_compartido_otras  = "Otras madres"),
         level = factor(level, levels = c("ASV", "Genus")),
         tiempo_bebe = factor(tiempo_bebe,
                              levels = c("Bebe_Basal", "Bebe_mes_6"),
                              labels = c("Meconio (basal)", "6 meses")))

p_share <- ggplot(plot_A2, aes(tipo, pct, fill = tipo)) +
  geom_boxplot(width = 0.55, outlier.size = 0.8, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 0.7, alpha = 0.45) +
  facet_grid(level ~ tiempo_bebe, scales = "free_y") +
  scale_fill_manual(values = c("Madre propia" = "#2E86AB",
                                "Otras madres" = "#E0E0E0")) +
  labs(title = "Compartición de taxa madre→hijo",
       subtitle = "% de taxa del bebé presentes también en la madre",
       x = NULL, y = "% de taxa compartidos") +
  theme_paper + theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "fig_A2_sharing.pdf"),
       p_share, width = 7.5, height = 7)

cat("\nFicheros generados:\n")
for (f in list.files(OUT_BASE, recursive = TRUE, full.names = TRUE)) {
  cat(sprintf("  %s\n", f))
}
cat("\n12_transmision_vertical.R completado.\n")
