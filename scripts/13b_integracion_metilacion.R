#!/usr/bin/env Rscript
# Integración exploratoria metilación (saliva), microbiota (heces) y adiposidad.


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
OUT_BASE    <- file.path(PROJECT_DIR, "R_output", "integracion_metilacion")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

# (1) Cargar metilación 
stopifnot(file.exists(MET_FILE))
meth <- read_excel(MET_FILE, sheet = "wide")
# La hoja "wide" nombra la primera columna como `x`; la renombramos a `Pareja`.
if (!"Pareja" %in% names(meth) && "x" %in% names(meth)) {
  meth <- dplyr::rename(meth, Pareja = x)
}
message(sprintf("Metilación: %d muestras × %d columnas",
                nrow(meth), ncol(meth)))

# Listado de genes (medias) y de CpGs individuales
gene_means <- grep("_mean$", names(meth), value = TRUE)
gene_cpgs  <- grep("_CpG\\d+$", names(meth), value = TRUE)
message(sprintf("  Genes (medias): %s", paste(gene_means, collapse = ", ")))

# (2) Cargar microbiota y residuos
ps <- readRDS(file.path(PROJECT_DIR, "R_output", "models", "ps_clean.rds"))
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
tt <- as.data.frame(tax_table(ps_genus), stringsAsFactors = FALSE)
genus_labels <- ifelse(is.na(tt$Genus) | tt$Genus == "",
                       paste0("UnknownG_", seq_len(nrow(tt))),
                       tt$Genus)
taxa_names(ps_genus) <- make.unique(genus_labels, sep = "_")

# Tabla de abundancias relativas por género
counts <- as(otu_table(ps_genus), "matrix")
if (taxa_are_rows(ps_genus)) counts <- t(counts)
rel <- sweep(counts, 1, rowSums(counts), "/")  # TSS por muestra
rel <- as.data.frame(rel)

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

# (3) Helpers 

# Mapear Origen_tiempo de microbiota a Tiempo de metilación
tiempo_map <- c(Madre_Basal = "m2", Bebe_Basal = "b2", Bebe_mes_6 = "b4")

# Prepara un data.frame con dos columnas: "micro" (abundancia relativa de un género) y "metilacion" (valor de metilación de un gen/CpG), alineados por Pareja, para un tiempo concreto de metilación.
prepare_pair <- function(micro_vec, micro_meta, gene_var, tiempo_meth, meth_df) {
  df_micro <- tibble(Pareja = micro_meta$Pareja, micro = micro_vec)
  df_meth  <- meth_df %>%
    filter(Tiempo == tiempo_meth) %>%
    select(Pareja, all_of(gene_var)) %>%
    rename(metilacion = all_of(gene_var))
  inner_join(df_micro, df_meth, by = "Pareja") %>%
    filter(!is.na(micro), !is.na(metilacion))
}

cor_one <- function(df) {
  if (nrow(df) < 6) return(c(rho = NA_real_, p = NA_real_, n = nrow(df)))
  ct <- suppressWarnings(cor.test(df$micro, df$metilacion, method = "spearman",
                                  exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value, n = nrow(df))
}

# Calcula la matriz Spearman entre un conjunto de features microbianas y los genes de metilación, dentro de un estrato concreto.
run_block <- function(feature_df,            # data.frame: cols = features, rows = muestras (alineado a micro_meta)
                      micro_meta,            # data.frame con $Pareja y $sample_id
                      strato_label,          # etiqueta "mothers_basal", "babies_6m", ...
                      tiempo_meth,           # tiempo en metilación a usar
                      vars_meth,             # vector con nombres de variables de metilación a probar
                      bloque,                # etiqueta de bloque ("B1_taxonomia", etc.)
                      tipo_feature,          # "genero" o "pathway"
                      meth_df) {             # tabla de metilación a usar
  out <- list()
  for (feat in colnames(feature_df)) {
    v_micro <- feature_df[[feat]]
    for (g in vars_meth) {
      df <- prepare_pair(v_micro, micro_meta, g, tiempo_meth, meth_df)
      r  <- cor_one(df)
      out[[length(out) + 1]] <- tibble(
        bloque       = bloque,
        estrato      = strato_label,
        tipo_feature = tipo_feature,
        feature      = feat,
        metilacion   = g,
        n            = r["n"],
        rho          = r["rho"],
        pval         = r["p"]
      )
    }
  }
  bind_rows(out)
}

# (4) Features priorizadas y preparación de datos comunes 
# Selección de géneros priorizados 
priorizados_madres <- c("GCA-900066575", "Sutterella", "Granulicatella",
                        "Turicibacter", "Lachnospira")
priorizados_b6m    <- c("Gemella", "UnknownG_442", "Clostridioides",
                        "Erysipelatoclostridium", "Eggerthella", "Blautia")

get_rel_subset <- function(rel_df, genus_vec, strato) {
  ids <- meta_micro$sample_id[meta_micro$Origen_tiempo == strato]
  ids <- intersect(ids, rownames(rel_df))
  available <- intersect(genus_vec, colnames(rel_df))
  missing   <- setdiff(genus_vec, available)
  if (length(missing) > 0) {
    message(sprintf("  [%s] géneros no encontrados: %s",
                    strato, paste(missing, collapse = ", ")))
  }
  list(
    feat = rel_df[ids, available, drop = FALSE],
    meta = meta_micro[ids, , drop = FALSE]
  )
}

# Microbiota funcional (MetaCyc)
PIC_DIR <- file.path(PROJECT_DIR, "qiime2_output", "picrust2_output")
pw_file <- file.path(PIC_DIR, "pathways_out", "path_abun_unstrat.tsv.gz")
prior_pw_madres <- c("FASYN-INITIAL-PWY", "PWY-5971", "PWY-5989",
                     "PWY-6282", "PWY-7664", "PWY0-862", "PWYG-321")
prior_pw_b2     <- c("RHAMCAT-PWY", "PWY-5741", "PWY-6731", "PWY-7347",
                     "SUCSYN-PWY", "PWY-722")
if (file.exists(pw_file)) {
  pw_raw <- read_tsv(pw_file, show_col_types = FALSE)
  pw_feat_names <- pw_raw[[1]]
  pw_mat <- as.data.frame(t(as.matrix(pw_raw[, -1])))
  colnames(pw_mat) <- pw_feat_names
  pw_mat <- pw_mat / rowSums(pw_mat)   # TSS
} else {
  warning("PICRUSt2 pathways no encontrado; se omite B2.")
  pw_mat <- NULL
}

# (5) Pipeline de correlaciones, parametrizado por la metilación 
run_all_blocks <- function(meth_df, gene_means) {

  # B1 · microbiota taxonómica × metilación (mismo individuo, mismo tiempo)
  b1 <- bind_rows(
    {
      s <- get_rel_subset(rel, priorizados_madres, "Madre_Basal")
      run_block(s$feat, s$meta, "mothers_basal", tiempo_map["Madre_Basal"],
                gene_means, "B1_taxonomia", "genero", meth_df)
    },
    {
      s <- get_rel_subset(rel, priorizados_b6m, "Bebe_mes_6")
      run_block(s$feat, s$meta, "babies_6m", tiempo_map["Bebe_mes_6"],
                gene_means, "B1_taxonomia", "genero", meth_df)
    }
  )

  # B2 · microbiota funcional (MetaCyc) × metilación
  b2 <- if (!is.null(pw_mat)) {
    bind_rows(
      {
        s <- get_rel_subset(pw_mat, prior_pw_madres, "Madre_Basal")
        run_block(s$feat, s$meta, "mothers_basal", tiempo_map["Madre_Basal"],
                  gene_means, "B2_funcional", "pathway", meth_df)
      },
      {
        s <- get_rel_subset(pw_mat, prior_pw_b2, "Bebe_Basal")
        run_block(s$feat, s$meta, "babies_basal", tiempo_map["Bebe_Basal"],
                  gene_means, "B2_funcional", "pathway", meth_df)
      }
    )
  } else tibble()

  # B3 · microbiota materna basal → metilación del bebé (eje vertical)
  run_vertical <- function(genus_vec, tiempo_bebe) {
    s <- get_rel_subset(rel, genus_vec, "Madre_Basal")
    run_block(s$feat, s$meta, paste0("madre_to_bebe_", tiempo_bebe),
              tiempo_bebe, gene_means, "B3_vertical", "genero", meth_df)
  }
  b3 <- bind_rows(
    run_vertical(priorizados_madres, "b2"),
    run_vertical(priorizados_madres, "b4")
  )

  # B4 · adiposidad × metilación, por tiempo — DOS enfoques en paralelo:
  #   B4_residuos   : residuos clínicos de Fase 1 × metilación
  #   B4_adiposidad : WLZ_6 observado del lactante × metilación (sin residuos)
  cor_meth_block <- function(adip_df, adip_var, bloque_lbl, tipo_lbl, feat_lbl) {
    map_df(c("m2", "b2", "b3", "b4"), function(tt) {
      df_meth_t <- meth_df %>% filter(Tiempo == tt) %>%
        select(Pareja, all_of(gene_means))
      joined <- inner_join(adip_df, df_meth_t, by = "Pareja")
      if (nrow(joined) < 6) return(tibble())
      map_df(gene_means, function(g) {
        v <- joined[[g]]
        r <- joined[[adip_var]]
        keep <- !is.na(v) & !is.na(r)
        if (sum(keep) < 6) return(tibble())
        ct <- suppressWarnings(cor.test(r[keep], v[keep],
                                        method = "spearman", exact = FALSE))
        tibble(bloque = bloque_lbl,
               estrato = paste0("metilacion_", tt),
               tipo_feature = tipo_lbl,
               feature = feat_lbl,
               metilacion = g,
               n = sum(keep),
               rho = unname(ct$estimate),
               pval = ct$p.value)
      })
    })
  }
  b4_res  <- cor_meth_block(residuals_df, "residuos", "B4_residuos",
                            "residuos_fase1", "residuos")
  b4_adip <- cor_meth_block(clin,         "WLZ_6",    "B4_adiposidad",
                            "WLZ_6", "WLZ_6")

  bind_rows(b1, b2, b3, b4_res, b4_adip) %>%
    filter(!is.na(pval)) %>%
    group_by(bloque) %>%
    mutate(qval = p.adjust(pval, method = "BH")) %>%
    ungroup() %>%
    arrange(bloque, qval, pval)
}

# (6) Ejecución principal: metilación CENSURADA 
cat("\n=== Correlaciones B1-B4 (metilación censurada) ===\n")
all_res <- run_all_blocks(meth, gene_means)

write_csv(all_res, file.path(OUT_BASE, "correlaciones_todas.csv"))

hits <- all_res %>% filter(qval < 0.25) %>% arrange(qval)
write_csv(hits, file.path(OUT_BASE, "hits_significativos.csv"))

# (9) Resumen consola 
cat("\n=== Resumen de hits significativos (q<0.25) por bloque ===\n")
print(all_res %>% filter(qval < 0.25) %>% count(bloque, estrato))

cat(sprintf("\nTotal correlaciones evaluadas: %d\n", nrow(all_res)))
cat(sprintf("Hits significativos (q<0.25): %d\n", nrow(hits)))
if (nrow(hits) > 0) {
  cat("\nTop 20:\n")
  print(hits %>% head(20))
}

#  (10) Sensibilidad a la censura de réplicas
cat("\n=== (10) Sensibilidad: metilación censurada vs sin censurar ===\n")

qc <- read_excel(MET_FILE, sheet = "qc_replicas")
# Reconstruye <gen>_mean SIN censura
meth_uncens <- qc %>%
  filter(!(Gen == "BCAT1" & CpG == "CpG4")) %>%
  group_by(Pareja, Tiempo, Gen) %>%
  summarise(valor = mean(media, na.rm = TRUE), .groups = "drop") %>%
  mutate(valor = ifelse(is.nan(valor), NA_real_, valor),
         col   = paste0(Gen, "_mean")) %>%
  select(Pareja, Tiempo, col, valor) %>%
  pivot_wider(names_from = col, values_from = valor)

all_res_uncens <- run_all_blocks(meth_uncens, gene_means)

key <- c("bloque", "estrato", "tipo_feature", "feature", "metilacion")
sens <- full_join(
  all_res        %>% select(all_of(key), n, rho, pval, qval),
  all_res_uncens %>% select(all_of(key), n, rho, pval, qval),
  by = key, suffix = c("_cens", "_uncens")
) %>%
  mutate(
    hit_cens   = !is.na(qval_cens)   & qval_cens   < 0.25,
    hit_uncens = !is.na(qval_uncens) & qval_uncens < 0.25,
    d_rho      = rho_cens - rho_uncens,
    cambia_hit = hit_cens != hit_uncens
  ) %>%
  arrange(desc(cambia_hit), pmin(qval_cens, qval_uncens, na.rm = TRUE))

write_csv(sens, file.path(OUT_BASE, "sensibilidad_censura_replicas.csv"))

cat(sprintf("Correlaciones comparadas        : %d\n", nrow(sens)))
cat(sprintf("Hits (q<0.25) solo CON censura   : %d\n",
            sum(sens$hit_cens & !sens$hit_uncens, na.rm = TRUE)))
cat(sprintf("Hits (q<0.25) solo SIN censura   : %d\n",
            sum(!sens$hit_cens & sens$hit_uncens, na.rm = TRUE)))
cat(sprintf("Hits que CAMBIAN de estado       : %d\n",
            sum(sens$cambia_hit, na.rm = TRUE)))
cat(sprintf("Máx |Δrho| entre versiones       : %.3f\n",
            max(abs(sens$d_rho), na.rm = TRUE)))
if (any(sens$cambia_hit, na.rm = TRUE)) {
  cat("\nCorrelaciones que cambian de hit al censurar:\n")
  print(sens %>% filter(cambia_hit) %>%
        select(all_of(key), rho_cens, qval_cens, rho_uncens, qval_uncens))
} else {
  cat("Ninguna correlación cambia de hit: la censura NO altera las conclusiones.\n")
}

cat("\nFicheros generados:\n")
cat(sprintf("  %s\n", file.path(OUT_BASE, "correlaciones_todas.csv")))
cat(sprintf("  %s\n", file.path(OUT_BASE, "hits_significativos.csv")))
cat(sprintf("  %s\n", file.path(OUT_BASE, "sensibilidad_censura_replicas.csv")))
cat("\n13b_integracion_metilacion.R completado.\n")
