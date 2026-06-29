#!/usr/bin/env Rscript
# Ensamblado de features de microbiota, funcionales y epigenéticas en formato wide
suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(phyloseq)
  library(vegan)
  library(tidyverse)
})

# Rutas
.file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
PROJECT_DIR <- if (length(.file_arg))
  normalizePath(file.path(dirname(sub("^--file=", "", .file_arg)), "..")) else normalizePath("..")
DATA_DIR    <- file.path(PROJECT_DIR, "data")
RES_DIR     <- file.path(PROJECT_DIR, "R_output")

WIDE_IN     <- file.path(DATA_DIR, "nemo_wide_clean_con_leyenda.xlsx")
WIDE_OUT    <- file.path(DATA_DIR, "nemo_wide_clean_con_leyenda_multiomics.xlsx")
MET_FILE    <- file.path(DATA_DIR, "pirosecuenciacion_metilacion_crudo_limpio.xlsx")
PS_FILE     <- file.path(RES_DIR, "models", "ps_clean.rds")
COR_FILE    <- file.path(RES_DIR, "integracion_metilacion", "correlaciones_todas.csv")
PIC_FILE    <- file.path(PROJECT_DIR, "qiime2_output", "picrust2_output",
                         "pathways_out", "path_abun_unstrat.tsv.gz")
VT_DIR      <- file.path(RES_DIR, "transmision_vertical")

stopifnot(file.exists(WIDE_IN), file.exists(MET_FILE), file.exists(PS_FILE))

# Mapa de estrato microbiota -> sufijo de columna
estr_suf <- c(Madre_Basal = "mad", Bebe_Basal = "bb", Bebe_mes_6 = "b6m")
# Mapa de estrato -> tiempo de metilación
estr_met <- c(Madre_Basal = "m2", Bebe_Basal = "b2", Bebe_mes_6 = "b4")

# dyad_id (NEMO_01) -> Pareja (N01) usado en microbiota y metilación
dyad_to_pareja <- function(x) sub("^NEMO_", "N", x)

sanitize <- function(x) gsub("(^_|_$)", "", gsub("[^A-Za-z0-9]+", "_", x))

# (1) Microbiota: género glomerado + etiquetas
ps <- readRDS(PS_FILE)
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
tt <- as.data.frame(tax_table(ps_genus), stringsAsFactors = FALSE)
genus_labels <- ifelse(is.na(tt$Genus) | tt$Genus == "",
                       paste0("UnknownG_", seq_len(nrow(tt))), tt$Genus)
taxa_names(ps_genus) <- make.unique(genus_labels, sep = "_")

counts <- as(otu_table(ps_genus), "matrix")
if (taxa_are_rows(ps_genus)) counts <- t(counts)        # muestras x géneros

meta <- as(sample_data(ps_genus), "data.frame")
meta$sample_id <- rownames(meta)

# CLR por muestra (pseudoconteo 0.5 sobre conteos) 
clr_mat <- {
  x <- counts + 0.5
  logx <- log(x)
  sweep(logx, 1, rowMeans(logx), "-")
}
clr_df <- as.data.frame(clr_mat)

# Diversidad alfa sobre conteos crudos
div_df <- data.frame(
  sample_id = rownames(counts),
  shannon   = vegan::diversity(counts, index = "shannon"),
  observed  = vegan::specnumber(counts),
  row.names = NULL
)

# Géneros priorizados por estrato
priorizados_madres <- c("GCA-900066575", "Sutterella", "Granulicatella",
                        "Turicibacter", "Lachnospira")
priorizados_b6m    <- c("Gemella", "UnknownG_442", "Clostridioides",
                        "Erysipelatoclostridium", "Eggerthella", "Blautia")

# Extrae, para un estrato, las columnas CLR de los géneros pedidos por Pareja
get_taxa_dyad <- function(feature_df, genus_vec, estrato) {
  ids <- meta$sample_id[meta$Origen_tiempo == estrato]
  avail <- intersect(genus_vec, colnames(feature_df))
  miss  <- setdiff(genus_vec, colnames(feature_df))
  if (length(miss)) message(sprintf("  [%s] géneros no hallados: %s",
                                    estrato, paste(miss, collapse = ", ")))
  suf <- estr_suf[[estrato]]
  out <- feature_df[ids, avail, drop = FALSE]
  out$Pareja <- meta$Pareja[match(ids, meta$sample_id)]
  names(out)[names(out) != "Pareja"] <-
    paste0("mb_", suf, "_", sanitize(avail))
  out
}

tax_mad <- get_taxa_dyad(clr_df, priorizados_madres, "Madre_Basal")
tax_b6m <- get_taxa_dyad(clr_df, priorizados_b6m,    "Bebe_mes_6")

# Diversidad alfa por estrato 
get_div_dyad <- function(estrato) {
  ids <- meta$sample_id[meta$Origen_tiempo == estrato]
  suf <- estr_suf[[estrato]]
  d <- div_df[match(ids, div_df$sample_id), c("shannon", "observed")]
  names(d) <- paste0("div_", suf, "_", names(d))
  d$Pareja <- meta$Pareja[match(ids, meta$sample_id)]
  d
}
div_mad <- get_div_dyad("Madre_Basal")
div_bb  <- get_div_dyad("Bebe_Basal")
div_b6m <- get_div_dyad("Bebe_mes_6")

# (2) Vías funcionales (PICRUSt2) en madre basal, CLR 
prior_pw_madres <- c("FASYN-INITIAL-PWY", "PWY-5971", "PWY-5989",
                     "PWY-6282", "PWY-7664", "PWY0-862", "PWYG-321")
pwy_mad <- NULL
if (file.exists(PIC_FILE)) {
  pw_raw <- read_tsv(PIC_FILE, show_col_types = FALSE)
  pw_names <- pw_raw[[1]]
  pw_mat <- as.matrix(pw_raw[, -1])              # vías x muestras
  pw_mat <- t(pw_mat)                            # muestras x vías
  colnames(pw_mat) <- pw_names
  pw_clr <- { x <- pw_mat + 0.5; lx <- log(x); sweep(lx, 1, rowMeans(lx), "-") }
  ids <- meta$sample_id[meta$Origen_tiempo == "Madre_Basal"]
  ids <- intersect(ids, rownames(pw_clr))
  avail <- intersect(prior_pw_madres, colnames(pw_clr))
  if (length(avail)) {
    pwy_mad <- as.data.frame(pw_clr[ids, avail, drop = FALSE])
    names(pwy_mad) <- paste0("pwy_mad_", sanitize(avail))
    pwy_mad$Pareja <- meta$Pareja[match(ids, meta$sample_id)]
  }
} else {
  warning("PICRUSt2 no encontrado; se omiten las vías funcionales.")
}

# (3) Transmisión vertical madre→bebé (métricas por díada) 
vt_list <- list()
sh_file <- file.path(VT_DIR, "asv_sharing.csv")
di_file <- file.path(VT_DIR, "distancias_intra_vs_inter.csv")
if (file.exists(sh_file)) {
  sh <- read_csv(sh_file, show_col_types = FALSE) %>%
    filter(level == "ASV") %>%
    select(Pareja, tiempo_bebe, pct_compartido_propia) %>%
    mutate(estr = recode(tiempo_bebe, "Bebe_Basal" = "bb", "Bebe_mes_6" = "b6m")) %>%
    transmute(Pareja, var = paste0("vt_", estr, "_pct_shared"),
              val = pct_compartido_propia) %>%
    pivot_wider(names_from = var, values_from = val)
  vt_list[["sharing"]] <- sh
}
if (file.exists(di_file)) {
  di <- read_csv(di_file, show_col_types = FALSE) %>%
    filter(metric == "bray") %>%
    select(Pareja, tiempo_bebe, dist_intra) %>%
    mutate(estr = recode(tiempo_bebe, "Bebe_Basal" = "bb", "Bebe_mes_6" = "b6m")) %>%
    transmute(Pareja, var = paste0("vt_", estr, "_bray_madre"),
              val = dist_intra) %>%
    pivot_wider(names_from = var, values_from = val)
  vt_list[["dist"]] <- di
}
vt_df <- if (length(vt_list)) reduce(vt_list, full_join, by = "Pareja") else NULL

# (4) Epigenética: genes implicados por estrato (|rho|>=0.30 en B1) 
meth <- read_excel(MET_FILE, sheet = "wide")

if (!"Pareja" %in% names(meth) && "x" %in% names(meth))
  meth <- dplyr::rename(meth, Pareja = x)

implicated <- list(Madre_Basal = character(0), Bebe_mes_6 = character(0))
if (file.exists(COR_FILE)) {
  cors <- read_csv(COR_FILE, show_col_types = FALSE)
  estr_cor <- c(Madre_Basal = "mothers_basal", Bebe_mes_6 = "babies_6m")
  for (e in names(estr_cor)) {
    g <- cors %>%
      filter(bloque == "B1_taxonomia", estrato == estr_cor[[e]],
             abs(rho) >= 0.30) %>%
      pull(metilacion) %>% unique()
    implicated[[e]] <- g
  }
}
# Forzar genes con señal robusta detectada FUERA del bloque B1, que el criterio
# automático (solo B1_taxonomia) no captaría:
#  - BCAT1: la correlación más fuerte de toda la integración (B2_funcional ×
#    metilación, q ≈ 0.11, única que sobrevive FDR junto a las vías de ácidos
#    grasos maternas). Se añade en madre basal (señal principal) y bebé 6 m.
forced_meth <- list(Madre_Basal = "BCAT1_mean", Bebe_mes_6 = "BCAT1_mean")
for (e in names(forced_meth))
  implicated[[e]] <- union(implicated[[e]], forced_meth[[e]])

message("Genes de metilación implicados:")
message("  Madre basal : ", paste(implicated$Madre_Basal, collapse = ", "))
message("  Bebé 6m     : ", paste(implicated$Bebe_mes_6,   collapse = ", "))

get_meth_dyad <- function(estrato) {
  genes <- implicated[[estrato]]
  if (!length(genes)) return(NULL)
  suf <- estr_suf[[estrato]]; tmp <- estr_met[[estrato]]
  m <- meth %>% filter(Tiempo == tmp) %>%
    select(Pareja, all_of(genes))
  names(m)[-1] <- paste0("met_", suf, "_", sanitize(sub("_mean$", "", genes)))
  m
}
met_mad <- get_meth_dyad("Madre_Basal")
met_b6m <- get_meth_dyad("Bebe_mes_6")

# NRF1 a 1 mes (b3): la metilación de NRF1 correlaciona con el WLZ_6 observado
# (B4_adiposidad, rho ≈ -0.51, q ≈ 0.15, sobrevive FDR). Ese tiempo (b3) no
# encaja en los estratos mad/b6m, así que se extrae aparte con sufijo b1m.
met_b1m <- meth %>%
  filter(Tiempo == "b3") %>%
  transmute(Pareja, met_b1m_NRF1 = NRF1_mean)

# (5) Ensamblar todo por Pareja y unir al wide por dyad_id 
feature_tabs <- Filter(Negate(is.null),
  list(tax_mad, tax_b6m, div_mad, div_bb, div_b6m,
       pwy_mad, vt_df, met_mad, met_b6m, met_b1m))
features <- reduce(feature_tabs, full_join, by = "Pareja")

base <- read_excel(WIDE_IN, sheet = "BASE DE DATOS")
base$Pareja <- dyad_to_pareja(base$dyad_id)

new_cols <- setdiff(names(features), "Pareja")
base_new <- base %>%
  left_join(features, by = "Pareja") %>%
  select(-Pareja)

# (6) Ampliar la hoja LEYENDA 
ley_raw <- read_excel(WIDE_IN, sheet = "LEYENDA", col_names = FALSE)
ley_cols <- c("Nombre en el Excel limpio",
              "Nombre original en el Excel crudo",
              "Descripción")
ley_body <- ley_raw
names(ley_body) <- ley_cols[seq_len(ncol(ley_body))]

# Diccionario de descripciones por prefijo/columna concreta
describe <- function(col) {
  if (grepl("^mb_mad_",  col)) return(paste0("NUEVA (microbiota). Abundancia CLR del género ",
                                             sub("^mb_mad_", "", col), " en MADRE basal (heces)."))
  if (grepl("^mb_b6m_",  col)) return(paste0("NUEVA (microbiota). Abundancia CLR del género ",
                                             sub("^mb_b6m_", "", col), " en BEBÉ 6 meses (heces)."))
  if (grepl("^div_mad_", col)) return(paste0("NUEVA (microbiota). Diversidad alfa (",
                                             sub("^div_mad_", "", col), ") en MADRE basal."))
  if (grepl("^div_bb_",  col)) return(paste0("NUEVA (microbiota). Diversidad alfa (",
                                             sub("^div_bb_", "", col), ") en BEBÉ basal."))
  if (grepl("^div_b6m_", col)) return(paste0("NUEVA (microbiota). Diversidad alfa (",
                                             sub("^div_b6m_", "", col), ") en BEBÉ 6 meses."))
  if (grepl("^pwy_mad_", col)) return(paste0("NUEVA (funcional PICRUSt2). Abundancia CLR de la vía MetaCyc ",
                                             sub("^pwy_mad_", "", col), " en MADRE basal (biosíntesis de ácidos grasos)."))
  if (grepl("^vt_.*_pct_shared", col)) return("NUEVA (transmisión vertical). % de ASVs del bebé compartidos con SU propia madre.")
  if (grepl("^vt_.*_bray_madre", col)) return("NUEVA (transmisión vertical). Distancia Bray-Curtis bebé↔su propia madre (menor = más similar).")
  if (grepl("^met_mad_", col)) return(paste0("NUEVA (epigenética). Metilación media (%) del gen ",
                                             sub("^met_mad_", "", col), " en MADRE basal (saliva, m2)."))
  if (grepl("^met_b6m_", col)) return(paste0("NUEVA (epigenética). Metilación media (%) del gen ",
                                             sub("^met_b6m_", "", col), " en BEBÉ 6 meses (saliva, b4)."))
  if (grepl("^met_b1m_", col)) return(paste0("NUEVA (epigenética). Metilación media (%) del gen ",
                                             sub("^met_b1m_", "", col), " en BEBÉ 1 mes (saliva, b3)."))
  "NUEVA. Variable multiómica."
}

ley_add <- tibble(
  `Nombre en el Excel limpio`         = c("VARIABLES MULTIÓMICAS — Microbiota 16S, funcional y epigenética (NUEVO)", new_cols),
  `Nombre original en el Excel crudo` = c(NA, rep("—", length(new_cols))),
  `Descripción`                       = c(NA, vapply(new_cols, describe, character(1)))
)
ley_new <- bind_rows(ley_body, ley_add)

# (7) Escribir Excel 
write_xlsx(list(`BASE DE DATOS` = base_new, `LEYENDA` = ley_new), WIDE_OUT)

# (8) Resumen
cat(sprintf("Columnas nuevas añadidas: %d\n", length(new_cols)))
cov <- sapply(new_cols, function(c) sum(!is.na(base_new[[c]])))
cat("\nCobertura (nº díadas con dato) por columna:\n")
print(data.frame(columna = new_cols, n_no_NA = as.integer(cov), row.names = NULL))
cat(sprintf("\nFichero escrito: %s\n", WIDE_OUT))
cat("\n14_exportar_features_a_wide.R completado.\n")
