#!/usr/bin/env Rscript
# Limpieza del fichero crudo de pirosequenciación de metilación (saliva).

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(tidyverse)
})

# Rutas
.file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
PROJECT_DIR <- if (length(.file_arg))
  normalizePath(file.path(dirname(sub("^--file=", "", .file_arg)), "..")) else normalizePath("..")
IN_FILE  <- file.path(PROJECT_DIR, "data",
                      "pirosecuenciacion_metilacion_crudo.xlsx")
OUT_FILE <- file.path(PROJECT_DIR, "data",
                      "pirosecuenciacion_metilacion_crudo_limpio.xlsx")

stopifnot(file.exists(IN_FILE))

# Posiciones que NO reflejan metilación (genotipo/SNP), fuera de las medias.
GENOTYPE_POS <- tibble(Gen = "BCAT1", CpG = "CpG4")

#  Si las dos réplicas de una misma CpG-muestra difieren MÁS que este valor, la medida se censura a NA antes de promediar. 
REPLICATE_RANGE_THRESHOLD <- 10

# Helpers 

# Convierte texto a numérico 
parse_meth <- function(x) {
  x <- as.character(x)
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- trimws(x)
  x[x %in% c("", "NA", "Failed", "FAILED", "failed", "-", "ND", "nd")] <- NA
  suppressWarnings(as.numeric(x))
}

# Normaliza formatos heredados de Sample ID antes de parsear.
normalize_sample_id <- function(id) {
  x <- trimws(as.character(id))
  # Sufijo de meses "6M" (con/ sin guion/espacio) -> código de tiempo b4
  x <- sub("(?i)[-_ ]*6\\s*M$", "-b4", x, perl = TRUE)
  x
}

# Parsea IDs
parse_sample <- function(id) {
  id <- as.character(id)
  m <- regmatches(id, regexec("^NEMO\\s*0*(\\d+)\\s*[-_]?\\s*([mb])\\s*(\\d+)$",
                              id, ignore.case = TRUE))
  Pareja <- character(length(id))
  Tiempo <- character(length(id))
  Origen <- character(length(id))
  for (i in seq_along(m)) {
    if (length(m[[i]]) == 4) {
      Pareja[i] <- sprintf("N%02d", as.integer(m[[i]][2]))
      letra     <- tolower(m[[i]][3])
      num       <- m[[i]][4]
      Tiempo[i] <- paste0(letra, num)            # m2, b2, b3, b4
      Origen[i] <- ifelse(letra == "m", "Madre", "Bebe")
    } else {
      Pareja[i] <- NA_character_
      Tiempo[i] <- NA_character_
      Origen[i] <- NA_character_
    }
  }
  tibble(Pareja = Pareja, Tiempo = Tiempo, Origen = Origen)
}

# Lee una hoja-gen del Excel y la devuelve en formato long limpio
read_gene_sheet <- function(file, sheet) {
  raw <- suppressMessages(read_excel(file, sheet = sheet, col_names = FALSE))
  headers <- as.character(unlist(raw[1, ]))
  pos_idx <- grep("^Pos\\.?\\s*\\d+$", headers)
  if (length(pos_idx) == 0) {
    warning(sprintf("Hoja %s: no se detectaron columnas Pos. X", sheet))
    return(NULL)
  }
  pos_labels <- paste0("CpG", seq_along(pos_idx))

  data <- raw[-c(1, 2), c(1, pos_idx), drop = FALSE]
  colnames(data) <- c("Sample_ID_raw", pos_labels)

  # Quitar filas vacías y controles de calibración 
  data <- data |>
    filter(!is.na(Sample_ID_raw) & Sample_ID_raw != "" &
             Sample_ID_raw != "Sample ID")

  # normalizar IDs heredados, parsear y AVISAR de los no parseables
  data <- data |> mutate(Sample_ID_norm = normalize_sample_id(Sample_ID_raw))
  meta <- parse_sample(data$Sample_ID_norm)

  no_parse <- unique(data$Sample_ID_raw[is.na(meta$Pareja)])
  no_parse <- no_parse[!grepl("^\\*|\\s.+\\s", no_parse)]
  if (length(no_parse) > 0) {
    warning(sprintf("Hoja %s: %d Sample ID NO parseables (muestras perdidas): %s",
                    sheet, length(no_parse), paste(no_parse, collapse = ", ")))
  }

  cpg_mat <- data |>
    select(all_of(pos_labels)) |>
    mutate(across(everything(), parse_meth))

  bind_cols(
    tibble(Sample_ID_raw = data$Sample_ID_raw),
    meta,
    cpg_mat
  ) |>
    mutate(Gen = sheet, .before = 1)
}

# (1) Leer todas las hojas 

sheets <- excel_sheets(IN_FILE)
message("Hojas detectadas: ", paste(sheets, collapse = ", "))

raw_list <- map(sheets, ~ read_gene_sheet(IN_FILE, .x))
names(raw_list) <- sheets

# (2) Construir long format con réplicas etiquetadas 

long_replicates <- raw_list |>
  map(~ .x |>
        pivot_longer(starts_with("CpG"), names_to = "CpG", values_to = "valor")) |>
  bind_rows() |>
  filter(!is.na(Pareja), !is.na(Tiempo)) |>
  group_by(Gen, Pareja, Tiempo, CpG) |>
  mutate(replica = row_number()) |>
  ungroup()

# (3) QC: diferencia entre réplicas 

qc_replicas <- long_replicates |>
  group_by(Gen, Pareja, Tiempo, CpG) |>
  summarise(
    n_replicas  = sum(!is.na(valor)),
    media       = mean(valor, na.rm = TRUE),   # media SIN censurar
    rango       = if (sum(!is.na(valor)) >= 2) diff(range(valor, na.rm = TRUE))
                  else NA_real_,
    .groups = "drop"
  ) |>
  mutate(
    media     = ifelse(is.nan(media), NA_real_, media),
    censurado = !is.na(rango) & rango > REPLICATE_RANGE_THRESHOLD,
    flag      = ifelse(censurado,
                       sprintf("Réplicas difieren > %g %% -> CENSURADO a NA",
                               REPLICATE_RANGE_THRESHOLD),
                       NA_character_)
  )

n_flag <- sum(qc_replicas$censurado)
message(sprintf("QC: %d CpG-muestra con réplicas que difieren > %g%% (censuradas a NA)",
                n_flag, REPLICATE_RANGE_THRESHOLD))

# (4) Long format definitivo (réplicas promediadas) 

long_summ <- long_replicates |>
  group_by(Gen, Sample_ID_raw, Pareja, Tiempo, Origen, CpG) |>
  summarise(
    rango = if (sum(!is.na(valor)) >= 2) diff(range(valor, na.rm = TRUE))
            else NA_real_,
    valor = mean(valor, na.rm = TRUE),
    .groups = "drop"
  ) |>
  # Si después de promediar queda NaN (todo NA), pasar a NA
  mutate(valor = ifelse(is.nan(valor), NA_real_, valor),
         #  > umbral es NA
         censurado = !is.na(rango) & rango > REPLICATE_RANGE_THRESHOLD)

n_cens <- sum(long_summ$censurado)
message(sprintf("CpG-muestra censuradas a NA (rango > %g %%): %d",
                REPLICATE_RANGE_THRESHOLD, n_cens))


long_clean <- long_summ |>
  mutate(valor = ifelse(censurado, NA_real_, valor)) |>
  select(Gen, Sample_ID_raw, Pareja, Tiempo, Origen, CpG, valor)

# (5) Wide format: 1 fila por muestra biológica 

# Etiqueta de columna
wide <- long_clean |>
  mutate(col = paste(Gen, CpG, sep = "_")) |>
  select(Pareja, Tiempo, Origen, col, valor) |>
  pivot_wider(names_from = col, values_from = valor) |>
  arrange(Pareja, Tiempo)

# medias por gen EXCLUYENDO posiciones de genotipo 
gene_means <- long_clean |>
  anti_join(GENOTYPE_POS, by = c("Gen", "CpG")) |>
  group_by(Pareja, Tiempo, Origen, Gen) |>
  summarise(media = mean(valor, na.rm = TRUE), .groups = "drop") |>
  mutate(media = ifelse(is.nan(media), NA_real_, media),
         col = paste0(Gen, "_mean")) |>
  select(Pareja, Tiempo, Origen, col, media) |>
  pivot_wider(names_from = col, values_from = media)

wide <- wide |>
  left_join(gene_means, by = c("Pareja", "Tiempo", "Origen")) |>
  mutate(sample_id_meth = paste0("NEMO", sub("^N", "", Pareja), "-", Tiempo),
         .after = Origen)

# marcar la columna de genotipo para que no se confunda con metilación
if ("BCAT1_CpG4" %in% names(wide)) {
  wide <- wide |> rename(BCAT1_CpG4_genotype = BCAT1_CpG4)
}

# Ordenar columnas
gene_order <- sheets
ordered_cols <- c("Pareja", "Tiempo", "Origen", "sample_id_meth")
for (g in gene_order) {
  cols_g <- c(paste0(g, "_mean"),
              grep(paste0("^", g, "_CpG"), names(wide), value = TRUE))
  ordered_cols <- c(ordered_cols, intersect(cols_g, names(wide)))
}
wide <- wide |> select(all_of(ordered_cols))

# (6) Hojas por gen (limpias) 

per_gene_sheets <- map(sheets, function(g) {
  out <- long_clean |>
    filter(Gen == g) |>
    pivot_wider(names_from = CpG, values_from = valor) |>
    arrange(Pareja, Tiempo) |>
    select(-Gen)
  # reflejar también aquí el renombrado del genotipo
  if (g == "BCAT1" && "CpG4" %in% names(out)) {
    out <- out |> rename(CpG4_genotype = CpG4)
  }
  out
})
names(per_gene_sheets) <- sheets

# (7) Guardar Excel multi-hoja 

out_sheets <- c(
  list(
    wide         = wide,
    long         = long_clean,
    qc_replicas  = qc_replicas
  ),
  per_gene_sheets
)

write_xlsx(out_sheets, OUT_FILE)

# (8) Resumen 

cat(sprintf("\nLimpieza completada.\n"))
cat(sprintf("  Muestras biológicas únicas: %d\n", nrow(wide)))
cat(sprintf("  Genes: %d (%s)\n", length(sheets), paste(sheets, collapse = ", ")))
cat(sprintf("  Total CpGs de metilación: %d\n",
            sum(grepl("_CpG\\d+$", names(wide)))))
cat("  BCAT1_mean calculado SIN CpG4 (genotipo); columna -> BCAT1_CpG4_genotype\n")
cat(sprintf("  %d CpG-muestra censuradas a NA por réplicas discordantes (> %g %%)\n",
            n_cens, REPLICATE_RANGE_THRESHOLD))
cat("\nDistribución por tiempo y origen:\n")
print(wide |> count(Origen, Tiempo))
cat("\nComprobación: ¿recuperada N04 a 6 meses (b4)?\n")
print(wide |> filter(Pareja == "N04") |> select(Pareja, Tiempo))
cat(sprintf("\nFichero generado:\n  %s\n", OUT_FILE))
