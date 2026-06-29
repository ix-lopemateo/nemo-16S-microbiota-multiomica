#!/usr/bin/env bash
# Predicción funcional con PICRUSt2 a partir de las ASVs filtradas

set -euo pipefail

# Rutas 
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
QIIME_DIR="${PROJECT_DIR}/qiime2_output"
PICRUST_DIR="${QIIME_DIR}/picrust2_output"
EXPORT_DIR="${QIIME_DIR}/picrust2_input"

REP_SEQS_QZA="${QIIME_DIR}/04_rep-seqs.qza"
TABLE_QZA="${QIIME_DIR}/04_table.qza"

# Threads
N_THREADS=4

# Imágenes Docker
QIIME_IMG="quay.io/qiime2/amplicon:2024.10"
PICRUST_IMG="quay.io/biocontainers/picrust2:2.5.3--pyhdfd78af_0"

# Validaciones 
command -v docker >/dev/null || { echo "ERROR: docker no encontrado." >&2; exit 1; }
[[ -f "$REP_SEQS_QZA" ]] || { echo "ERROR: $REP_SEQS_QZA no existe." >&2; exit 1; }
[[ -f "$TABLE_QZA" ]]    || { echo "ERROR: $TABLE_QZA no existe." >&2; exit 1; }

mkdir -p "$EXPORT_DIR"

# Paso 1: exportar rep-seqs y tabla desde QIIME2 a FASTA y BIOM 
if [[ ! -f "${EXPORT_DIR}/dna-sequences.fasta" ]]; then
  echo "[1/4] Exportando rep-seqs.qza → FASTA..."
  docker run --rm --platform linux/amd64 \
    -v "${PROJECT_DIR}:/data" \
    "$QIIME_IMG" \
    qiime tools export \
      --input-path "/data/qiime2_output/04_rep-seqs.qza" \
      --output-path "/data/qiime2_output/picrust2_input"
else
  echo "[1/4] FASTA ya existe — omitido."
fi

if [[ ! -f "${EXPORT_DIR}/feature-table.biom" ]]; then
  echo "[2/4] Exportando table.qza → BIOM..."
  docker run --rm --platform linux/amd64 \
    -v "${PROJECT_DIR}:/data" \
    "$QIIME_IMG" \
    qiime tools export \
      --input-path "/data/qiime2_output/04_table.qza" \
      --output-path "/data/qiime2_output/picrust2_input"
else
  echo "[2/4] BIOM ya existe — omitido."
fi

# Paso 2: PICRUSt2 full pipeline 
# Borramos picrust2_output si existe para evitar error 
if [[ -d "$PICRUST_DIR" ]]; then
  echo "[3/4] Directorio picrust2_output existe — eliminando para re-ejecutar..."
  rm -rf "$PICRUST_DIR"
fi

echo "[3/4] Ejecutando picrust2_pipeline.py (esto puede tardar 20-45 min con SEPP)..."
docker run --rm --platform linux/amd64 \
  -v "${PROJECT_DIR}:/data" \
  "$PICRUST_IMG" \
  picrust2_pipeline.py \
    -s "/data/qiime2_output/picrust2_input/dna-sequences.fasta" \
    -i "/data/qiime2_output/picrust2_input/feature-table.biom" \
    -o "/data/qiime2_output/picrust2_output" \
    -p $N_THREADS \
    --placement_tool sepp \
    --hsp_method emp_prob \
    --verbose

# Paso 3: añadir descripciones a KO y pathways 
echo "[4/4] Añadiendo descripciones a KO y pathways..."
docker run --rm --platform linux/amd64 \
  -v "${PROJECT_DIR}:/data" \
  "$PICRUST_IMG" bash -c "
    add_descriptions.py \
      -i /data/qiime2_output/picrust2_output/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \
      -m KO \
      -o /data/qiime2_output/picrust2_output/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz &&
    add_descriptions.py \
      -i /data/qiime2_output/picrust2_output/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
      -m EC \
      -o /data/qiime2_output/picrust2_output/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz &&
    add_descriptions.py \
      -i /data/qiime2_output/picrust2_output/pathways_out/path_abun_unstrat.tsv.gz \
      -m METACYC \
      -o /data/qiime2_output/picrust2_output/pathways_out/path_abun_unstrat_descrip.tsv.gz
  "

echo "PICRUSt2 completado. Salidas en ${PICRUST_DIR}/ (KO, EC, pathways con descripciones; marker_predicted_and_nsti para QC)."
