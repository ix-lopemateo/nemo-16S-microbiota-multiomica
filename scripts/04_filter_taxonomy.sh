#!/usr/bin/env bash
# Filtra ASVs problemáticas tras la asignación taxonómica
# Entrada: 01_table.qza, 01_rep-seqs.qza, 03_taxonomy.qza
# Salida:  04_table.qza, 04_rep-seqs.qza  (versiones limpias para análisis downstream)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_DIR}/qiime2_output"
METADATA="${PROJECT_DIR}/metadata/sample-metadata-basal-6m.tsv"

TABLE="${OUTPUT_DIR}/01_table.qza"
REP_SEQS="${OUTPUT_DIR}/01_rep-seqs.qza"
TAXONOMY="${OUTPUT_DIR}/03_taxonomy.qza"

if ! command -v qiime &> /dev/null; then
    echo "ERROR: qiime no encontrado. Activa el entorno conda de QIIME2."
    exit 1
fi

[[ -f "$TABLE" ]]    || { echo "ERROR: No existe $TABLE"; exit 1; }
[[ -f "$REP_SEQS" ]] || { echo "ERROR: No existe $REP_SEQS"; exit 1; }
[[ -f "$TAXONOMY" ]] || { echo "ERROR: No existe $TAXONOMY — ejecuta antes 03_taxonomia.sh"; exit 1; }
[[ -f "$METADATA" ]] || { echo "ERROR: No existe $METADATA"; exit 1; }

# Paso 1: filtrar tabla por taxonomía 
echo "[1/4] Filtrando feature table (mito/cloro/sin Phylum)..."
qiime taxa filter-table \
    --i-table "$TABLE" \
    --i-taxonomy "$TAXONOMY" \
    --p-include "p__" \
    --p-exclude "mitochondria,chloroplast" \
    --o-filtered-table "${OUTPUT_DIR}/04_table.qza"

# Paso 2: filtrar rep-seqs 
echo "[2/4] Filtrando rep-seqs..."
qiime taxa filter-seqs \
    --i-sequences "$REP_SEQS" \
    --i-taxonomy "$TAXONOMY" \
    --p-include "p__" \
    --p-exclude "mitochondria,chloroplast" \
    --o-filtered-sequences "${OUTPUT_DIR}/04_rep-seqs.qza"

# Paso 3: visualizaciones 
echo "[3/4] Generando resúmenes..."
qiime feature-table summarize \
    --i-table "${OUTPUT_DIR}/04_table.qza" \
    --o-visualization "${OUTPUT_DIR}/04_table.qzv"

qiime feature-table tabulate-seqs \
    --i-data "${OUTPUT_DIR}/04_rep-seqs.qza" \
    --o-visualization "${OUTPUT_DIR}/04_rep-seqs.qzv"

# Barplot taxonómico tras el filtrado
echo "[4/4] Generando barplot taxonómico limpio..."
qiime taxa barplot \
    --i-table "${OUTPUT_DIR}/04_table.qza" \
    --i-taxonomy "$TAXONOMY" \
    --m-metadata-file "$METADATA" \
    --o-visualization "${OUTPUT_DIR}/04_taxa_barplot.qzv"

echo "Salidas (limpias) en ${OUTPUT_DIR}: 04_table, 04_rep-seqs, 04_taxa_barplot. Los pasos siguientes usan los 04_*."