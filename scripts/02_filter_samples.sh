#!/usr/bin/env bash
# Filtra la feature table completa (todos los tiempos) a basal y mes_6
# Entrada:  01_table_all.qza, 01_rep-seqs_all.qza
# Salida:   01_table.qza, 01_rep-seqs.qza  
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_DIR}/qiime2_output"
METADATA="${PROJECT_DIR}/metadata/sample-metadata-basal-6m.tsv"

if ! command -v qiime &> /dev/null; then
    echo "ERROR: qiime no encontrado. Activa el entorno conda de QIIME2."
    exit 1
fi

[[ -f "${OUTPUT_DIR}/01_table_all.qza" ]]    || { echo "ERROR: No existe 01_table_all.qza — ejecuta primero 01_dada2.sh"; exit 1; }
[[ -f "${OUTPUT_DIR}/01_rep-seqs_all.qza" ]] || { echo "ERROR: No existe 01_rep-seqs_all.qza"; exit 1; }
[[ -f "$METADATA" ]]                         || { echo "ERROR: No existe $METADATA"; exit 1; }

echo "Filtrando feature table a muestras basal + mes_6..."

qiime feature-table filter-samples \
    --i-table "${OUTPUT_DIR}/01_table_all.qza" \
    --m-metadata-file "$METADATA" \
    --o-filtered-table "${OUTPUT_DIR}/01_table.qza"

echo "Filtrando rep-seqs a features presentes en la tabla filtrada..."

qiime feature-table filter-seqs \
    --i-data "${OUTPUT_DIR}/01_rep-seqs_all.qza" \
    --i-table "${OUTPUT_DIR}/01_table.qza" \
    --o-filtered-data "${OUTPUT_DIR}/01_rep-seqs.qza"

echo "Generando visualizaciones..."

qiime feature-table summarize \
    --i-table "${OUTPUT_DIR}/01_table.qza" \
    --o-visualization "${OUTPUT_DIR}/01_table.qzv"

qiime feature-table tabulate-seqs \
    --i-data "${OUTPUT_DIR}/01_rep-seqs.qza" \
    --o-visualization "${OUTPUT_DIR}/01_rep-seqs.qzv"

echo "Salidas en ${OUTPUT_DIR}: 01_table, 01_rep-seqs (basal + mes_6)."
