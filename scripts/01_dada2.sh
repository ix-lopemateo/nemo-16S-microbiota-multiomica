#!/usr/bin/env bash
# Denoising paired-end con DADA2
# Proyecto NEMO V3-V4, MiSeq 300×2 sobre TODOS los tiempos
set -euo pipefail

# Rutas y variables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
INPUT_QZA="${PROJECT_DIR}/qiime2_output/00b_reads_trimmed.qza"
OUTPUT_DIR="${PROJECT_DIR}/qiime2_output"

# Validaciones previas 
if ! command -v qiime &> /dev/null; then
    echo "ERROR: qiime no encontrado. Activa el entorno conda de QIIME2."
    exit 1
fi

[[ -f "$INPUT_QZA" ]] || { echo "ERROR: No existe $INPUT_QZA — ejecuta primero 00_importar_reads.sh"; exit 1; }

mkdir -p "$OUTPUT_DIR"

# Longitudes de truncado fijadas a partir del perfil de calidad (resumen de
# 00b_cutadapt.sh). trim-left = 0 porque los primers ya se quitaron con cutadapt.
# n-threads = 1 para reproducibilidad.
echo "Ejecutando DADA2 paired-end denoising..."

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$INPUT_QZA" \
    --p-trunc-len-f 260 \
    --p-trunc-len-r 200 \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-max-ee-f 2 \
    --p-max-ee-r 2 \
    --p-chimera-method consensus \
    --p-n-reads-learn 100000 \
    --p-n-threads 1 \
    --o-table              "${OUTPUT_DIR}/01_table_all.qza" \
    --o-representative-sequences "${OUTPUT_DIR}/01_rep-seqs_all.qza" \
    --o-denoising-stats    "${OUTPUT_DIR}/01_denoising-stats_all.qza" \
    --verbose

echo "DADA2 completado."

# Visualizaciones 
echo "Generando visualizaciones..."

qiime feature-table summarize \
    --i-table "${OUTPUT_DIR}/01_table_all.qza" \
    --o-visualization "${OUTPUT_DIR}/01_table_all.qzv"

qiime feature-table tabulate-seqs \
    --i-data "${OUTPUT_DIR}/01_rep-seqs_all.qza" \
    --o-visualization "${OUTPUT_DIR}/01_rep-seqs_all.qzv"

qiime metadata tabulate \
    --m-input-file "${OUTPUT_DIR}/01_denoising-stats_all.qza" \
    --o-visualization "${OUTPUT_DIR}/01_denoising-stats_all.qzv"

echo "Salidas en ${OUTPUT_DIR}: 01_table_all, 01_rep-seqs_all, 01_denoising-stats_all (.qza/.qzv)."