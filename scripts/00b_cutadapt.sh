#!/usr/bin/env bash
# Recorte de primers V3-V4 con cutadapt
# Primers V3-V4 
# Reads sin primer se descartan.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_DIR}/qiime2_output"
INPUT_QZA="${OUTPUT_DIR}/00_reads_all.qza"
OUTPUT_QZA="${OUTPUT_DIR}/00b_reads_trimmed.qza"

# Validaciones previas
if ! command -v qiime &> /dev/null; then
    echo "ERROR: qiime no encontrado. Activa el entorno conda de QIIME2."
    exit 1
fi

[[ -f "$INPUT_QZA" ]] || { echo "ERROR: No existe $INPUT_QZA — ejecuta primero 00_importar_reads.sh"; exit 1; }

echo "Recortando primers V3-V4 con cutadapt..."

# Primers V3-V4 
qiime cutadapt trim-paired \
    --i-demultiplexed-sequences "$INPUT_QZA" \
    --p-front-f CCTACGGGNGGCWGCAG \
    --p-front-r GACTACHVGGTATCTAATCC \
    --p-discard-untrimmed \
    --p-cores 1 \
    --o-trimmed-sequences "$OUTPUT_QZA" \
    --verbose

#Sacar resumen de reads recortados
echo "Generando resumen de reads recortados..."
qiime demux summarize \
    --i-data "$OUTPUT_QZA" \
    --o-visualization "${OUTPUT_DIR}/00b_reads_trimmed_summary.qzv"

echo "Salidas: ${OUTPUT_QZA} y 00b_reads_trimmed_summary.qzv. Revisar el perfil de calidad antes de fijar --p-trunc-len en 01_dada2.sh."