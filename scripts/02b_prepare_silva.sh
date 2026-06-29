#!/usr/bin/env bash
# Descarga las referencias SILVA 138 y extrae la región V3-V4 
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SILVA_DIR="${PROJECT_DIR}/qiime2_output/silva138"
mkdir -p "$SILVA_DIR"

# URLs oficiales de QIIME2 (release 2024.10)
URL_SEQS="https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza"
URL_TAX="https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza"

REF_SEQS="${SILVA_DIR}/silva-138-99-seqs.qza"
REF_TAX="${SILVA_DIR}/silva-138-99-tax.qza"
EXT_SEQS="${SILVA_DIR}/silva-138-99-seqs-V3V4.qza"

if ! command -v qiime &> /dev/null; then
    echo "ERROR: qiime no encontrado. Activa el entorno conda de QIIME2."
    exit 1
fi

# Paso 1: descargar referencias SILVA 138 
if [[ ! -f "$REF_SEQS" ]]; then
    echo "[1/2] Descargando secuencias de referencia SILVA 138..."
    curl -L -o "$REF_SEQS" "$URL_SEQS"
else
    echo "[1/2] SILVA 138 seqs ya presente, omitido."
fi

if [[ ! -f "$REF_TAX" ]]; then
    echo "      Descargando taxonomía de referencia SILVA 138..."
    curl -L -o "$REF_TAX" "$URL_TAX"
else
    echo "      SILVA 138 tax ya presente, omitido."
fi

# Paso 2: extraer región V3-V4 con los primers reales 
if [[ ! -f "$EXT_SEQS" ]]; then
    echo "[2/2] Extrayendo región V3-V4 (primers 341F/805R) de la referencia..."
    qiime feature-classifier extract-reads \
        --i-sequences "$REF_SEQS" \
        --p-f-primer CCTACGGGNGGCWGCAG \
        --p-r-primer GACTACHVGGTATCTAATCC \
        --p-min-length 250 \
        --p-max-length 500 \
        --p-n-jobs 1 \
        --o-reads "$EXT_SEQS" \
        --verbose
else
    echo "[2/2] Región V3-V4 ya extraída, omitido."
fi

echo ""
echo "Referencias SILVA preparadas:"
echo "  $EXT_SEQS  — secuencias V3-V4"
echo "  $REF_TAX   — taxonomía"
echo ""
echo "Siguiente paso: bash scripts/03_taxonomia.sh"