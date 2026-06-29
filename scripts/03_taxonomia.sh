#!/usr/bin/env bash
# Asignación taxonómica de las ASVs por consenso con vsearch contra la referencia SILVA 138 V3-V4 
# Entrada: 01_rep-seqs.qza, silva-138-99-seqs-V3V4.qza, silva-138-99-tax.qza
# Salida:  03_taxonomy.qza, 03_taxonomy.qzv, 03_taxa_barplot.qzv
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_DIR}/qiime2_output"

REP_SEQS="${OUTPUT_DIR}/01_rep-seqs.qza"
TABLE="${OUTPUT_DIR}/01_table.qza"
REF_SEQS="${OUTPUT_DIR}/silva138/silva-138-99-seqs-V3V4.qza"
REF_TAX="${OUTPUT_DIR}/silva138/silva-138-99-tax.qza"
METADATA="${PROJECT_DIR}/metadata/sample-metadata-basal-6m.tsv"

if ! command -v qiime &> /dev/null; then
    echo "ERROR: qiime no encontrado. Activa el entorno conda de QIIME2."
    exit 1
fi

[[ -f "$REP_SEQS" ]] || { echo "ERROR: No existe $REP_SEQS — ejecuta antes 02_filter_samples.sh"; exit 1; }
[[ -f "$TABLE" ]]    || { echo "ERROR: No existe $TABLE"; exit 1; }
[[ -f "$REF_SEQS" ]] || { echo "ERROR: No existe $REF_SEQS — ejecuta antes 02b_prepare_silva.sh"; exit 1; }
[[ -f "$REF_TAX" ]]  || { echo "ERROR: No existe $REF_TAX"; exit 1; }
[[ -f "$METADATA" ]] || { echo "ERROR: No existe $METADATA"; exit 1; }

# Paso 1: clasificación por consenso vsearch 
echo "[1/3] Clasificando ASVs por consenso vsearch contra SILVA 138 V3-V4..."
qiime feature-classifier classify-consensus-vsearch \
    --i-query "$REP_SEQS" \
    --i-reference-reads "$REF_SEQS" \
    --i-reference-taxonomy "$REF_TAX" \
    --p-perc-identity 0.80 \
    --p-maxaccepts 10 \
    --p-min-consensus 0.51 \
    --p-threads 1 \
    --o-classification "${OUTPUT_DIR}/03_taxonomy.qza" \
    --o-search-results "${OUTPUT_DIR}/03_vsearch_hits.qza" \
    --verbose

# Paso 2: tabular taxonomía 
echo "[2/3] Generando visualización tabular de la taxonomía..."
qiime metadata tabulate \
    --m-input-file "${OUTPUT_DIR}/03_taxonomy.qza" \
    --o-visualization "${OUTPUT_DIR}/03_taxonomy.qzv"

# Paso 3: barplot por muestra 
echo "[3/3] Generando barplot taxonómico interactivo..."
qiime taxa barplot \
    --i-table "$TABLE" \
    --i-taxonomy "${OUTPUT_DIR}/03_taxonomy.qza" \
    --m-metadata-file "$METADATA" \
    --o-visualization "${OUTPUT_DIR}/03_taxa_barplot.qzv"

echo "Salidas en ${OUTPUT_DIR}: 03_taxonomy, 03_taxa_barplot, 03_vsearch_hits. Revisar el barplot (mito/cloroplastos, ASVs sin Phylum) antes del filtrado."