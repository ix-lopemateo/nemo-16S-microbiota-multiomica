#!/usr/bin/env bash
# Cálculo de diversidad alpha y beta a partir de la tabla filtrada (04) y del árbol filogenético enraizado (05).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_DIR}/qiime2_output"
METADATA="${PROJECT_DIR}/metadata/sample-metadata-basal-6m.tsv"

TABLE="${OUTPUT_DIR}/04_table.qza"
TREE="${OUTPUT_DIR}/05_rooted-tree.qza"

SAMPLING_DEPTH="${SAMPLING_DEPTH:-2000}"
MAX_DEPTH_RAREF="${MAX_DEPTH_RAREF:-30000}"
CORE_DIR="${OUTPUT_DIR}/06_core-metrics-d${SAMPLING_DEPTH}"

if ! command -v qiime &> /dev/null; then
    echo "ERROR: qiime no encontrado. Activa el entorno conda de QIIME2."
    exit 1
fi

[[ -f "$TABLE" ]]    || { echo "ERROR: No existe $TABLE — ejecuta antes 04_filter_taxonomy.sh"; exit 1; }
[[ -f "$TREE" ]]     || { echo "ERROR: No existe $TREE — ejecuta antes 05_arbol_filogenetico.sh"; exit 1; }
[[ -f "$METADATA" ]] || { echo "ERROR: No existe $METADATA"; exit 1; }

# ── Paso 1: curvas de rarefacción (diagnóstico) ───────────────────────────────
if [[ ! -f "${OUTPUT_DIR}/06_alpha-rarefaction.qzv" ]]; then
    echo "[1/2] Generando curvas de rarefacción hasta profundidad ${MAX_DEPTH_RAREF}..."
    qiime diversity alpha-rarefaction \
        --i-table "$TABLE" \
        --i-phylogeny "$TREE" \
        --p-max-depth "$MAX_DEPTH_RAREF" \
        --p-min-depth 100 \
        --p-steps 20 \
        --p-iterations 10 \
        --m-metadata-file "$METADATA" \
        --o-visualization "${OUTPUT_DIR}/06_alpha-rarefaction.qzv"
else
    echo "[1/2] Curvas de rarefacción ya generadas, omitido."
fi

# ── Paso 2: core-metrics-phylogenetic ─────────────────────────────────────────
if [[ ! -d "$CORE_DIR" ]]; then
    echo "[2/2] Ejecutando core-metrics-phylogenetic a profundidad ${SAMPLING_DEPTH}..."
    qiime diversity core-metrics-phylogenetic \
        --i-table "$TABLE" \
        --i-phylogeny "$TREE" \
        --p-sampling-depth "$SAMPLING_DEPTH" \
        --m-metadata-file "$METADATA" \
        --p-n-jobs-or-threads 1 \
        --output-dir "$CORE_DIR" \
        --verbose
else
    echo "[2/2] core-metrics ya existe en ${CORE_DIR}, omitido."
    echo "      (borra el directorio si quieres recalcular)"
fi

echo "Salidas: 06_alpha-rarefaction.qzv y ${CORE_DIR}/ (métricas alpha y beta, PCoA y Emperor plots)."
