#!/usr/bin/env bash
# Tests estadísticos sobre los resultados de diversidad del paso 06.
# Salidas: Visualizaciones .qzv con tests y boxplots.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_DIR}/qiime2_output"
METADATA="${PROJECT_DIR}/metadata/sample-metadata-basal-6m.tsv"
CORE_DIR="${OUTPUT_DIR}/06_core-metrics-d2000"
STATS_DIR="${OUTPUT_DIR}/07_stats"
mkdir -p "$STATS_DIR"

if ! command -v qiime &> /dev/null; then
    echo "ERROR: qiime no encontrado. Activa el entorno conda de QIIME2."
    exit 1
fi

[[ -d "$CORE_DIR" ]] || { echo "ERROR: No existe $CORE_DIR — ejecuta antes 06_diversidad.sh"; exit 1; }
[[ -f "$METADATA" ]] || { echo "ERROR: No existe $METADATA"; exit 1; }

# ── Bloque 1: alpha-group-significance ────────────────────────────────────────
echo "Tests de diversidad alpha (Kruskal-Wallis)"
ALPHA_METRICS=(shannon faith_pd observed_features evenness)

for metric in "${ALPHA_METRICS[@]}"; do
    OUT="${STATS_DIR}/alpha_${metric}_significance.qzv"
    if [[ ! -f "$OUT" ]]; then
        echo "  [α] $metric ..."
        qiime diversity alpha-group-significance \
            --i-alpha-diversity "${CORE_DIR}/${metric}_vector.qza" \
            --m-metadata-file "$METADATA" \
            --o-visualization "$OUT"
    else
        echo "  [α] $metric — ya existe, omitido"
    fi
done

# ── Bloque 2: PERMANOVA beta ──────────────────────────────────────────────────
echo ""
echo "Tests de diversidad beta (PERMANOVA, 999 permutaciones)"
BETA_METRICS=(bray_curtis jaccard weighted_unifrac unweighted_unifrac)
# Variables a testear (deben existir como columnas en el metadata)
VARIABLES=(Origen Tiempo Origen_tiempo Peso_Pre Sexo_bebe)

for metric in "${BETA_METRICS[@]}"; do
    DIST="${CORE_DIR}/${metric}_distance_matrix.qza"
    for var in "${VARIABLES[@]}"; do
        OUT="${STATS_DIR}/beta_${metric}_${var}.qzv"
        if [[ ! -f "$OUT" ]]; then
            echo "  [β] $metric × $var ..."
            qiime diversity beta-group-significance \
                --i-distance-matrix "$DIST" \
                --m-metadata-file "$METADATA" \
                --m-metadata-column "$var" \
                --p-method permanova \
                --p-permutations 999 \
                --p-pairwise \
                --o-visualization "$OUT" 2>/dev/null \
                || echo "    (saltado: posible columna ausente o n insuficiente)"
        else
            echo "  [β] $metric × $var — ya existe, omitido"
        fi
    done
done

echo "Salidas en ${STATS_DIR}/ (tests alpha Kruskal-Wallis y beta PERMANOVA)."
