#!/usr/bin/env bash
# Construye el árbol filogenético de las ASVs limpias mediante MAFFT (alineamiento) → máscara de
# posiciones variables → FastTree (árbol de máxima verosimilitud aproximada) → enraizado por midpoint.
# Entrada: 04_rep-seqs.qza
# Salida:  05_aligned-rep-seqs.qza, 05_masked-aligned-rep-seqs.qza, 05_unrooted-tree.qza,05_rooted-tree.qza  
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_DIR}/qiime2_output"
REP_SEQS="${OUTPUT_DIR}/04_rep-seqs.qza"

if ! command -v qiime &> /dev/null; then
    echo "ERROR: qiime no encontrado. Activa el entorno conda de QIIME2."
    exit 1
fi

[[ -f "$REP_SEQS" ]] || { echo "ERROR: No existe $REP_SEQS — ejecuta antes 04_filter_taxonomy.sh"; exit 1; }

echo "Construyendo árbol filogenético (MAFFT → mask → FastTree → midpoint root)..."
echo "Esto puede tardar 20-60 minutos para ~3000 ASVs."

qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "$REP_SEQS" \
    --p-n-threads 1 \
    --o-alignment              "${OUTPUT_DIR}/05_aligned-rep-seqs.qza" \
    --o-masked-alignment       "${OUTPUT_DIR}/05_masked-aligned-rep-seqs.qza" \
    --o-tree                   "${OUTPUT_DIR}/05_unrooted-tree.qza" \
    --o-rooted-tree            "${OUTPUT_DIR}/05_rooted-tree.qza" \
    --verbose

echo "Salidas en ${OUTPUT_DIR}: 05_aligned-rep-seqs, 05_masked-aligned-rep-seqs, 05_unrooted-tree, 05_rooted-tree."
