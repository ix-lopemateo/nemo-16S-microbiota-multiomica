#!/usr/bin/env bash
# Importación de reads paired-end a QIIME2 (Casava 1.8 demultiplexed)
# Proyecto NEMO a TODOS los tiempos (B2,B3,B4,B5,M1,M2)
set -euo pipefail

# Rutas
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
R1_DIR="${PROJECT_DIR}/Reads/R1"
R2_DIR="${PROJECT_DIR}/Reads/R2"
METADATA="${PROJECT_DIR}/metadata/sample-metadata-full.tsv"
CASAVA_DIR="${PROJECT_DIR}/Reads/casava_input"
OUTPUT_DIR="${PROJECT_DIR}/qiime2_output"
OUTPUT_QZA="${OUTPUT_DIR}/00_reads_all.qza"

# Validaciones previas 
if ! command -v qiime &> /dev/null; then
    echo "ERROR: qiime no encontrado. Activa el entorno conda de QIIME2."
    exit 1
fi

[[ -d "$R1_DIR" ]] || { echo "ERROR: No existe $R1_DIR"; exit 1; }
[[ -d "$R2_DIR" ]] || { echo "ERROR: No existe $R2_DIR"; exit 1; }
[[ -f "$METADATA" ]] || { echo "ERROR: No existe $METADATA"; exit 1; }

mkdir -p "$CASAVA_DIR" "$OUTPUT_DIR"

# Crear symlinks con nomenclatura Casava 
echo "Creando directorio Casava con symlinks..."

# Limpiar symlinks previos para que runs no se mezclen si se ejecuta varias veces 
find "$CASAVA_DIR" -type l -delete

sample_counter=1
missing_samples=()

# Leer SampleIDs del metadata (saltando la cabecera #SampleID)
while IFS=$'\t' read -r sample_id _rest; do
    [[ "$sample_id" == "#SampleID" ]] && continue
    [[ -z "$sample_id" ]] && continue

    # Buscar archivo R1 cuyo nombre empiece por el SampleID seguido de _
    r1_file=$(find "$R1_DIR" -maxdepth 1 -name "${sample_id}_*_R1.fastq.gz" | head -1)
    r2_file=$(find "$R2_DIR" -maxdepth 1 -name "${sample_id}_*_R2.fastq.gz" | head -1)

    if [[ -z "$r1_file" || -z "$r2_file" ]]; then
        missing_samples+=("$sample_id")
        continue
    fi

    s_tag="S${sample_counter}"
    ln -sf "$r1_file" "${CASAVA_DIR}/${sample_id}_${s_tag}_L001_R1_001.fastq.gz"
    ln -sf "$r2_file" "${CASAVA_DIR}/${sample_id}_${s_tag}_L001_R2_001.fastq.gz"

    (( sample_counter++ ))
done < "$METADATA"

# Sacar informe de muestras no encontradas
if [[ ${#missing_samples[@]} -gt 0 ]]; then
    echo "AVISO: No se encontraron archivos FASTQ para las siguientes muestras:"
    printf '  %s\n' "${missing_samples[@]}"
fi

echo "Symlinks creados: $((sample_counter - 1)) muestras"

# Importar a QIIME2 
echo "Importando reads a QIIME2..."

qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$CASAVA_DIR" \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path "$OUTPUT_QZA"

echo "Importación completada: ${OUTPUT_QZA}"

#  Sacar resumen de reads importados 
echo "Generando resumen de reads..."
qiime demux summarize \
    --i-data "$OUTPUT_QZA" \
    --o-visualization "${OUTPUT_DIR}/00_reads_all_summary.qzv"

echo "Resumen: ${OUTPUT_DIR}/00_reads_all_summary.qzv"
