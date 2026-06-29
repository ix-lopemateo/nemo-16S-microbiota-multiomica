# NEMO — Microbiota 16S rRNA e integración multiómica

Código de análisis de la cohorte prospectiva **NEMO** (díadas madre–lactante), que
estudia la microbiota intestinal (16S rRNA, región V3–V4), la metilación de genes
candidatos y su relación con la adiposidad infantil temprana.

Este repositorio acompaña al Trabajo de Fin de Máster como **material para la
reproducibilidad computacional**. Contiene el pipeline completo, desde las lecturas
crudas hasta los modelos de regresión, organizado en scripts numerados.

> **Datos no incluidos.** Por tratarse de una cohorte humana sujeta a protección de
> datos, el repositorio **no contiene** las lecturas crudas (FASTQ) ni los datos
> clínicos, antropométricos o de metilación. Solo se publica el código. Los datos
> quedan disponibles bajo petición razonada y acceso restringido (ver
> [`data/README.md`](data/README.md)).

## Estructura

```
nemo-16S-microbiota-multiomica/
├── scripts/     pipeline numerado (00–17): .sh = QIIME2, .R = análisis estadístico
├── env/         entornos: imágenes Docker, conda de QIIME2 y versiones de R
├── data/        entradas restringidas (no incluidas); ver data/README.md
├── LICENSE
└── README.md
```

Los scripts esperan, en la raíz del repositorio, las carpetas que el usuario debe
crear al disponer de los datos: `data/` (entradas), `metadata/`, `Reads/` (FASTQ),
`qiime2_output/` (artefactos QIIME2) y `R_output/` (resultados, se crea sola). Todas
están en `.gitignore`.

## Requisitos de software

QIIME2 y PICRUSt2 se ejecutaron mediante **Docker** (plataforma `linux/amd64`); las
imágenes y el entorno de R están documentados en [`env/`](env/) (`docker-images.txt`,
`r-session.txt` y el YAML de conda).

- **QIIME2** 2024.10 — imagen `quay.io/qiime2/amplicon:2024.10`. Como alternativa sin
  Docker, el entorno conda se reconstruye con
  `conda env create -f env/qiime2-amplicon-2024.10-py310-osx-conda.yml`.
- **PICRUSt2** 2.5.3 — imagen `quay.io/biocontainers/picrust2:2.5.3` (ver `10_picrust2.sh`).
- **R** 4.5.2 con los paquetes (versiones exactas en `env/r-session.txt`):
  - Microbioma: `phyloseq`, `qiime2R`, `microbiome`, `vegan`, `Maaslin2`.
  - Modelos y diagnóstico: `broom`, `car`, `sandwich`, `lmtest`, `logistf`, `boot`, `pROC`.
  - Manejo de datos y figuras: `tidyverse` (`dplyr`, `tidyr`, `ggplot2`…), `readxl`,
    `writexl`, `corrplot`, `RColorBrewer`, `viridis`, `scales`, `patchwork`,
    `knitr`, `kableExtra`.

## Orden de ejecución

Los scripts de R se ejecutan desde la carpeta `scripts/` (las rutas son relativas a
ella). El pipeline tiene una dependencia iterativa entre el bloque clínico y el
ómico a través de los **residuos clínicos**, de modo que el orden de ejecución no
coincide con el orden numérico de los ficheros.

1. **QIIME2 (terminal, entorno conda activo):**
   `00 → 00b → 01 → 02 → 02b → 03 → 04 → 05 → 06_diversidad → 07`
2. **Microbiota en R:** `06_phyloseq_import.R → 08_composicion.R` (genera
   `R_output/models/ps_clean.rds`, objeto central del análisis).
3. **Fase 1 clínica:** `14b_modelos_clinico_univariante.R →
   14c_modelos_clinico_multivariante.R`. Produce `residuos_clinicos_fase1.rds`, que
   **debe existir antes** de los bloques ómicos siguientes.
4. **Abundancia diferencial:** `09_maaslin2.R`, `09b`, `09c`.
5. **Función predicha:** `10_picrust2.sh` (Docker) → `11_funcional.R`.
6. **Transmisión vertical:** `12_transmision_vertical.R`.
7. **Epigenética e integración:** `13a → 13b → 13c`.
8. **Ensamblado multiómico:** `14_exportar_features_a_wide.R`.
9. **Fase 2 ómica:** `15_modelos_univariante.R → 16_modelos_multivariante.R`.
10. **Riesgo de sobrepeso:** `17_regresion_logistica.R`.

## Descripción de los scripts

| Script | Qué hace |
|---|---|
| `00_importar_reads.sh` | Importa los FASTQ pareados a QIIME2 (Casava). |
| `00b_cutadapt.sh` | Recorta los cebadores V3–V4 (341F/805R). |
| `01_dada2.sh` | Denoising con DADA2 (ASVs). |
| `02_filter_samples.sh` | Retiene basal + 6 meses; descarta otros tiempos y controles. |
| `02b_prepare_silva.sh` | Prepara la referencia SILVA 138 para la región V3–V4. |
| `03_taxonomia.sh` | Asignación taxonómica por consenso (vsearch). |
| `04_filter_taxonomy.sh` | Elimina mitocondrias, cloroplastos y ASVs sin filo. |
| `05_arbol_filogenetico.sh` | Alineamiento (MAFFT) y árbol (FastTree) enraizado. |
| `06_diversidad.sh` | Rarefacción y métricas de diversidad alfa/beta. |
| `07_estadistica_diversidad.sh` | Kruskal–Wallis (alfa) y PERMANOVA (beta). |
| `06_phyloseq_import.R` | Importa los artefactos QIIME2 a phyloseq. |
| `08_composicion.R` | Composición, figuras descriptivas y `ps_clean.rds`. |
| `09_maaslin2.R` | Abundancia diferencial de géneros (MaAsLin2) por estrato. |
| `09b_maaslin2_prev10.R` | Igual con prevalencia mínima del 10 %. |
| `09c_validar_hits.R` | Validación de los hits a prevalencia 15 % y 20 %. |
| `10_picrust2.sh` | Predicción funcional con PICRUSt2 (Docker). |
| `11_funcional.R` | Abundancia diferencial funcional (KO y vías MetaCyc). |
| `12_transmision_vertical.R` | Transmisión vertical madre→hijo (distancias y ASVs compartidas). |
| `13a_limpiar_metilacion.R` | Depuración de la pirosecuenciación de metilación. |
| `13b_integracion_metilacion.R` | Correlaciones microbiota × metilación × adiposidad. |
| `13c_plots_integracion.R` | Figuras de la integración. |
| `14_exportar_features_a_wide.R` | Ensambla las features ómicas a una fila por díada. |
| `14b_modelos_clinico_univariante.R` | Fase 1 clínica: cribado univariante. |
| `14c_modelos_clinico_multivariante.R` | Fase 1 clínica: multivariante; genera los residuos. |
| `15_modelos_univariante.R` | Cribado univariante de candidatas (clínicas y ómicas). |
| `16_modelos_multivariante.R` | Modelo multivariante de la adiposidad infantil. |
| `17_regresion_logistica.R` | Regresión logística penalizada (Firth) del riesgo de sobrepeso. |

## Datos y reproducibilidad

Con los datos restringidos colocados en sus carpetas (ver `data/README.md`), el
pipeline se reproduce íntegramente siguiendo el orden anterior. Sin ellos, el
repositorio documenta de forma completa el flujo de análisis y los parámetros
empleados.

