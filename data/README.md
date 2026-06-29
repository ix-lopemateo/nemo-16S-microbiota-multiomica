# Datos (no incluidos)

Por confidencialidad, este repositorio **no contiene datos**. Al tratarse de una
cohorte humana sujeta a protección de datos y consentimiento informado, las
lecturas crudas y los datos clínicos, antropométricos y de metilación quedan
disponibles **bajo petición razonada y acceso restringido**.

Para reproducir el análisis, coloca aquí (y en las carpetas hermanas) los ficheros
de entrada que esperan los scripts:

- `data/nemo_wide_clean_con_leyenda.xlsx` — datos clínicos y antropométricos (entrada
  de la Fase 1 clínica y del ensamblado multiómico).
- `data/nemo_wide_clean_con_leyenda_multiomics.xlsx` — wide clínico+ómico (lo genera
  `14_exportar_features_a_wide.R`).
- `data/pirosecuenciacion_metilacion_crudo.xlsx` — pirosecuenciación de metilación
  (entrada de `13a`).
- `../metadata/sample-metadata-basal-6m.tsv` — metadatos de muestra para QIIME2 y phyloseq.
- `../Reads/R1/`, `../Reads/R2/` — FASTQ pareados (entrada de `00_importar_reads.sh`).

Las carpetas `qiime2_output/` y `R_output/` se generan al ejecutar el pipeline.
