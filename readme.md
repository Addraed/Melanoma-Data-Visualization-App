# 🧬 Melanoma Subtype Classifier ( TFM Bioinformática y Bioestadística )

> Desarrollo de una aplicación web para la identificación de síntomas específicos de distintos subtipos de cáncer de melanoma para una mejor profilaxis y un tratamiento más especializado.

**Autor:** Manuel Rosario Marín Fernández  
**Máster:** MU Bioinformática y Bioestadística / Universitat Oberta de Catalunya (UOC) / Universitat de Barcelona  
**Tutora:** Romina Astrid Rebrij  
**Fecha:** Enero 2023  

---

## 📋 Descripción

Este proyecto aplica técnicas de **aprendizaje automático no supervisado** sobre datos genómicos del proyecto [TCGA-SKCM](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) para identificar patrones y dependencias entre variables clínicas y genómicas en pacientes con cáncer de melanoma cutáneo.

A partir de los patrones extraídos, se desarrolló una **aplicación web interactiva con R Shiny** orientada a profesionales sanitarios sin conocimientos avanzados de bioinformática, facilitando la exploración visual de las reglas de asociación obtenidas.

---

## 🎯 Objetivos

- Extraer y analizar patrones e interacciones relevantes en datos multi-ómicos de melanoma
- Identificar perfiles clínicos para los subtipos de mutación: **BRAF**, **RAS** y **Triple WT**
- Desarrollar una aplicación web accesible para la visualización de los resultados

---

## 🗂️ Estructura del repositorio

```
├── data/
│   ├── reglas_asoc.csv               # Reglas de asociación extraídas
│   ├── lista_items.csv               # Ítems (antecedentes y consecuentes)
│   └── matriz_adyacencia_porgrupos.csv  # Matriz para el diagrama de cuerdas
│
├── notebook/
│   └── TFM_code_raw_ML_Melanoma.ipynb   # Pipeline completo de análisis en Python
│
├── shiny_app/
│   └── app.R                         # Código fuente de la aplicación Shiny
│
└── memoria/
    └── TFM_ManuelRosarioMarinFernandez.pdf  # Memoria del TFM
```

---

## ⚙️ Metodología

### Datos
- **Fuente:** Proyecto TCGA-SKCM — 331 pacientes, datos de expresión génica (RNA-seq), epigenómica y variables clínico-patológicas
- **Preprocesamiento:** Normalización log2, centrado de mediana, filtrado de genes planos por varianza
- **Variables de interés:** `MUTATIONSUBTYPES`, `MIRCluster`, `UV_signature`, `LYMPHOCYTE.SCORE`, `RNASEQ-CLUSTER_CONSENHIER`, `MethTypes201408`

### Análisis
- Algoritmo **FP-Growth** (mlxtend) con soporte mínimo de 0.015
- Extracción de **965 itemsets frecuentes** → **317 reglas de asociación** (confianza ≥ 0.85)
- Filtrado por **Lift > 1** y **Leverage > 0** → **11 reglas de interés clínico**

### Resultados destacados

| Subtipo | Características asociadas |
|---|---|
| **BRAF** | Cluster MITF-low, LScore bajo, hiper/hipometilación CpG, MIR type 1 y 4 |
| **RAS** | LScore bajo, islas CpG metiladas, sobreexpresión genes inmunes |
| **Triple WT** | Ausencia firma UV, hipermetilación CpG, cluster queratina, MIR type 3 |

---

## 🖥️ Aplicación web

La app desarrollada con **R Shiny** permite explorar interactivamente las reglas de asociación extraídas. Está dividida en 4 secciones:

- **Inicio — TFM:** Descripción del proyecto y contexto
- **Diagrama de cuerdas:** Visualización interactiva de dependencias entre variables
- **Variables explicadas:** Descripción de variables y medidas de calidad
- **Reglas desglosadas:** Tabla filtrable de reglas por subtipo de mutación

🔗 **App desplegada:** [ShinyApps.io](http://0dt2j0-manuel-mar0n0fern0ndez.shinyapps.io/TFMVisualizing_ShinyApp)

---

## 🛠️ Stack tecnológico

| Categoría | Herramientas |
|---|---|
| **Lenguajes** | Python, R |
| **ML / Análisis** | mlxtend (FP-Growth), pandas, NumPy |
| **Visualización** | matplotlib, plotly, chorddiag |
| **App web** | R Shiny, semantic.dashboard, rpy2 |
| **Entornos** | Google Colab, RStudio |
| **Datos** | TCGAbiolinks, RSEM |

---

## 📦 Instalación y uso

### Requisitos Python
```bash
pip install pandas numpy matplotlib plotly mlxtend
```

### Requisitos R
```r
install.packages(c("shiny", "shinyWidgets", "semantic.dashboard", 
                   "chorddiag", "DT", "dplyr"))
```

### Ejecutar la app localmente
```r
shiny::runApp("shiny_app/app.R")
```

---

## 📄 Referencia principal

Cancer Genome Atlas Network. *Genomic Classification of Cutaneous Melanoma.* Cell, vol. 161, Issue 7, (2015) 1681–1696. DOI: [10.1016/j.cell.2015.05.044](https://doi.org/10.1016/j.cell.2015.05.044)

---

## 📜 Licencia

Este trabajo está sujeto a una licencia **CC BY-NC-ND 3.0 ES** (Reconocimiento - NoComercial - SinObraDerivada).  
© Manuel Rosario Marín Fernández, 2023.
