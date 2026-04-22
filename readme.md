# 🧬 Melanoma Subtype Classifier

> Identificación de perfiles moleculares específicos de los subtipos de melanoma cutáneo mediante reglas de asociación sobre datos multi-ómicos TCGA-SKCM.

**Autor:** Manuel Rosario Marín Fernández  
**Máster:** MU Bioinformática y Bioestadística · UOC / Universitat de Barcelona  
**Tutora:** Romina Astrid Rebrij  
**TFM original:** Enero 2023  

---

## 📋 Descripción

Análisis multi-ómico sobre el dataset **TCGA-SKCM** (331 pacientes, melanoma cutáneo) aplicando el algoritmo FP-Growth para extraer reglas de asociación que identifican perfiles clínicos y moleculares por subtipo de mutación: **BRAF**, **RAS**, **NF1** y **Triple WT**.

El proyecto ha evolucionado desde una aplicación R Shiny de visualización estática hacia una **aplicación web completa** con pipeline ejecutable en tiempo real, visualización interactiva de reglas y documentación clínica integrada.

---

## 🗂️ Estructura del repositorio

```
├── notebook/
│   └── TFM_code_raw_ML_Melanoma.ipynb    # Pipeline original (Google Colab)
│
├── shiny_app/
│   └── app.R                             # App Shiny original del TFM (2023)
│
├── web_app/                              # ← Aplicación web actual
│   ├── main.py                           #   Backend FastAPI
│   ├── requirements.txt
│   ├── Dockerfile
│   └── static/
│       └── index.html                    #   Frontend interactivo
│
├── TFM_ManuelRosarioMarinFernandez.pdf   # Memoria del TFM
├── render.yaml                           # Despliegue en Render
└── README.md
```

---

## 🔄 Evolución del proyecto

### V1. App Shiny (TFM 2023)

La aplicación original se desarrolló en R con Shiny como entregable del TFM. Permitía explorar visualmente los resultados pre-calculados del análisis: las 317 reglas de asociación, la tabla de los 11 perfiles clínicos y un diagrama de cuerdas por subtipos.

**Limitaciones:** los resultados estaban pre-calculados y embebidos por lo que no era posible modificar parámetros, ejecutar el pipeline ni cargar datos nuevos desde la interfaz.

🔗 [ShinyApps.io](http://0dt2j0-manuel-mar0n0fern0ndez.shinyapps.io/TFMVisualizing_ShinyApp)

---

### V2. Web App completa (FastAPI + HTML)

Rediseño completo de la aplicación con el objetivo de hacer el pipeline **ejecutable en tiempo real** desde el navegador, sin depender de R ni de resultados pre-calculados.

**Qué cambia:**

- El pipeline FP-Growth se ejecuta en el servidor (Python/mlxtend) con parámetros configurables
- Los datos se cargan desde la misma fuente del TFM original (Google Drive) y se cachean
- Los resultados calculados se comparan automáticamente con las 11 reglas de la Tabla 3 de la memoria
- Se añade una visualización de grafo de red force-directed interactivo
- Se incluye documentación clínica y estadística integrada en la propia app

🔗 [melanoma-pipeline.onrender.com](https://melanoma-pipeline.onrender.com)

---

## ⚙️ Pipeline

```
Datos TCGA-SKCM (331 pacientes)
    → Limpieza y normalización
    → One-hot encoding de 6 variables multi-ómicas
    → FP-Growth (mlxtend, min_support=0.015)
    → Reglas de asociación (min_confidence=0.85) → 317 reglas
    → Filtrado: lift>1 · leverage>0 · consecuente=MUTATIONSUBTYPES
    → 11 reglas clínicas (Tabla 3, memoria TFM)
```

### Variables del análisis

| Variable | Valores | Descripción |
|---|---|---|
| `MUTATIONSUBTYPES` | BRAF · RAS · NF1 · Triple WT | Subtipo de mutación ( variable objetivo) |
| `UV-signature` | UV signature · not UV | Firma mutacional por exposición UV |
| `RNASEQ-CLUSTER_CONSENHIER` | MITF-low · keratin · immune | Clúster de expresión génica (RNA-seq) |
| `MethTypes.201408` | normal-like · CpG island · hyper · hypo | Tipo de metilación CpG |
| `MIRCluster` | MIR.type.1–4 | Clúster de expresión microRNA |
| `LYMPHOCYTE.SCORE` | 0.0 – 6.0 | Densidad de linfocitos infiltrantes (TIL) |

---

## 🖥️ Funcionalidades de la web app

| Paso | Descripción |
|---|---|
| **1 · Carga de datos** | Datos TFM (Google Drive) · CSV propio · Datos sintéticos calibrados |
| **2 · Preprocesamiento** | One-hot encoding · LYMPHOCYTE.SCORE como categórica 0.0–6.0 |
| **3–5 · Pipeline** | FP-Growth + reglas + filtrado clínico, parámetros configurables |
| **6 · Resultados** | Scatter plot · Reglas TFM · Reglas calculadas · Coincidencias |
| **7 · Grafo** | Red force-directed interactiva · filtros por lift y subtipo |
| **8 · Documentación** | Variables clínicas · métricas estadísticas · contexto del TFM |

---

## 📦 Stack tecnológico

| Capa | V1 (Shiny) | V2 (Web App) |
|---|---|---|
| Datos | TCGAbiolinks (R) | Google Drive (URL TFM) + caché CSV |
| Pipeline ML | mlxtend en Colab (pre-calculado) | mlxtend ejecutado en servidor en tiempo real |
| Backend | R Shiny Server | FastAPI (Python) |
| Frontend | Shiny UI (R) | HTML/JS vanilla |
| Visualización | Diagrama de cuerdas (R) | Scatter plot + Grafo force-directed (Canvas) |
| Despliegue | ShinyApps.io | Render (Docker) |

---

## 🚀 Despliegue local

```bash
git clone https://github.com/Addraed/Melanoma-Data-Visualization-App.git
cd Melanoma-Data-Visualization-App/web_app

pip install -r requirements.txt
uvicorn main:app --reload --port 8000
# → http://localhost:8000
```

---

## 📄 Referencia principal

Cancer Genome Atlas Network. *Genomic Classification of Cutaneous Melanoma.* Cell, vol. 161, Issue 7, 1681–1696, 2015. DOI: [10.1016/j.cell.2015.05.044](https://doi.org/10.1016/j.cell.2015.05.044)

---

## 📜 Licencia

**CC BY-NC-ND 3.0 ES** : Reconocimiento · NoComercial · SinObraDerivada  
© Manuel Rosario Marín Fernández, 2023.
