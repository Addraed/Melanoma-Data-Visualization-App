"""
Melanoma Subtype Classifier — Backend FastAPI
Ejecuta R via subprocess (Rscript) + datos sintéticos como fallback
"""

import os
import random
import subprocess
import logging
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse, JSONResponse
from pydantic import BaseModel

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(title="Melanoma Pipeline API", version="1.0.0")

CACHE_CSV  = Path("/data_cache/skcm_subtipos.csv")
SCRIPT_DIR = Path("/scripts")
SCRIPT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# SCRIPT R — exporta subtipos SKCM a CSV con mapeo flexible
# ============================================================
R_EXPORT_SCRIPT = SCRIPT_DIR / "export_skcm.R"
R_EXPORT_SCRIPT.write_text(r"""
library(TCGAbiolinks)
skcm <- TCGAquery_subtype(tumor = 'skcm')

cat('=== Columnas disponibles en TCGAbiolinks ===\n')
cat(paste(colnames(skcm), collapse=', '), '\n\n')

# Mapeo flexible: nombre estándar TFM -> posibles aliases en TCGAbiolinks
col_map <- list(
  'MUTATIONSUBTYPES'          = c('MUTATIONSUBTYPES','Mutation.Subtype','mutation_subtype'),
  'UV-signature'              = c('UV-signature','UV.signature','UV_signature','uv_signature',
                                   'UV.Signature','UVsignature'),
  'RNASEQ-CLUSTER_CONSENHIER' = c('RNASEQ-CLUSTER_CONSENHIER','RNASEQ.CLUSTER_CONSENHIER',
                                   'RNASEQ_CLUSTER_CONSENHIER','RNAseq_cluster',
                                   'RNASEQ.CLUSTER.CONSENHIER'),
  'MethTypes.201408'          = c('MethTypes.201408','MethTypes201408','MethTypes_201408',
                                   'Meth.type','meth_type'),
  'MIRCluster'                = c('MIRCluster','MIR.cluster','miRNA_cluster','MIR_cluster'),
  'LYMPHOCYTE.SCORE'          = c('LYMPHOCYTE.SCORE','Lymphocyte.Score','LYMPHOCYTE_SCORE',
                                   'lymphocyte_score','LymphocyteScore')
)

# Busqueda insensible a separadores como ultimo recurso
norm_name <- function(s) tolower(gsub('[-._]','',s))

out <- data.frame(row.names=seq_len(nrow(skcm)))
for (standard in names(col_map)) {
  found <- NULL
  # Buscar por alias exacto
  for (alias in col_map[[standard]]) {
    if (alias %in% colnames(skcm)) { found <- alias; break }
  }
  # Buscar por nombre normalizado si no se encontro
  if (is.null(found)) {
    for (col in colnames(skcm)) {
      if (norm_name(col) == norm_name(standard)) { found <- col; break }
    }
  }
  if (!is.null(found)) {
    out[[standard]] <- skcm[[found]]
    cat('OK:', found, '->', standard, '\n')
  } else {
    out[[standard]] <- NA
    cat('FALTA:', standard, '(no encontrada)\n')
  }
}

dir.create('/data_cache', showWarnings=FALSE, recursive=TRUE)
write.csv(out, '/data_cache/skcm_subtipos.csv', row.names=FALSE)
cat('\nExportado:', nrow(out), 'pacientes\n')
cat('Columnas:', paste(colnames(out), collapse=', '), '\n')
""")

# ============================================================
# REGLAS TFM — Tabla 3 de la memoria
# ============================================================
TFM_RULES = [
    {"id":1,"ant":["LYMPHOCYTE.SCORE=0.0","MethTypes.201408=hypo-methylated","RNASEQ-CLUSTER_CONSENHIER=MITF-low"],"cons":["MUTATIONSUBTYPES=BRAF_Hotspot_Mutants"],"conf":1.00,"sup":0.02,"lift":2.22,"lvg":0.01},
    {"id":2,"ant":["MIRCluster=MIR.type.1","RNASEQ-CLUSTER_CONSENHIER=MITF-low"],"cons":["MUTATIONSUBTYPES=BRAF_Hotspot_Mutants"],"conf":0.86,"sup":0.04,"lift":1.90,"lvg":0.02},
    {"id":3,"ant":["MIRCluster=MIR.type.1","MethTypes.201408=hypo-methylated","RNASEQ-CLUSTER_CONSENHIER=MITF-low"],"cons":["MUTATIONSUBTYPES=BRAF_Hotspot_Mutants"],"conf":1.00,"sup":0.02,"lift":2.22,"lvg":0.01},
    {"id":4,"ant":["MIRCluster=MIR.type.4","MethTypes.201408=hyper-methylated","RNASEQ-CLUSTER_CONSENHIER=MITF-low"],"cons":["MUTATIONSUBTYPES=BRAF_Hotspot_Mutants"],"conf":1.00,"sup":0.03,"lift":2.22,"lvg":0.02},
    {"id":5,"ant":["LYMPHOCYTE.SCORE=2.0","MethTypes.201408=hyper-methylated","RNASEQ-CLUSTER_CONSENHIER=MITF-low"],"cons":["MUTATIONSUBTYPES=BRAF_Hotspot_Mutants"],"conf":0.86,"sup":0.02,"lift":1.94,"lvg":0.01},
    {"id":6,"ant":["LYMPHOCYTE.SCORE=0.0","MethTypes.201408=CpG island-methylated","RNASEQ-CLUSTER_CONSENHIER=immune"],"cons":["MUTATIONSUBTYPES=RAS_Hotspot_Mutants"],"conf":1.00,"sup":0.02,"lift":3.62,"lvg":0.02},
    {"id":7,"ant":["LYMPHOCYTE.SCORE=0.0","MIRCluster=MIR.type.4","MethTypes.201408=CpG island-methylated"],"cons":["MUTATIONSUBTYPES=RAS_Hotspot_Mutants"],"conf":0.86,"sup":0.02,"lift":3.10,"lvg":0.01},
    {"id":8,"ant":["LYMPHOCYTE.SCORE=0.0","MethTypes.201408=hyper-methylated","UV-signature=not UV"],"cons":["MUTATIONSUBTYPES=Triple_WT"],"conf":0.86,"sup":0.02,"lift":6.20,"lvg":0.02},
    {"id":9,"ant":["MIRCluster=MIR.type.3","RNASEQ-CLUSTER_CONSENHIER=keratin","UV-signature=not UV"],"cons":["MUTATIONSUBTYPES=Triple_WT"],"conf":0.86,"sup":0.02,"lift":6.33,"lvg":0.02},
    {"id":10,"ant":["MIRCluster=MIR.type.3","MethTypes.201408=hyper-methylated","UV-signature=not UV"],"cons":["MUTATIONSUBTYPES=Triple_WT"],"conf":1.00,"sup":0.02,"lift":7.24,"lvg":0.01},
    {"id":11,"ant":["LYMPHOCYTE.SCORE=6.0","UV-signature=not UV"],"cons":["MUTATIONSUBTYPES=Triple_WT"],"conf":0.86,"sup":0.02,"lift":6.33,"lvg":0.02},
]

# ============================================================
# DATOS SINTÉTICOS — basados en distribuciones del TFM + 11 reglas
# ============================================================
def generate_synthetic_tcga(n: int = 331, noise_pct: float = 0.08) -> pd.DataFrame:
    """
    Genera datos sintéticos TCGA-SKCM con las distribuciones del paper Cell 2015
    y correlaciones de las 11 reglas del TFM.
    """
    profiles = [
        # BRAF — MITF-low + hypo-methylated + MIR.type.1 + LScore=0 (reglas 1,3)
        {"w":0.12,"d":{"MUTATIONSUBTYPES":"BRAF_Hotspot_Mutants","UV-signature":"UV signature","RNASEQ-CLUSTER_CONSENHIER":"MITF-low","MethTypes.201408":"hypo-methylated","MIRCluster":"MIR.type.1","LYMPHOCYTE.SCORE":[0,0]}},
        # BRAF — MITF-low + hyper-methylated + MIR.type.4 (regla 4)
        {"w":0.10,"d":{"MUTATIONSUBTYPES":"BRAF_Hotspot_Mutants","UV-signature":"UV signature","RNASEQ-CLUSTER_CONSENHIER":"MITF-low","MethTypes.201408":"hyper-methylated","MIRCluster":"MIR.type.4","LYMPHOCYTE.SCORE":[2,2]}},
        # BRAF — MITF-low + hyper-methylated + LScore=2 (regla 5)
        {"w":0.08,"d":{"MUTATIONSUBTYPES":"BRAF_Hotspot_Mutants","UV-signature":"UV signature","RNASEQ-CLUSTER_CONSENHIER":"MITF-low","MethTypes.201408":"hyper-methylated","MIRCluster":"MIR.type.1","LYMPHOCYTE.SCORE":[2,2]}},
        # BRAF — keratin variado
        {"w":0.22,"d":{"MUTATIONSUBTYPES":"BRAF_Hotspot_Mutants","UV-signature":"UV signature","RNASEQ-CLUSTER_CONSENHIER":"keratin","MethTypes.201408":"normal-like","MIRCluster":"MIR.type.1","LYMPHOCYTE.SCORE":[0,3]}},
        # RAS — immune + CpG island + LScore=0 (regla 6)
        {"w":0.12,"d":{"MUTATIONSUBTYPES":"RAS_Hotspot_Mutants","UV-signature":"UV signature","RNASEQ-CLUSTER_CONSENHIER":"immune","MethTypes.201408":"CpG island-methylated","MIRCluster":"MIR.type.2","LYMPHOCYTE.SCORE":[0,0]}},
        # RAS — CpG island + MIR.type.4 + LScore=0 (regla 7)
        {"w":0.08,"d":{"MUTATIONSUBTYPES":"RAS_Hotspot_Mutants","UV-signature":"UV signature","RNASEQ-CLUSTER_CONSENHIER":"keratin","MethTypes.201408":"CpG island-methylated","MIRCluster":"MIR.type.4","LYMPHOCYTE.SCORE":[0,0]}},
        # NF1
        {"w":0.14,"d":{"MUTATIONSUBTYPES":"NF1_Any_Mutants","UV-signature":"UV signature","RNASEQ-CLUSTER_CONSENHIER":"immune","MethTypes.201408":"hyper-methylated","MIRCluster":"MIR.type.2","LYMPHOCYTE.SCORE":[3,6]}},
        # Triple WT — not UV + hyper-methylated + LScore=0 (regla 8)
        {"w":0.04,"d":{"MUTATIONSUBTYPES":"Triple_WT","UV-signature":"not UV","RNASEQ-CLUSTER_CONSENHIER":"keratin","MethTypes.201408":"hyper-methylated","MIRCluster":"MIR.type.3","LYMPHOCYTE.SCORE":[0,0]}},
        # Triple WT — not UV + keratin + MIR.type.3 (regla 9)
        {"w":0.04,"d":{"MUTATIONSUBTYPES":"Triple_WT","UV-signature":"not UV","RNASEQ-CLUSTER_CONSENHIER":"keratin","MethTypes.201408":"CpG island-methylated","MIRCluster":"MIR.type.3","LYMPHOCYTE.SCORE":[1,4]}},
        # Triple WT — not UV + hyper-methylated + MIR.type.3 (regla 10)
        {"w":0.04,"d":{"MUTATIONSUBTYPES":"Triple_WT","UV-signature":"not UV","RNASEQ-CLUSTER_CONSENHIER":"keratin","MethTypes.201408":"hyper-methylated","MIRCluster":"MIR.type.3","LYMPHOCYTE.SCORE":[1,3]}},
        # Triple WT — not UV + LScore=6 (regla 11)
        {"w":0.02,"d":{"MUTATIONSUBTYPES":"Triple_WT","UV-signature":"not UV","RNASEQ-CLUSTER_CONSENHIER":"immune","MethTypes.201408":"hyper-methylated","MIRCluster":"MIR.type.3","LYMPHOCYTE.SCORE":[6,6]}},
    ]

    all_vals = {
        "MUTATIONSUBTYPES":          ["BRAF_Hotspot_Mutants","RAS_Hotspot_Mutants","Triple_WT","NF1_Any_Mutants"],
        "UV-signature":              ["UV signature","not UV"],
        "RNASEQ-CLUSTER_CONSENHIER": ["keratin","immune","MITF-low"],
        "MethTypes.201408":          ["normal-like","CpG island-methylated","hypo-methylated","hyper-methylated"],
        "MIRCluster":                ["MIR.type.1","MIR.type.2","MIR.type.3","MIR.type.4"],
    }

    # Normalizar pesos
    total_w = sum(p["w"] for p in profiles)
    rows = []
    for _ in range(n):
        rv, acc = random.random() * total_w, 0
        prof = profiles[-1]
        for p in profiles:
            acc += p["w"]
            if rv < acc:
                prof = p
                break

        row = {}
        noisy = random.random() < noise_pct
        for k, v in prof["d"].items():
            if k == "LYMPHOCYTE.SCORE":
                lo, hi = v
                row[k] = lo if lo == hi else random.randint(lo, hi)
                if noisy:
                    row[k] = random.randint(0, 6)
            else:
                row[k] = random.choice(all_vals[k]) if noisy else v
        rows.append(row)

    return pd.DataFrame(rows)


# ============================================================
# UTILIDADES
# ============================================================
def run_rscript(script_path: Path, timeout: int = 120):
    rscript_bin = "Rscript"
    for p in ["/usr/local/bin/Rscript", "/usr/bin/Rscript"]:
        if Path(p).exists():
            rscript_bin = p
            break
    try:
        result = subprocess.run(
            [rscript_bin, "--vanilla", str(script_path)],
            capture_output=True, text=True, timeout=timeout,
            env={**os.environ, "R_HOME": "/usr/local/lib/R"}
        )
        return result.returncode == 0, result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        return False, "Timeout"
    except FileNotFoundError:
        return False, "Rscript no encontrado"
    except Exception as e:
        return False, str(e)


# URL del CSV original del TFM (Google Drive del notebook)
GDRIVE_URL = "https://drive.google.com/uc?export=download&id=1VF3MqQ3J7GBc527ClL2IgJ23XTc5-hX7"

# Columnas requeridas del TFM
REQUIRED_COLS = {"MUTATIONSUBTYPES","UV-signature","RNASEQ-CLUSTER_CONSENHIER",
                 "MethTypes.201408","MIRCluster","LYMPHOCYTE.SCORE"}


def clean_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Limpia el DataFrame: elimina guiones, NAs y filas sin subtipo."""
    # Reemplazar "-", "NA", "" por NaN
    df = df.replace({'-': None, 'NA': None, '': None})
    # Eliminar filas sin subtipo de mutación (variable objetivo)
    df = df.dropna(subset=['MUTATIONSUBTYPES'])
    # LYMPHOCYTE.SCORE: asegurar numérico
    if 'LYMPHOCYTE.SCORE' in df.columns:
        df['LYMPHOCYTE.SCORE'] = pd.to_numeric(df['LYMPHOCYTE.SCORE'], errors='coerce').fillna(0)
    logger.info(f"DataFrame limpiado: {len(df)} filas válidas")
    return df.reset_index(drop=True)


def load_skcm_data() -> pd.DataFrame:
    """
    Carga los datos TCGA-SKCM del TFM.
    Orden de prioridad:
    1. Caché CSV local (si existe y está completo)
    2. Google Drive — URL del notebook original del TFM
    """
    # 1. Intentar caché local
    if CACHE_CSV.exists():
        df = pd.read_csv(CACHE_CSV)
        missing = REQUIRED_COLS - set(df.columns)
        if not missing:
            df = clean_dataframe(df)
            logger.info(f"CSV cargado desde caché: {len(df)} pacientes")
            return df
        logger.warning(f"Caché incompleto (faltan: {missing}), descargando de nuevo...")
        CACHE_CSV.unlink()

    # 2. Descargar desde Google Drive (misma URL que el notebook del TFM)
    logger.info(f"Descargando datos desde Google Drive (URL del TFM)...")
    try:
        import urllib.request
        CACHE_CSV.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(GDRIVE_URL, str(CACHE_CSV))
        df = pd.read_csv(CACHE_CSV, na_values=['-', 'NA', 'nan', ''])
        missing = REQUIRED_COLS - set(df.columns)
        if missing:
            CACHE_CSV.unlink()
            raise RuntimeError(f"CSV descargado pero faltan columnas: {missing}")
        df = clean_dataframe(df)
        df.to_csv(CACHE_CSV, index=False)  # Guardar versión limpia
        logger.info(f"Descargado OK: {len(df)} pacientes")
        return df
    except Exception as e:
        raise RuntimeError(f"No se pudieron cargar los datos del TFM: {e}")


def run_fpgrowth_pipeline(df: pd.DataFrame,
                           min_support: float = 0.015,
                           min_confidence: float = 0.85,
                           min_lift: float = 1.0,
                           min_leverage: float = 0.0) -> dict:
    from mlxtend.frequent_patterns import fpgrowth, association_rules

    logger.info(f"Pipeline con {len(df)} pacientes, columnas: {list(df.columns)}")

    cat_cols = [c for c in ["MUTATIONSUBTYPES","UV-signature","RNASEQ-CLUSTER_CONSENHIER",
                             "MethTypes.201408","MIRCluster"] if c in df.columns]
    if not cat_cols:
        raise ValueError(f"Sin columnas categóricas. Disponibles: {list(df.columns)}")

    transacciones = pd.get_dummies(df[cat_cols], prefix_sep="=")

    if "LYMPHOCYTE.SCORE" in df.columns:
        lscore = pd.to_numeric(df["LYMPHOCYTE.SCORE"], errors="coerce").fillna(0)
        transacciones["LYMPHOCYTE.SCORE=0.0"] = (lscore == 0).astype(int)
        transacciones["LYMPHOCYTE.SCORE=low"]  = ((lscore >= 1) & (lscore <= 3)).astype(int)
        transacciones["LYMPHOCYTE.SCORE=high"] = (lscore >= 4).astype(int)

    transacciones = transacciones.astype(bool)
    logger.info(f"Transacciones: {transacciones.shape}")

    itemsets = fpgrowth(transacciones, min_support=min_support,
                        use_colnames=True, max_len=4)
    if len(itemsets) == 0:
        return {"n_itemsets":0,"n_rules_total":0,"n_rules_filtered":0,
                "n_rules_clinical":0,"rules_all":[],"rules_clinical":[],
                "params":{"min_support":min_support,"n_patients":len(df)}}

    reglas = association_rules(itemsets, metric="confidence",
                               min_threshold=min_confidence)

    reglas_filtradas = reglas[
        (reglas["lift"] > min_lift) & (reglas["leverage"] > min_leverage)
    ].copy()

    reglas_clinicas = reglas_filtradas[
        reglas_filtradas["consequents"].apply(
            lambda c: any("MUTATIONSUBTYPES" in str(x) for x in c)
        )
    ].copy()

    def to_dict(df_r):
        if len(df_r) == 0:
            return []
        r = df_r.copy()
        r["antecedents"] = r["antecedents"].apply(list)
        r["consequents"] = r["consequents"].apply(list)
        records = r[["antecedents","consequents","support","confidence",
                     "lift","leverage","conviction"]].to_dict(orient="records")
        # Limpiar inf/NaN (conviction = inf cuando confidence = 1.0)
        def clean(v):
            if isinstance(v, float) and (v != v or abs(v) == float('inf')):
                return None
            if isinstance(v, float):
                return round(v, 4)
            return v
        return [{k: clean(val) if not isinstance(val, list) else val
                 for k, val in rec.items()} for rec in records]

    return {
        "n_itemsets": len(itemsets),
        "n_rules_total": len(reglas),
        "n_rules_filtered": len(reglas_filtradas),
        "n_rules_clinical": len(reglas_clinicas),
        "rules_clinical": to_dict(reglas_clinicas.sort_values("lift", ascending=False)),
        "params": {
            "min_support": min_support,
            "min_confidence": min_confidence,
            "min_lift": min_lift,
            "min_leverage": min_leverage,
            "n_patients": len(df),
            "n_cols": len(transacciones.columns),
        }
    }


def match_tfm_rule(rule: dict) -> Optional[int]:
    ant_set = set(rule["antecedents"])
    cons_set = set(rule["consequents"])
    for tfm in TFM_RULES:
        if set(tfm["ant"]) == ant_set and set(tfm["cons"]) == cons_set:
            return tfm["id"]
    return None


# ============================================================
# ENDPOINTS
# ============================================================

class PipelineParams(BaseModel):
    min_support: float = 0.015
    min_confidence: float = 0.85
    min_lift: float = 1.0
    min_leverage: float = 0.0


@app.get("/health")
def health():
    return {"status": "ok", "version": "1.0.0"}


@app.get("/api/tcga-status")
def tcga_status():
    rscript_ok = any(Path(p).exists() for p in ["/usr/local/bin/Rscript", "/usr/bin/Rscript"])
    csv_info = {}
    if CACHE_CSV.exists():
        try:
            df = pd.read_csv(CACHE_CSV)
            csv_info = {"rows": len(df), "columns": list(df.columns)}
        except Exception as e:
            csv_info = {"error": str(e)}
    return {
        "cache_csv": CACHE_CSV.exists(),
        "rscript_available": rscript_ok,
        "csv_info": csv_info,
    }


@app.get("/api/tcga-data")
def get_tcga_data():
    try:
        df = load_skcm_data()
        return JSONResponse({
            "source": "tcgabiolinks_r",
            "n_patients": len(df),
            "columns": list(df.columns),
            "data": df.fillna("").to_dict(orient="records"),
        })
    except Exception as e:
        raise HTTPException(status_code=503, detail=str(e))


@app.get("/api/synthetic-data")
def get_synthetic_data(n: int = 331, noise: float = 0.08):
    """Genera datos sintéticos TCGA-SKCM basados en distribuciones del TFM."""
    df = generate_synthetic_tcga(n=min(n, 1000), noise_pct=min(noise, 0.5))
    return JSONResponse({
        "source": "synthetic",
        "n_patients": len(df),
        "columns": list(df.columns),
        "data": df.to_dict(orient="records"),
    })


@app.post("/api/run-pipeline")
def run_pipeline(params: PipelineParams):
    try:
        df = load_skcm_data()
    except Exception as e:
        raise HTTPException(status_code=503, detail=f"Datos TCGA no disponibles: {e}")

    try:
        result = run_fpgrowth_pipeline(df, **params.model_dump())
    except Exception as e:
        logger.error(f"Error pipeline: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))

    for rule in result["rules_clinical"]:
        rule["tfm_match"] = match_tfm_rule(rule)

    result["tfm_rules"] = TFM_RULES
    matched = sum(1 for r in result["rules_clinical"] if r.get("tfm_match"))
    result["tfm_match_stats"] = {
        "tfm_rules_total": len(TFM_RULES),
        "computed_rules_total": result["n_rules_clinical"],
        "matched": matched,
    }
    return JSONResponse(result)


@app.post("/api/run-pipeline-synthetic")
def run_pipeline_synthetic(params: PipelineParams, n: int = 331, noise: float = 0.08):
    """Ejecuta el pipeline sobre datos sintéticos."""
    df = generate_synthetic_tcga(n=n, noise_pct=noise)
    try:
        result = run_fpgrowth_pipeline(df, **params.model_dump())
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

    for rule in result["rules_clinical"]:
        rule["tfm_match"] = match_tfm_rule(rule)
    result["tfm_rules"] = TFM_RULES
    matched = sum(1 for r in result["rules_clinical"] if r.get("tfm_match"))
    result["tfm_match_stats"] = {
        "tfm_rules_total": len(TFM_RULES),
        "computed_rules_total": result["n_rules_clinical"],
        "matched": matched,
    }
    result["source"] = "synthetic"
    return JSONResponse(result)


@app.post("/api/refresh-cache")
def refresh_cache():
    """Borra el CSV cacheado y lo regenera via TCGAbiolinks."""
    if CACHE_CSV.exists():
        CACHE_CSV.unlink()
        logger.info("Caché eliminado")
    try:
        df = load_skcm_data()
        return JSONResponse({"status": "ok", "rows": len(df), "columns": list(df.columns)})
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/tfm-rules")
def get_tfm_rules():
    return JSONResponse({"rules": TFM_RULES})


# Frontend estático
app.mount("/static", StaticFiles(directory="static"), name="static")

@app.get("/graph", response_class=HTMLResponse)
def serve_graph():
    graph = Path("static/graph.html")
    if graph.exists():
        return HTMLResponse(graph.read_text())
    return HTMLResponse("<h1>Grafo no encontrado</h1>")

@app.get("/", response_class=HTMLResponse)
def serve_frontend():
    index = Path("static/index.html")
    if index.exists():
        return HTMLResponse(index.read_text())
    return HTMLResponse("<h1>Frontend no encontrado</h1>")
