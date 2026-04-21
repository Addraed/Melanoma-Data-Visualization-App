"""
Melanoma Subtype Classifier — Backend FastAPI
Ejecuta R via subprocess (Rscript) en lugar de rpy2
"""

import os
import json
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

# Script R para exportar subtipos SKCM
R_EXPORT_SCRIPT = SCRIPT_DIR / "export_skcm.R"
R_EXPORT_SCRIPT.write_text("""
library(TCGAbiolinks)
skcm <- TCGAquery_subtype(tumor = 'skcm')

cat('Columnas en TCGAbiolinks:\n')
cat(paste(colnames(skcm), collapse=', '), '\n')

col_map <- list(
  'MUTATIONSUBTYPES'          = c('MUTATIONSUBTYPES','Mutation.Subtype'),
  'UV-signature'              = c('UV-signature','UV.signature','UV_signature'),
  'RNASEQ-CLUSTER_CONSENHIER' = c('RNASEQ-CLUSTER_CONSENHIER','RNASEQ.CLUSTER_CONSENHIER','RNASEQ_CLUSTER_CONSENHIER'),
  'MethTypes.201408'          = c('MethTypes.201408','MethTypes201408','MethTypes_201408'),
  'MIRCluster'                = c('MIRCluster','MIR.cluster'),
  'LYMPHOCYTE.SCORE'          = c('LYMPHOCYTE.SCORE','Lymphocyte.Score','LYMPHOCYTE_SCORE')
)

out <- data.frame(row.names=seq_len(nrow(skcm)))
for (standard in names(col_map)) {
  found <- NULL
  for (alias in col_map[[standard]]) {
    if (alias %in% colnames(skcm)) { found <- alias; break }
  }
  if (!is.null(found)) {
    out[[standard]] <- skcm[[found]]
    cat('OK:', found, '->', standard, '\n')
  } else {
    out[[standard]] <- NA
    cat('FALTA:', standard, '\n')
  }
}

dir.create('/data_cache', showWarnings=FALSE, recursive=TRUE)
write.csv(out, '/data_cache/skcm_subtipos.csv', row.names=FALSE)
cat('Exportado:', nrow(out), 'pacientes\n')
""")

# Reglas exactas del TFM — Tabla 3 de la memoria
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


def run_rscript(script_path: Path, timeout: int = 120) -> tuple[bool, str]:
    """Ejecuta un script R via Rscript y retorna (éxito, output)."""
    # Encontrar Rscript en rutas conocidas de rocker
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
        output = result.stdout + result.stderr
        return result.returncode == 0, output
    except subprocess.TimeoutExpired:
        return False, "Timeout ejecutando Rscript"
    except FileNotFoundError:
        return False, "Rscript no encontrado en PATH"
    except Exception as e:
        return False, str(e)


def load_skcm_data() -> pd.DataFrame:
    """Carga el CSV cacheado o lo genera via Rscript."""
    if CACHE_CSV.exists():
        logger.info("Cargando desde caché CSV")
        return pd.read_csv(CACHE_CSV)

    logger.info("Generando datos via Rscript...")
    # Asegurar path absoluto de Rscript (rocker lo instala en /usr/local/bin)
    for rscript_path in ["/usr/local/bin/Rscript", "/usr/bin/Rscript", "Rscript"]:
        if Path(rscript_path).exists() or rscript_path == "Rscript":
            break
    ok, output = run_rscript(R_EXPORT_SCRIPT, timeout=180)
    logger.info(f"Rscript output: {output[:500]}")

    if ok and CACHE_CSV.exists():
        return pd.read_csv(CACHE_CSV)

    raise RuntimeError(f"No se pudo generar el CSV: {output[:200]}")


def run_fpgrowth_pipeline(df: pd.DataFrame,
                           min_support: float = 0.015,
                           min_confidence: float = 0.85,
                           min_lift: float = 1.0,
                           min_leverage: float = 0.0) -> dict:
    """Pipeline FP-Growth equivalente al notebook del TFM."""
    from mlxtend.frequent_patterns import fpgrowth, association_rules

    # Mapear nombres de columna — TCGAbiolinks puede exportar con distintos separadores
    col_aliases = {
        "MUTATIONSUBTYPES":          ["MUTATIONSUBTYPES", "Mutation.Subtype"],
        "UV-signature":              ["UV-signature", "UV.signature", "UV_signature"],
        "RNASEQ-CLUSTER_CONSENHIER": ["RNASEQ-CLUSTER_CONSENHIER", "RNASEQ.CLUSTER_CONSENHIER",
                                      "RNASEQ_CLUSTER_CONSENHIER"],
        "MethTypes.201408":          ["MethTypes.201408", "MethTypes201408", "MethTypes_201408"],
        "MIRCluster":                ["MIRCluster", "MIR.cluster"],
        "LYMPHOCYTE.SCORE":          ["LYMPHOCYTE.SCORE", "Lymphocyte.Score", "LYMPHOCYTE_SCORE"],
    }

    def find_col(target, df_cols):
        if target in df_cols:
            return target
        for alias in col_aliases.get(target, []):
            if alias in df_cols:
                return alias
        norm = lambda s: s.lower().replace("-","").replace(".","").replace("_","")
        for col in df_cols:
            if norm(col) == norm(target):
                return col
        return None

    rename_map = {}
    for standard in col_aliases:
        found = find_col(standard, df.columns)
        if found and found != standard:
            rename_map[found] = standard
    if rename_map:
        df = df.rename(columns=rename_map)
        logger.info(f"Columnas renombradas: {rename_map}")

    logger.info(f"Columnas disponibles: {list(df.columns)}")

    cat_cols = ["MUTATIONSUBTYPES", "UV-signature", "RNASEQ-CLUSTER_CONSENHIER",
                "MethTypes.201408", "MIRCluster"]
    num_col = "LYMPHOCYTE.SCORE"

    cat_cols = [c for c in cat_cols if c in df.columns]
    if not cat_cols:
        raise ValueError(f"No se encontraron columnas. Disponibles: {list(df.columns)}")

    # One-hot encoding
    transacciones = pd.get_dummies(df[cat_cols], prefix_sep="=")
    lscore_col = find_col("LYMPHOCYTE.SCORE", df.columns)
    lscore = pd.to_numeric(df[lscore_col] if lscore_col else pd.Series([0]*len(df)), errors="coerce").fillna(0)
    transacciones["LYMPHOCYTE.SCORE=0.0"] = (lscore == 0).astype(int)
    transacciones["LYMPHOCYTE.SCORE=low"]  = ((lscore >= 1) & (lscore <= 3)).astype(int)
    transacciones["LYMPHOCYTE.SCORE=high"] = (lscore >= 4).astype(int)
    transacciones = transacciones.astype(bool)

    # FP-Growth
    itemsets = fpgrowth(transacciones, min_support=min_support,
                        use_colnames=True, max_len=4)
    reglas = association_rules(itemsets, metric="confidence",
                               min_threshold=min_confidence)

    # Filtrar
    reglas_filtradas = reglas[
        (reglas["lift"] > min_lift) & (reglas["leverage"] > min_leverage)
    ].copy()

    def has_mutation_cons(cons):
        return any("MUTATIONSUBTYPES" in str(c) for c in cons)

    reglas_clinicas = reglas_filtradas[
        reglas_filtradas["consequents"].apply(has_mutation_cons)
    ].copy()

    def fs_to_list(series):
        return series.apply(lambda x: list(x))

    def rules_to_dict(df_r):
        if len(df_r) == 0:
            return []
        r = df_r.copy()
        r["antecedents"] = fs_to_list(r["antecedents"])
        r["consequents"] = fs_to_list(r["consequents"])
        return r[["antecedents","consequents","support","confidence",
                   "lift","leverage","conviction"]].round(4).to_dict(orient="records")

    return {
        "n_itemsets": len(itemsets),
        "n_rules_total": len(reglas),
        "n_rules_filtered": len(reglas_filtradas),
        "n_rules_clinical": len(reglas_clinicas),
        "rules_all": rules_to_dict(reglas_filtradas.head(500)),
        "rules_clinical": rules_to_dict(
            reglas_clinicas.sort_values("lift", ascending=False)
        ),
        "params": {
            "min_support": min_support,
            "min_confidence": min_confidence,
            "min_lift": min_lift,
            "min_leverage": min_leverage,
            "n_patients": len(df),
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
    # Buscar Rscript en ubicaciones conocidas de rocker
    rscript_paths = ["/usr/local/bin/Rscript", "/usr/bin/Rscript"]
    rscript_ok = any(Path(p).exists() for p in rscript_paths)

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
            "source": "cache_csv" if not CACHE_CSV.exists() else "tcgabiolinks_r",
            "n_patients": len(df),
            "columns": list(df.columns),
            "data": df.fillna("").to_dict(orient="records"),
        })
    except Exception as e:
        raise HTTPException(status_code=503, detail=str(e))


@app.post("/api/run-pipeline")
def run_pipeline(params: PipelineParams):
    try:
        df = load_skcm_data()
    except Exception as e:
        raise HTTPException(status_code=503, detail=f"Datos no disponibles: {e}")

    try:
        result = run_fpgrowth_pipeline(
            df,
            min_support=params.min_support,
            min_confidence=params.min_confidence,
            min_lift=params.min_lift,
            min_leverage=params.min_leverage,
        )
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

    # Limpiar inf y NaN antes de serializar (conviction puede ser inf cuando conf=1.0)
    def clean_floats(obj):
        if isinstance(obj, float):
            if obj != obj or obj == float('inf') or obj == float('-inf'):
                return None
            return round(obj, 4)
        if isinstance(obj, dict):
            return {k: clean_floats(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [clean_floats(i) for i in obj]
        return obj

    result = clean_floats(result)
    return JSONResponse(result)


@app.get("/api/tfm-rules")
def get_tfm_rules():
    return JSONResponse({"rules": TFM_RULES})



@app.post("/api/refresh-cache")
def refresh_cache():
    """Fuerza la regeneración del CSV borrando el caché actual."""
    if CACHE_CSV.exists():
        CACHE_CSV.unlink()
        logger.info("Caché CSV eliminado")
    try:
        df = load_skcm_data()
        return JSONResponse({
            "status": "ok",
            "rows": len(df),
            "columns": list(df.columns)
        })
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# Frontend estático
app.mount("/static", StaticFiles(directory="static"), name="static")

@app.get("/", response_class=HTMLResponse)
def serve_frontend():
    index = Path("static/index.html")
    if index.exists():
        return HTMLResponse(index.read_text())
    return HTMLResponse("<h1>Frontend no encontrado</h1>")
