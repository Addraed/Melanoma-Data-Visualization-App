"""
Melanoma Subtype Classifier — Backend FastAPI
Equivalente web a TCGAbiolinks::TCGAquery_subtype("skcm") + pipeline FP-Growth
"""

import os
import json
import time
import logging
import asyncio
from pathlib import Path
from typing import Optional
from functools import lru_cache

import pandas as pd
import numpy as np
from fastapi import FastAPI, HTTPException, BackgroundTasks
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse, JSONResponse, StreamingResponse
from pydantic import BaseModel

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(title="Melanoma Pipeline API", version="1.0.0")

CACHE_PATH = Path("/data_cache/skcm_subtipos.rds")
CACHE_CSV = Path("/data_cache/skcm_subtipos.csv")

# ============================================================
# UTILIDADES R / rpy2
# ============================================================

def get_r_interface():
    """Inicializa rpy2 lazy — solo cuando se necesita."""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr
        pandas2ri.activate()
        return ro, pandas2ri
    except Exception as e:
        logger.error(f"rpy2 no disponible: {e}")
        return None, None


def load_skcm_from_r() -> pd.DataFrame:
    """
    Carga los subtipos SKCM via TCGAbiolinks::TCGAquery_subtype('skcm').
    Equivalente exacto al comando del TFM original.
    Usa caché RDS si está disponible.
    """
    ro, pandas2ri = get_r_interface()
    if ro is None:
        raise RuntimeError("rpy2 no disponible")

    # Intentar cargar desde caché RDS primero
    if CACHE_PATH.exists():
        logger.info("Cargando desde caché RDS...")
        r_code = f"""
        skcm <- readRDS('{CACHE_PATH}')
        skcm
        """
        result = ro.r(r_code)
        df = pandas2ri.rpy2py(result)
        logger.info(f"Cache cargado: {len(df)} pacientes")
        return df

    # Si no hay caché, descargar via TCGAbiolinks
    logger.info("Descargando via TCGAbiolinks::TCGAquery_subtype('skcm')...")
    r_code = f"""
    suppressMessages(library(TCGAbiolinks))
    skcm <- TCGAquery_subtype(tumor = 'skcm')
    dir.create('{CACHE_PATH.parent}', showWarnings=FALSE, recursive=TRUE)
    saveRDS(skcm, '{CACHE_PATH}')
    skcm
    """
    result = ro.r(r_code)
    df = pandas2ri.rpy2py(result)
    logger.info(f"Descargado y cacheado: {len(df)} pacientes")
    return df


def prepare_pipeline_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Selecciona y limpia las 6 variables del TFM.
    Maneja distintos nombres de columna según versión de TCGAbiolinks.
    """
    # Mapeo de posibles nombres de columna a los del TFM
    col_map = {
        "MUTATIONSUBTYPES": ["MUTATIONSUBTYPES", "Mutation.Subtype", "mutation_subtype"],
        "UV-signature": ["UV-signature", "UV.signature", "UV_signature"],
        "RNASEQ-CLUSTER_CONSENHIER": ["RNASEQ-CLUSTER_CONSENHIER", "RNASEQ.CLUSTER_CONSENHIER", "RNAseq_cluster"],
        "MethTypes.201408": ["MethTypes.201408", "MethTypes201408", "Meth.type"],
        "MIRCluster": ["MIRCluster", "MIR.cluster", "miRNA_cluster"],
        "LYMPHOCYTE.SCORE": ["LYMPHOCYTE.SCORE", "Lymphocyte.Score", "lymphocyte_score"],
    }

    result = pd.DataFrame()
    available = list(df.columns)

    for target_col, candidates in col_map.items():
        found = None
        for c in candidates:
            if c in available:
                found = c
                break
        if found:
            result[target_col] = df[found]
        else:
            logger.warning(f"Columna '{target_col}' no encontrada. Cols disponibles: {available[:10]}")
            result[target_col] = None

    # Limpiar valores nulos
    result = result.dropna(subset=["MUTATIONSUBTYPES"])
    result["LYMPHOCYTE.SCORE"] = pd.to_numeric(result["LYMPHOCYTE.SCORE"], errors="coerce").fillna(0)

    return result.reset_index(drop=True)


# ============================================================
# FP-GROWTH + REGLAS DE ASOCIACIÓN (Python, equivalente al TFM)
# ============================================================

def run_fpgrowth_pipeline(df: pd.DataFrame, min_support: float = 0.015,
                           min_confidence: float = 0.85,
                           min_lift: float = 1.0,
                           min_leverage: float = 0.0) -> dict:
    """
    Pipeline completo FP-Growth equivalente al notebook del TFM.
    Retorna itemsets, todas las reglas, y las reglas filtradas.
    """
    from mlxtend.frequent_patterns import fpgrowth, association_rules

    # One-hot encoding — igual que el TFM
    cat_cols = ["MUTATIONSUBTYPES", "UV-signature", "RNASEQ-CLUSTER_CONSENHIER",
                "MethTypes.201408", "MIRCluster"]
    num_col = "LYMPHOCYTE.SCORE"

    transacciones = pd.get_dummies(df[cat_cols], prefix_sep="=")

    # Discretizar LYMPHOCYTE.SCORE en bins: 0.0, low(1-3), high(4+)
    lscore = pd.to_numeric(df[num_col], errors="coerce").fillna(0)
    transacciones["LYMPHOCYTE.SCORE=0.0"] = (lscore == 0).astype(int)
    transacciones["LYMPHOCYTE.SCORE=low"]  = ((lscore >= 1) & (lscore <= 3)).astype(int)
    transacciones["LYMPHOCYTE.SCORE=high"] = (lscore >= 4).astype(int)

    transacciones = transacciones.astype(bool)

    # FP-Growth
    itemsets = fpgrowth(transacciones, min_support=min_support,
                        use_colnames=True, max_len=4)
    logger.info(f"FP-Growth: {len(itemsets)} itemsets frecuentes")

    # Reglas de asociación
    reglas = association_rules(itemsets, metric="confidence",
                               min_threshold=min_confidence)
    logger.info(f"Reglas generadas: {len(reglas)}")

    # Filtrar por lift y leverage
    reglas_filtradas = reglas[
        (reglas["lift"] > min_lift) & (reglas["leverage"] > min_leverage)
    ].copy()

    # Filtrar por consecuente = MUTATIONSUBTYPES (como en el TFM)
    def has_mutation_cons(cons):
        return any("MUTATIONSUBTYPES" in str(c) for c in cons)

    reglas_clinicas = reglas_filtradas[
        reglas_filtradas["consequents"].apply(has_mutation_cons)
    ].copy()

    logger.info(f"Reglas clínicas (MUTATIONSUBTYPES consecuente): {len(reglas_clinicas)}")

    # Convertir frozensets a listas para serialización JSON
    def fs_to_list(series):
        return series.apply(lambda x: list(x))

    def rules_to_dict(df_rules):
        if len(df_rules) == 0:
            return []
        r = df_rules.copy()
        r["antecedents"] = fs_to_list(r["antecedents"])
        r["consequents"] = fs_to_list(r["consequents"])
        return r[["antecedents", "consequents", "support", "confidence",
                   "lift", "leverage", "conviction"]].round(4).to_dict(orient="records")

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
            "n_encoded_cols": len(transacciones.columns),
        }
    }


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


def match_tfm_rule(computed_rule: dict) -> Optional[int]:
    """
    Comprueba si una regla calculada coincide con alguna regla del TFM.
    Retorna el ID de la regla del TFM o None.
    """
    ant_set = set(computed_rule["antecedents"])
    cons_set = set(computed_rule["consequents"])

    for tfm in TFM_RULES:
        if set(tfm["ant"]) == ant_set and set(tfm["cons"]) == cons_set:
            return tfm["id"]
    return None


# ============================================================
# ENDPOINTS API
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
    """Informa si el caché SKCM está disponible."""
    rds_ok = CACHE_PATH.exists()
    csv_ok = CACHE_CSV.exists()
    return {
        "cache_rds": rds_ok,
        "cache_csv": csv_ok,
        "rpy2_available": get_r_interface()[0] is not None,
    }


@app.get("/api/tcga-data")
def get_tcga_data():
    """
    Retorna el dataset TCGA-SKCM completo con las 6 variables del TFM.
    Equivalente a TCGAbiolinks::TCGAquery_subtype('skcm').
    Usa caché si está disponible.
    """
    # Intentar desde caché CSV primero (más rápido)
    if CACHE_CSV.exists():
        df = pd.read_csv(CACHE_CSV)
        return JSONResponse({
            "source": "cache_csv",
            "n_patients": len(df),
            "columns": list(df.columns),
            "data": df.to_dict(orient="records"),
        })

    # Intentar via R / TCGAbiolinks
    try:
        df_raw = load_skcm_from_r()
        df = prepare_pipeline_data(df_raw)
        # Guardar CSV para futuras peticiones
        CACHE_CSV.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(CACHE_CSV, index=False)
        return JSONResponse({
            "source": "tcgabiolinks_r",
            "n_patients": len(df),
            "columns": list(df.columns),
            "data": df.to_dict(orient="records"),
        })
    except Exception as e:
        logger.error(f"Error cargando datos TCGA: {e}")
        raise HTTPException(status_code=503, detail=str(e))


@app.post("/api/run-pipeline")
def run_pipeline(params: PipelineParams):
    """
    Ejecuta el pipeline completo FP-Growth sobre los datos TCGA-SKCM.
    Retorna itemsets, reglas calculadas, reglas clínicas y comparativa con TFM.
    """
    # Cargar datos
    if CACHE_CSV.exists():
        df = pd.read_csv(CACHE_CSV)
    else:
        try:
            df_raw = load_skcm_from_r()
            df = prepare_pipeline_data(df_raw)
            CACHE_CSV.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(CACHE_CSV, index=False)
        except Exception as e:
            raise HTTPException(status_code=503, detail=f"Datos no disponibles: {e}")

    # Ejecutar pipeline
    try:
        result = run_fpgrowth_pipeline(
            df,
            min_support=params.min_support,
            min_confidence=params.min_confidence,
            min_lift=params.min_lift,
            min_leverage=params.min_leverage,
        )
    except Exception as e:
        logger.error(f"Error en pipeline: {e}")
        raise HTTPException(status_code=500, detail=str(e))

    # Marcar qué reglas calculadas coinciden con las del TFM
    for rule in result["rules_clinical"]:
        tfm_id = match_tfm_rule(rule)
        rule["tfm_match"] = tfm_id

    # Añadir las reglas del TFM siempre
    result["tfm_rules"] = TFM_RULES

    # Stats de coincidencia
    matched = sum(1 for r in result["rules_clinical"] if r.get("tfm_match"))
    result["tfm_match_stats"] = {
        "tfm_rules_total": len(TFM_RULES),
        "computed_rules_total": result["n_rules_clinical"],
        "matched": matched,
    }

    return JSONResponse(result)


@app.post("/api/upload-and-run")
async def upload_and_run(params: PipelineParams):
    """Ejecuta el pipeline sobre datos subidos manualmente."""
    pass  # Implementado en el frontend via /api/run-pipeline con datos en memoria


@app.get("/api/tfm-rules")
def get_tfm_rules():
    """Retorna las 11 reglas exactas del TFM (Tabla 3 de la memoria)."""
    return JSONResponse({"rules": TFM_RULES, "source": "TFM Tabla 3, Cell 2015"})


# ============================================================
# FRONTEND ESTÁTICO
# ============================================================

app.mount("/static", StaticFiles(directory="static"), name="static")


@app.get("/", response_class=HTMLResponse)
def serve_frontend():
    index = Path("static/index.html")
    if index.exists():
        return HTMLResponse(index.read_text())
    return HTMLResponse("<h1>Frontend no encontrado</h1>")
