# Base: R + Python en Ubuntu
FROM rocker/r-ver:4.3.3

# Instalar Python y dependencias del sistema
RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Instalar BiocManager y TCGAbiolinks en R
RUN R -e "install.packages('BiocManager', repos='https://cran.r-project.org')"
RUN R -e "BiocManager::install('TCGAbiolinks', ask=FALSE, update=FALSE)"

# Crear directorios
RUN mkdir -p /data_cache /scripts

# Script R que exporta los datos SKCM a CSV
RUN echo 'library(TCGAbiolinks)\n\
skcm <- TCGAquery_subtype(tumor = "skcm")\n\
cols <- c("MUTATIONSUBTYPES","UV-signature","RNASEQ-CLUSTER_CONSENHIER",\n\
          "MethTypes.201408","MIRCluster","LYMPHOCYTE.SCORE")\n\
available <- cols[cols %in% colnames(skcm)]\n\
out <- skcm[, available, drop=FALSE]\n\
write.csv(out, "/data_cache/skcm_subtipos.csv", row.names=FALSE)\n\
cat("Exportado:", nrow(out), "pacientes\\n")\n\
' > /scripts/export_skcm.R

# Pre-cachear datos SKCM como CSV en build time
RUN Rscript /scripts/export_skcm.R || echo "Cache fallido, se hara en runtime"

WORKDIR /app

# Instalar dependencias Python (sin rpy2)
COPY web_app/requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

# Copiar la app
COPY web_app/main.py .
COPY web_app/static/ ./static/

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
