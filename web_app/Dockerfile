# Usamos ubuntu con R ya instalado via rocker, pero añadimos rpy2 del sistema
FROM rocker/r-ver:4.3.3

# Instalar Python, dependencias y rpy2 del sistema (pre-compilado, sin necesidad de compilar)
RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-dev \
    python3-rpy2 \
    gcc g++ make \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# Instalar BiocManager y TCGAbiolinks en R
RUN R -e "install.packages('BiocManager', repos='https://cran.r-project.org')"
RUN R -e "BiocManager::install('TCGAbiolinks', ask=FALSE, update=FALSE)"

# Crear directorio de caché
RUN mkdir -p /data_cache /repo_data

# Pre-cachear datos SKCM en build time
RUN R -e "\
  tryCatch({\
    library(TCGAbiolinks);\
    skcm <- TCGAquery_subtype('skcm');\
    saveRDS(skcm, '/data_cache/skcm_subtipos.rds');\
    message('Cache OK: ', nrow(skcm), ' pacientes');\
  }, error=function(e) message('Cache build failed, se hara en runtime: ', e\$message))\
"

WORKDIR /app

ENV R_HOME=/usr/local/lib/R

# Instalar dependencias Python (sin rpy2 — viene del sistema)
COPY web_app/requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

# Copiar la app
COPY web_app/main.py .
COPY web_app/static/ ./static/

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
