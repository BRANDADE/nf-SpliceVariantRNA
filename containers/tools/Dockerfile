# Image des scripts de post-traitement (bin/sl_filter.py, sl_recap.py, sl_sashimi.py) et de ggsashimi.
#
# Construction (machine avec accès à Internet), puis conversion en .sif pour le cluster :
#   docker build -t splicevariant-tools:1.0 containers/tools
#   singularity build splicevariant-tools_1.0.sif docker-daemon://splicevariant-tools:1.0
# et renseigner --tools_container /chemin/absolu/splicevariant-tools_1.0.sif
#
# rocker/r-ver fige les paquets CRAN à la date de sortie de la version de R (dépôt Posit PPM).
FROM rocker/r-ver:4.4.2

LABEL description="nf-SpliceVariantRNA : post-traitement SpliceLauncher (Python) + ggsashimi 1.1.5" \
      version="1.0"

# ggsashimi : R (ggplot2, data.table, gridExtra, svglite) + Python (pysam)
RUN install2.r --error --skipinstalled ggplot2 data.table gridExtra svglite \
    && rm -rf /tmp/downloaded_packages

RUN apt-get update \
    && apt-get install -y --no-install-recommends python3 python3-venv ca-certificates curl procps \
    && rm -rf /var/lib/apt/lists/*

# Versions testées (tests/python, CI)
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir \
        pandas==3.0.6 pysam==0.24.1 reportlab==5.0.1 pypdf==6.19.0 pdfplumber==0.11.10
ENV PATH="/opt/venv/bin:${PATH}"

ARG GGSASHIMI_VERSION=v1.1.5
RUN curl -fsSL -o /usr/local/bin/ggsashimi.py \
        https://raw.githubusercontent.com/guigolab/ggsashimi/${GGSASHIMI_VERSION}/ggsashimi.py \
    && sed -i '1s|.*|#!/opt/venv/bin/python3|' /usr/local/bin/ggsashimi.py \
    && chmod +x /usr/local/bin/ggsashimi.py
