#!/usr/bin/env bash
# Download supplementary file from GEO (GSE236754, ~7.8 MB)
set -euo pipefail
mkdir -p data/raw
cd data/raw
curl -L -o GSE236754_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE236754&format=file"
mkdir -p GSE236754_RAW && tar -xf GSE236754_RAW.tar -C GSE236754_RAW
ls -lh GSE236754_RAW