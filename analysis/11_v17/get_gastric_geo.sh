#!/usr/bin/env bash
# Download additional public gastric cancer cohorts for a locked external validation of
# NetImmune / CAF_ECM (and an independent single-cell atlas for the attribution check).
# Run on a node WITH internet access (the compute node used for the analyses has none):
#   bash /groups/adv2105_gp/yichen/Yi/com/gastric/2/data_download/get_gastric_geo.sh
# Series matrices and small supplementary files only; large RAW tarballs are skipped.
set -euo pipefail
OUT=/groups/adv2105_gp/yichen/Yi/com/gastric/2/data_download/geo
mkdir -p "$OUT"; cd "$OUT"

# bulk cohorts with survival (series matrix + non-RAW supplementary files, which carry
# outcome tables for some series). GSE34942/GSE35809 are Singapore series that may share
# patients with GSE15459 -- de-duplicate by sample before pooling. GSE51105 (Busuttil
# 2014) is the stromal/macrophage prognosis cohort; GSE14210 has chemotherapy response.
BULK="GSE15459 GSE26253 GSE26901 GSE13861 GSE29272 GSE34942 GSE57303 GSE35809 GSE51105 GSE22377 GSE14210"
# oesophageal adenocarcinoma checkpoint cohort (PERFECT, van den Ende 2021): the compendium copy kept
# only 3 of the 48 module genes, so take the deposited matrix
BULK="$BULK GSE165252"
# independent gastric single-cell atlases (supplementary count matrices can be large)
SC="GSE163558 GSE134520"

stub() { local g=$1; echo "${g:0:${#g}-3}nnn"; }

for g in $BULK; do
  d=$(stub "$g"); mkdir -p "$g"
  echo "== $g"
  wget -q -nc -e netrc=off -P "$g" "https://ftp.ncbi.nlm.nih.gov/geo/series/$d/$g/matrix/" -r -np -nd -A "*series_matrix.txt.gz" -e robots=off || echo "  matrix: failed"
  wget -q -nc -e netrc=off -P "$g" "https://ftp.ncbi.nlm.nih.gov/geo/series/$d/$g/suppl/" -r -np -nd -R "*RAW.tar,index.html*" -e robots=off || echo "  suppl: none or failed"
done

# single-cell supplementary files are large; run with DO_SC=1 to fetch them
[ "${DO_SC:-0}" = "1" ] || SC=""
for g in $SC; do
  d=$(stub "$g"); mkdir -p "$g"
  echo "== $g (single-cell)"
  wget -q -nc -e netrc=off -P "$g" "https://ftp.ncbi.nlm.nih.gov/geo/series/$d/$g/matrix/" -r -np -nd -A "*series_matrix.txt.gz" -e robots=off || true
  wget -q -nc -e netrc=off -P "$g" "https://ftp.ncbi.nlm.nih.gov/geo/series/$d/$g/suppl/" -r -np -nd -R "index.html*" -e robots=off || true
done

# METABRIC (cBioPortal study bundle) for breast prognosis paired with the I-SPY2 trial
# analysis; the earlier S3 link in public_imc/ returned AccessDenied. Tries the current
# datahub host first, then the old S3 bucket.
mkdir -p METABRIC
( cd METABRIC && { wget -q -nc -e netrc=off https://datahub.assets.cbioportal.org/brca_metabric.tar.gz \
    || wget -q -nc -e netrc=off https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz \
    || echo "  METABRIC: failed (download manually from cbioportal.org/study/summary?id=brca_metabric)"; } )

echo; du -sh "$OUT"/* 2>/dev/null
echo "Done. Tell Claude the files are in $OUT."
