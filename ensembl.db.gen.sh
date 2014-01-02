#! /bin/bash
TARGET_GENOME=${1}
DBLIST=$2
GENOME=`grep ${TARGET_GENOME} ${DBLIST} | cut -f2,2`
SPECIES=`grep ${TARGET_GENOME} ${DBLIST} | cut -f3,3`
SPECIES_NAME="$(tr '[:lower:]' '[:upper:]' <<< ${SPECIES:0:1})${SPECIES:1}"
echo ${SPECIES_NAME}
VERSION=`grep ${TARGET_GENOME} ${DBLIST} | cut -f4,4`
echo "wget -q -O ./tmp/${SPECIES_NAME}.${GENOME}.${VERSION}.gtf.gz ftp://ftp.ensembl.org/pub/release-${VERSION}/gtf/${SPECIES}/${SPECIES_NAME}.${GENOME}.${VERSION}.gtf.gz"
wget -q -O ./tmp/${SPECIES_NAME}.${GENOME}.${VERSION}.gtf.gz ftp://ftp.ensembl.org/pub/release-${VERSION}/gtf/${SPECIES}/${SPECIES_NAME}.${GENOME}.${VERSION}.gtf.gz
cd tmp
gunzip ${SPECIES_NAME}.${GENOME}.${VERSION}.gtf.gz
mv ${SPECIES_NAME}.${GENOME}.${VERSION}.gtf ${TARGET_GENOME}.ensembl.gtf
cd ..
echo "./gtf2txt_plot.pl ensembl ./tmp/${TARGET_GENOME}.ensembl.gtf ./tmp/${TARGET_GENOME}.ensembl.txt"
./gtf2txt_plot.pl ensembl ./tmp/${TARGET_GENOME}.ensembl.gtf ./tmp/${TARGET_GENOME}.ensembl.txt
echo "python filter_ensembl.py ./tmp/${TARGET_GENOME}.ensembl.gtf  ./tmp/${TARGET_GENOME}.ensembl.txt ./tmp/${TARGET_GENOME}.ensembl.biotype.txt ${TARGET_GENOME}"
python filter_ensembl.py ./tmp/${TARGET_GENOME}.ensembl.gtf  ./tmp/${TARGET_GENOME}.ensembl.txt ./tmp/${TARGET_GENOME}.ensembl.biotype.txt ${TARGET_GENOME}
echo "Rscript genDB.R ${TARGET_GENOME} ensembl ./tmp/${TARGET_GENOME}.ensembl.biotype.txt"
Rscript genDB.R ${TARGET_GENOME} ensembl ./tmp/${TARGET_GENOME}.ensembl.biotype.txt
