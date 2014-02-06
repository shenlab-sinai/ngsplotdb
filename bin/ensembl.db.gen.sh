#! /bin/bash
TARGET_GENOME=$1
GENOME=$2
SPECIES=$3
SPECIES_NAME="$(tr '[:lower:]' '[:upper:]' <<< ${SPECIES:0:1})${SPECIES:1}"
echo ${SPECIES_NAME}
VERSION=$4
ENS_PATH=$5
PIPELINE_JSON=$6
PIPELINE_TYPE=$( basename ${PIPELINE_JSON} .json )
i=$7
DBLIST=$8
PRO_BIN_PATH=$9
NP_DB_PATH=${10}


case ${PIPELINE_TYPE} in 
	"animal" )
		GFF=`echo ${ENS_PATH} | sed -e "s/\*SPECIES_NAME\*/${SPECIES_NAME}/g" \
			-e "s/\*SPECIES\*/${SPECIES}/g" \
			-e "s/\*VERSION\*/${VERSION}/g" \
			-e "s/\*GENOME\*/${GENOME}/g"`
			;;
	"plant" )
		GFF=`echo ${ENS_PATH} | sed -e "s/\*SPECIES_NAME\*/${SPECIES_NAME}/g" \
			-e "s/\*SPECIES\*/${SPECIES}/g" \
			-e "s/\*VERSION\*/${VERSION}/g" \
			-e "s/\*GENOME\*/${GENOME}/g"`
			;;
	"bacteria" )
		BAC_COLLECTION_COL=`cat ${PIPELINE_JSON} | ${PRO_BIN_PATH}/jq -r '.Bacteria_Collection'`
		BAC_COLLECTS=( $(cut -f"${BAC_COLLECTION_COL}" ${DBLIST}) )
		BAC_COLLECT=${BAC_COLLECTS[${i}]}
		GFF=`echo ${ENS_PATH} | sed -e "s/\*SPECIES_NAME\*/${SPECIES_NAME}/g" \
			-e "s/\*BACTERIA_COLLECTION\*/${BAC_COLLECT}/g" \
			-e "s/\*SPECIES\*/${SPECIES}/g" \
			-e "s/\*VERSION\*/${VERSION}/g" \
			-e "s/\*GENOME\*/${GENOME}/g"`
			;;
	*) echo "Unknown PIPELINE_TYPE!"
esac

echo -e "wget -q -O ${NP_DB_PATH}/tmp/${SPECIES_NAME}.${GENOME}.${VERSION}.gtf.gz ${GFF}"
wget -q -O ${NP_DB_PATH}/tmp/${SPECIES_NAME}.${GENOME}.${VERSION}.gtf.gz ${GFF}
gzip -dc < ${NP_DB_PATH}/tmp/${SPECIES_NAME}.${GENOME}.${VERSION}.gtf.gz \
	> ${NP_DB_PATH}/tmp/${TARGET_GENOME}.ensembl.gtf
echo -e "perl -I ${PRO_BIN_PATH} ${PRO_BIN_PATH}/gtf2txt_plot.pl ensembl \
	${NP_DB_PATH}/tmp/${TARGET_GENOME}.ensembl.gtf \
	${NP_DB_PATH}/tmp/${TARGET_GENOME}.ensembl.txt"
perl -I ${PRO_BIN_PATH} ${PRO_BIN_PATH}/gtf2txt_plot.pl ensembl \
	${NP_DB_PATH}/tmp/${TARGET_GENOME}.ensembl.gtf \
	${NP_DB_PATH}/tmp/${TARGET_GENOME}.ensembl.txt
echo -e "${PRO_BIN_PATH}/filter_ensembl.py ${NP_DB_PATH}/tmp/${TARGET_GENOME}.ensembl.gtf \
	${NP_DB_PATH}/tmp/${TARGET_GENOME}.ensembl.txt \
	${NP_DB_PATH}/tmp/${TARGET_GENOME}.ensembl.biotype.txt ${TARGET_GENOME}"
${PRO_BIN_PATH}/filter_ensembl.py ${NP_DB_PATH}/tmp/${TARGET_GENOME}.ensembl.gtf  \
	${NP_DB_PATH}/tmp/${TARGET_GENOME}.ensembl.txt \
	${NP_DB_PATH}/tmp/${TARGET_GENOME}.ensembl.biotype.txt ${TARGET_GENOME}
echo -e "Rscript ${PRO_BIN_PATH}/genDB.R ${TARGET_GENOME} ensembl \
	${NP_DB_PATH}/tmp/${TARGET_GENOME}.ensembl.biotype.txt"
Rscript ${PRO_BIN_PATH}/genDB.R ${TARGET_GENOME} ensembl \
	${NP_DB_PATH}/tmp/${TARGET_GENOME}.ensembl.biotype.txt
