#! /bin/bash
TARGET_GENOME=$1
GENOME=$2
SPECIES=$3
SPECIES_NAME="$(tr '[:lower:]' '[:upper:]' <<< ${SPECIES:0:1})${SPECIES:1}"
echo ${SPECIES_NAME}
VERSION=$4
ENS_PATH=$5
PIPELINE_TYPE=$6
i=$7
DBLIST=$8

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
		BAC_COLLECTION_COL=`cat json/"${PIPELINE_TYPE}".json | ./jq -r '.Bacteria_Collection'`
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

echo -e "wget -q -O ./tmp/${SPECIES_NAME}.${GENOME}.${VERSION}.gtf.gz ${GFF}"
wget -q -O ./tmp/${SPECIES_NAME}.${GENOME}.${VERSION}.gtf.gz ${GFF}
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
