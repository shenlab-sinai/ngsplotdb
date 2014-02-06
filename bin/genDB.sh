#!/bin/bash

NPVer="3.00"

# relative path of script
export PRO_BIN_PATH="`dirname \"$0\"`"
export NP_DB_PATH="./output"

# jq is needed for parsing json file
# $ wget http://stedolan.github.io/jq/download/linux32/jq (32-bit system)
# $ wget http://stedolan.github.io/jq/download/linux64/jq (64-bit system)

## if jq not in the path, just download it
if [ ! -f "${PRO_BIN_PATH}/jq" ]; then
        wget -q -O ${PRO_BIN_PATH}/jq http://stedolan.github.io/jq/download/linux64/jq
        chmod u+x ${PRO_BIN_PATH}/jq
fi

if [ $# -lt 2 ]; then
	echo "Usage: $0 db_list.txt pipeline_type.json"
	exit 0
fi

CORES=4
DBLIST=$1
PIPELINE_JSON=$2

# $GENOME, UCSC genome id from dblist based on the definition of pipeline json.
ID_COL=`cat ${PIPELINE_JSON} | ${PRO_BIN_PATH}/jq -r '.ID'`
GENOMES=( $(cut -f"${ID_COL}" ${DBLIST}) )
AM_COL=`cat ${PIPELINE_JSON} | ${PRO_BIN_PATH}/jq -r '.Assembly'`
AMS=( $(cut -f"${AM_COL}" ${DBLIST}) )
SP_COL=`cat ${PIPELINE_JSON} | ${PRO_BIN_PATH}/jq -r '.Species'`
SPS=( $(cut -f"${SP_COL}" ${DBLIST}) )
EN_COL=`cat ${PIPELINE_JSON} | ${PRO_BIN_PATH}/jq -r '.EnsVer'`
ENS=( $(cut -f"${EN_COL}" ${DBLIST}) )

ENS_PATH=`cat ${PIPELINE_JSON} | ${PRO_BIN_PATH}/jq -r '.EnsPath'`

No_DBLIST=${#GENOMES[@]}

mkdir -p ${NP_DB_PATH}/database
mkdir -p ${NP_DB_PATH}/tmp

function ngsplotdb(){
	i=$1
	NPVer=$2
	DBLIST=$3
	ENS_PATH=$4
	PIPELINE_JSON=$5
	PRO_BIN_PATH=$6
	## To get data folder under project, Mac OSX doesn't support readlink
	## PRO_DATA_PATH=$(readlink -m ${PRO_BIN_PATH}/../data)
	PRO_DATA_PATH=$(cd ${PRO_BIN_PATH}/../data; pwd)
	NP_DB_PATH=$7

	ID_COL=`cat ${PIPELINE_JSON} | ${PRO_BIN_PATH}/jq -r '.ID'`
	GENOMES=( $(cut -f"${ID_COL}" ${DBLIST}) )
	AM_COL=`cat ${PIPELINE_JSON} | ${PRO_BIN_PATH}/jq -r '.Assembly'`
	AMS=( $(cut -f"${AM_COL}" ${DBLIST}) )
	SP_COL=`cat ${PIPELINE_JSON} | ${PRO_BIN_PATH}/jq -r '.Species'`
	SPS=( $(cut -f"${SP_COL}" ${DBLIST}) )
	EN_COL=`cat ${PIPELINE_JSON} | ${PRO_BIN_PATH}/jq -r '.EnsVer'`
	ENS=( $(cut -f"${EN_COL}" ${DBLIST}) )

	echo -e "${GENOMES[${i}]}"
	echo "Processing ensembl annotation..."
	echo -e "${PRO_BIN_PATH}/ensembl.db.gen.sh ${GENOMES[${i}]} ${AMS[${i}]} ${SPS[${i}]} \
		${ENS[${i}]} ${ENS_PATH} ${PIPELINE_JSON} ${i} ${DBLIST} ${PRO_BIN_PATH} ${NP_DB_PATH}"
	${PRO_BIN_PATH}/ensembl.db.gen.sh ${GENOMES[${i}]} ${AMS[${i}]} ${SPS[${i}]} \
		${ENS[${i}]} ${ENS_PATH} ${PIPELINE_JSON} ${i} ${DBLIST} ${PRO_BIN_PATH} ${NP_DB_PATH}
	echo "Processing refseq annotation..."
	${PRO_BIN_PATH}/refseq.db.gen.sh ${GENOMES[${i}]} ${PRO_BIN_PATH} ${NP_DB_PATH}
	mkdir -p ${NP_DB_PATH}/database/${GENOMES[${i}]}/
	sed -e "s/\*ID\*/${GENOMES[${i}]}/g" \
		-e "s/\*ASSEMBLY\*/${AMS[${i}]}/g" \
		-e "s/\*SPECIES\*/${SPS[${i}]}/g" \
		-e "s/\*ENSVER\*/${ENS[${i}]}/g" \
		-e "s/\*NPVER\*/${NPVer}/g" \
		${PRO_DATA_PATH}/metainfo_template.txt > ${NP_DB_PATH}/database/${GENOMES[${i}]}/.metainfo
	cut -f1 ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.ensembl.biotype.txt|sort -u > \
		${NP_DB_PATH}/database/${GENOMES[${i}]}/.chrnames.ensembl
	if [ -f ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.refseq.biotype.txt ]
		then
		cut -f1 ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.refseq.biotype.txt|sort -u > \
			${NP_DB_PATH}/database/${GENOMES[${i}]}/.chrnames.refseq
	fi

	## process the data needed for region_analysis
	if [ -f ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.refseq.gp ]
		then
		echo "Processing annotation package for Region_Analysis..."
		echo -e "${PRO_BIN_PATH}/gen_spec_anno.sh ${GENOMES[${i}]} ${NPVer} \
			${ENS[${i}]} ${SPS[${i}]} ${AMS[${i}]} ${PRO_BIN_PATH} ${NP_DB_PATH}"
		${PRO_BIN_PATH}/gen_spec_anno.sh ${GENOMES[${i}]} ${NPVer} \
			${ENS[${i}]} ${SPS[${i}]} ${AMS[${i}]} ${PRO_BIN_PATH} ${NP_DB_PATH}
	fi

	# add CpG islands annotations
	# "select * " will interrupt the output of MySQL, because such as panTro4 and rheMac2 will
	# return different colunms.
	echo -e "mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N \
	-e select chrom, chromStart as start, chromEnd as end, name from cpgIslandExt;" ${GENOMES[${i}]} \
	"> ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.gp"
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N \
	-e "select chrom, chromStart as start, chromEnd as end, name from cpgIslandExt;" ${GENOMES[${i}]} \
	> ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.gp
	if [ -s "${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.gp" ]
		then
		echo "Processing CGI information..."
		cut -f1,2,3 ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.gp > ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.bed
		echo -e "region_analysis.py -d ensembl -g ${GENOMES[${i}]} -i ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.bed"
		region_analysis.py -d ensembl -g ${GENOMES[${i}]} -i ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.bed
		mv ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.bed.annotated ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.ensembl.txt
		echo -e "region_analysis.py -d refseq -g ${GENOMES[${i}]} -i ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.bed"
		region_analysis.py -d refseq -g ${GENOMES[${i}]} -i ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.bed
		mv ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.bed.annotated ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.refseq.txt

		echo -e "${PRO_BIN_PATH}/genDB_region_analysis_feature.R ${GENOMES[${i}]} cgi \
			ensembl ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.ensembl.txt \
			${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.ensembl.txt"
		Rscript ${PRO_BIN_PATH}/genDB_region_analysis_feature.R ${GENOMES[${i}]} cgi \
			ensembl ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.ensembl.txt \
			${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.ensembl.txt
		echo -e "genDB_region_analysis_feature.R ${GENOMES[${i}]} cgi \
			refseq ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.refseq.txt \
			${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.refseq.txt"
		Rscript ${PRO_BIN_PATH}/genDB_region_analysis_feature.R ${GENOMES[${i}]} cgi \
			refseq ${NP_DB_PATH}/tmp/${GENOMES[${i}]}.refseq.txt \
			${NP_DB_PATH}/tmp/${GENOMES[${i}]}.cgi.refseq.txt
	else
		echo "No CGI annotation of ${GENOMES[${i}]} genome in UCSC! Skip the CGI annotation step!"
	fi

	mv ${GENOMES[${i}]}*.RData ${NP_DB_PATH}/database/${GENOMES[${i}]}/
	cd ${NP_DB_PATH}/database/
	tar czvf ngsplotdb_${GENOMES[${i}]}_${ENS[${i}]}_${NPVer}.tar.gz ${GENOMES[${i}]}
	return 0
}
export -f ngsplotdb

seq 0 `echo ${No_DBLIST}-1 | bc` | xargs -n 1 -I {} -P ${CORES} bash -c \
	"ngsplotdb {} ${NPVer} ${DBLIST} ${ENS_PATH} ${PIPELINE_JSON} ${PRO_BIN_PATH} ${NP_DB_PATH}"
