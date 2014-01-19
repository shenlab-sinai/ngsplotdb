#!/bin/bash

NPVer="3.00"

# jq is needed for parsing json file
# $ wget http://stedolan.github.io/jq/download/linux32/jq (32-bit system)
# $ wget http://stedolan.github.io/jq/download/linux64/jq (64-bit system)

## if jq not in the path, just download it
if [ ! -f "jq" ]; then
        wget http://stedolan.github.io/jq/download/linux64/jq
        chmod u+x jq
fi

if [ $# -lt 2 ]; then
	echo "Usage: $0 db_list.txt pipeline_type"
	exit 0
fi

CORES=4
DBLIST=$1
PIPELINE_TYPE=$2

# $GENOME, UCSC genome id from dblist based on the definition of pipeline json.
ID_COL=`cat json/"${PIPELINE_TYPE}".json | ./jq -r '.ID'`
GENOMES=( $(cut -f"${ID_COL}" ${DBLIST}) )
AM_COL=`cat json/"${PIPELINE_TYPE}".json | ./jq -r '.Assembly'`
AMS=( $(cut -f"${AM_COL}" ${DBLIST}) )
SP_COL=`cat json/"${PIPELINE_TYPE}".json | ./jq -r '.Species'`
SPS=( $(cut -f"${SP_COL}" ${DBLIST}) )
EN_COL=`cat json/"${PIPELINE_TYPE}".json | ./jq -r '.EnsVer'`
ENS=( $(cut -f"${EN_COL}" ${DBLIST}) )

ENS_PATH=`cat json/"${PIPELINE_TYPE}".json | ./jq -r '.EnsPath'`

No_DBLIST=${#GENOMES[@]}

mkdir -p database
mkdir -p tmp

function ngsplotdb(){
	i=$1
	NPVer=$2
	DBLIST=$3
	ENS_PATH=$4
	PIPELINE_TYPE=$5
	ID_COL=`cat json/"${PIPELINE_TYPE}".json | ./jq -r '.ID'`
	GENOMES=( $(cut -f"${ID_COL}" ${DBLIST}) )
	AM_COL=`cat json/"${PIPELINE_TYPE}".json | ./jq -r '.Assembly'`
	AMS=( $(cut -f"${AM_COL}" ${DBLIST}) )
	SP_COL=`cat json/"${PIPELINE_TYPE}".json | ./jq -r '.Species'`
	SPS=( $(cut -f"${SP_COL}" ${DBLIST}) )
	EN_COL=`cat json/"${PIPELINE_TYPE}".json | ./jq -r '.EnsVer'`
	ENS=( $(cut -f"${EN_COL}" ${DBLIST}) )

	echo -e "${GENOMES[${i}]}"
	echo "Processing ensembl annotation..."
	./ensembl.db.gen.sh ${GENOMES[${i}]} ${AMS[${i}]} ${SPS[${i}]} \
		${ENS[${i}]} ${ENS_PATH} ${PIPELINE_TYPE}
	echo "Processing refseq annotation..."
	./refseq.db.gen.sh ${GENOMES[${i}]}
	mkdir -p ./database/${GENOMES[${i}]}/
	sed -e "s/\*ID\*/${GENOMES[${i}]}/g" \
		-e "s/\*ASSEMBLY\*/${AMS[${i}]}/g" \
		-e "s/\*SPECIES\*/${SPS[${i}]}/g" \
		-e "s/\*ENSVER\*/${ENS[${i}]}/g" \
		-e "s/\*NPVER\*/${NPVer}/g" \
		metainfo_template.txt > ./database/${GENOMES[${i}]}/.metainfo
	cut -f1 ./tmp/${GENOMES[${i}]}.ensembl.biotype.txt|sort -u > \
		./database/${GENOMES[${i}]}/.chrnames.ensembl
	if [ -f ./tmp/${GENOMES[${i}]}.refseq.biotype.txt ]
		then
		cut -f1 ./tmp/${GENOMES[${i}]}.refseq.biotype.txt|sort -u > \
			./database/${GENOMES[${i}]}/.chrnames.refseq
	fi

	## process the data needed for region_analysis
	if [ -f ./tmp/${GENOMES[${i}]}.refseq.gp ]
		then
		echo "Processing annotation package for Region_Analysis..."
		echo -e "./gen_spec_anno.sh ${GENOMES[${i}]} ${NPVer} ${ENS[${i}]} ${SPS[${i}]} ${AMS[${i}]}"
		./gen_spec_anno.sh ${GENOMES[${i}]} ${NPVer} ${ENS[${i}]} ${SPS[${i}]} ${AMS[${i}]}
	fi

	# add CpG islands annotations
	# "select * " will interrupt the output of MySQL, because such as panTro4 and rheMac2 will
	# return different colunms.
	echo -e "mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N \
	-e select chrom, chromStart as start, chromEnd as end, name from cpgIslandExt;" ${GENOMES[${i}]} \
	"> ./tmp/${GENOMES[${i}]}.cgi.gp"
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N \
	-e "select chrom, chromStart as start, chromEnd as end, name from cpgIslandExt;" ${GENOMES[${i}]} \
	> ./tmp/${GENOMES[${i}]}.cgi.gp
	if [ -s "./tmp/${GENOMES[${i}]}.cgi.gp" ]
		then
		echo "Processing CGI information..."
		cut -f1,2,3 ./tmp/${GENOMES[${i}]}.cgi.gp > ./tmp/${GENOMES[${i}]}.cgi.bed
		echo -e "region_analysis.py -d ensembl -g ${GENOMES[${i}]} -i ./tmp/${GENOMES[${i}]}.cgi.bed"
		region_analysis.py -d ensembl -g ${GENOMES[${i}]} -i ./tmp/${GENOMES[${i}]}.cgi.bed
		mv ./tmp/${GENOMES[${i}]}.cgi.bed.annotated ./tmp/${GENOMES[${i}]}.cgi.ensembl.txt
		echo -e "region_analysis.py -d refseq -g ${GENOMES[${i}]} -i ./tmp/${GENOMES[${i}]}.cgi.bed"
		region_analysis.py -d refseq -g ${GENOMES[${i}]} -i ./tmp/${GENOMES[${i}]}.cgi.bed
		mv ./tmp/${GENOMES[${i}]}.cgi.bed.annotated ./tmp/${GENOMES[${i}]}.cgi.refseq.txt

		echo -e "genDB_region_analysis_feature.R ${GENOMES[${i}]} cgi ensembl ./tmp/${GENOMES[${i}]}.ensembl.biotype.txt"
		Rscript genDB_region_analysis_feature.R ${GENOMES[${i}]} cgi ensembl ./tmp/${GENOMES[${i}]}.ensembl.txt
		echo -e "genDB_region_analysis_feature.R ${GENOMES[${i}]} cgi refseq ./tmp/${GENOMES[${i}]}.refseq.biotype.txt"
		Rscript genDB_region_analysis_feature.R ${GENOMES[${i}]} cgi refseq ./tmp/${GENOMES[${i}]}.refseq.txt
	else
		echo "No CGI annotation of ${GENOMES[${i}]} genome in UCSC! Skip the CGI annotation step!"
	fi

	mv ${GENOMES[${i}]}*.RData ./database/${GENOMES[${i}]}/
	cd ./database/
	tar czvf ngsplotdb_${GENOMES[${i}]}_${ENS[${i}]}_${NPVer}.tar.gz ${GENOMES[${i}]}
	cd ..
	return 0
}
export -f ngsplotdb

seq 0 `echo ${No_DBLIST}-1 | bc` | xargs -n 1 -I {} -P ${CORES} bash -c \
	"ngsplotdb {} ${NPVer} ${DBLIST} ${ENS_PATH} ${PIPELINE_TYPE}"
