#!/bin/bash

if [ $# -lt 1 ]; then
	echo "Usage: $0 db_list.txt"
	exit 0
fi

DBLIST=$1

# $GENOME, UCSC genome id
GENOMES=`cut -f1,1 ${DBLIST}`
mkdir -p database
mkdir -p tmp

for GENOME in ${GENOMES}; do

	## database processing
	./ensembl.db.gen.sh ${GENOME} ${DBLIST}
	./refseq.db.gen.sh ${GENOME}

	mkdir -p ./database/${GENOME}/
	## mv ${GENOME}*.RData ./database/${GENOME}/

	grep ${GENOME} ${DBLIST}|tr "\t" "\n" > vals.tmp
	# echo ${NPVER} >> vals.tmp
	echo -e "ID\nAssembly\nSpecies\nEnsVer\nNPVer" > vars.tmp
	paste vars.tmp vals.tmp > ./database/${GENOME}/.metainfo
	rm -f vals.tmp vars.tmp

	cut -f1 ./tmp/${GENOME}.ensembl.biotype.txt|sort -u > ./database/${GENOME}/.chrnames.ensembl
	cut -f1 ./tmp/${GENOME}.refseq.biotype.txt|sort -u > ./database/${GENOME}/.chrnames.refseq

	ENSVER=`grep ${GENOME} ${DBLIST}|cut -f4`
	NPVER=`grep ${GENOME} ${DBLIST}|cut -f5`

	## process the data needed for region_analysis
	./gen_spec_anno.sh ${GENOME}

	## add CpG islands annotations
	## "select * " will interrupt the output of MySQL, because such as panTro4 and rheMac2 will
	## return different colunms.
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N \
	-e "select chrom, chromStart as start, chromEnd as end, name from cpgIslandExt;" ${GENOME} \
	> ./tmp/${GENOME}.cgi.gp
	if [ -s "./tmp/${GENOME}.cgi.gp" ]; then
		cut -f1,2,3 ./tmp/${GENOME}.cgi.gp > ./tmp/${GENOME}.cgi.bed
		region_analysis.py -d ensembl -g ${GENOME} -i ./tmp/${GENOME}.cgi.bed
		mv ./tmp/${GENOME}.cgi.bed.annotated ./tmp/${GENOME}.cgi.ensembl.txt
		region_analysis.py -d refseq -g ${GENOME} -i ./tmp/${GENOME}.cgi.bed
		mv ./tmp/${GENOME}.cgi.bed.annotated ./tmp/${GENOME}.cgi.refseq.txt

		Rscript genDB_region_analysis_feature.R ${GENOME} cgi ensembl ./tmp/${GENOME}.ensembl.biotype.txt
		Rscript genDB_region_analysis_feature.R ${GENOME} cgi refseq ./tmp/${GENOME}.refseq.biotype.txt
	fi

	mv ${GENOME}*.RData ./database/${GENOME}/
	cd ./database/
	tar czvf ngsplotdb_${GENOME}_${ENSVER}_${NPVER}.tar.gz ${GENOME}
	cd ..
done

