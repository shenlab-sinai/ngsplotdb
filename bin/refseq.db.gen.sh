#! /bin/bash

PRO_BIN_PATH=$2
NP_DB_PATH=$3

echo -e "mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e 'select * from refGene;' ${1} \
	> ${NP_DB_PATH}/tmp/${1}.refseq.gp"
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from refGene;" ${1} \
	> ${NP_DB_PATH}/tmp/${1}.refseq.gp
if [ -s "${NP_DB_PATH}/tmp/${1}.refseq.gp" ]; then
	cut -f 2- ${NP_DB_PATH}/tmp/${1}.refseq.gp | \
		${PRO_BIN_PATH}/genePredToGtf file stdin ${NP_DB_PATH}/tmp/${1}.refseq.gtf
	echo -e "perl -I ${PRO_BIN_PATH} ${PRO_BIN_PATH}/gtf2txt_plot.pl refseq \
		${NP_DB_PATH}/tmp/${1}.refseq.gtf ${NP_DB_PATH}/tmp/${1}.refseq.txt"
	perl -I ${PRO_BIN_PATH} ${PRO_BIN_PATH}/gtf2txt_plot.pl refseq \
		${NP_DB_PATH}/tmp/${1}.refseq.gtf ${NP_DB_PATH}/tmp/${1}.refseq.txt
	echo -e "${PRO_BIN_PATH}/filter_refseq.py ${NP_DB_PATH}/tmp/${1}.refseq.txt \
		${NP_DB_PATH}/tmp/${1}.refseq.biotype.txt ${1}"
	${PRO_BIN_PATH}/filter_refseq.py ${NP_DB_PATH}/tmp/${1}.refseq.txt \
		${NP_DB_PATH}/tmp/${1}.refseq.biotype.txt ${1}
	echo -e "Rscript ${PRO_BIN_PATH}/genDB.R ${1} refseq \
		${NP_DB_PATH}/tmp/${1}.refseq.biotype.txt"
	Rscript ${PRO_BIN_PATH}/genDB.R ${1} refseq \
		${NP_DB_PATH}/tmp/${1}.refseq.biotype.txt
else
	echo -e "No ${1} genome annotation in UCSC!"
	rm ${NP_DB_PATH}/tmp/${1}.refseq.gp
fi
