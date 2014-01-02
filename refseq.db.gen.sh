#! /bin/bash

echo "mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e 'select * from refGene;' ${1} > ./tmp/${1}.refseq.gp"
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from refGene;" ${1} > ./tmp/${1}.refseq.gp
if [ -s "./tmp/${1}.refseq.gp" ]; then
	cut -f 2- ./tmp/${1}.refseq.gp | ./genePredToGtf file stdin ./tmp/${1}.refseq.gtf
	echo "./gtf2txt_plot.pl refseq ./tmp/${1}.refseq.gtf ./tmp/${1}.refseq.txt"
	./gtf2txt_plot.pl refseq ./tmp/${1}.refseq.gtf ./tmp/${1}.refseq.txt
	echo "python filter_refseq.py ./tmp/${1}.refseq.txt ./tmp/${1}.refseq.biotype.txt ${1}"
	python filter_refseq.py ./tmp/${1}.refseq.txt ./tmp/${1}.refseq.biotype.txt ${1}
	echo "Rscript genDB.R ${1} refseq ./tmp/${1}.refseq.biotype.txt"
	Rscript genDB.R ${1} refseq ./tmp/${1}.refseq.biotype.txt
else
	echo "No ${1} genome annotation in UCSC!"
	rm ./tmp/${1}.refseq.gp
fi

