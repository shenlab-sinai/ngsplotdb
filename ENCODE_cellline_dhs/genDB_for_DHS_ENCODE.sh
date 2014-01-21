#! /bin/bash

NPVer="3.00"
GENOME="hg19"
FEATURE="dhs"

CORES=4

mkdir -p data

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseUniform/files.txt

cut -f1,1 files.txt > list.txt
for BED in `cat list.txt`; do
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseUniform/${BED}
done
mv *.gz data/
cd data
gunzip *.gz
cd ..

mkdir -p ./tmp
cp ${GENOME}.*.biotype.txt ./tmp
cd data
for i in *.narrowPeak; do
	echo $i
	cut -f 1,2,3 ${i} > ../tmp/${i/narrowPeak/bed}
done
cd ..

cd tmp
for i in wgEncodeAwgDnaseUwduke*.bed; do
	mv ${i} ${i/wgEncodeAwgDnaseUwduke/}
done

for i in wgEncodeAwgDnaseDuke*.bed; do
	mv ${i} ${i/wgEncodeAwgDnaseDuke/}
done

for i in wgEncodeAwgDnaseUw*.bed; do
	mv ${i} ${i/wgEncodeAwgDnaseUw/}
done

for i in *UniPk.bed; do
	mv ${i} ${i/UniPk/}
done

CELLLINES=`ls *.bed`
CELLLINES=${CELLLINES//\.bed/}

function Bed2DB(){
	CELLLINE=$1
	GENOME=$2
	FEATURE=$3

	echo -e "region_analysis.py -d ensembl -g ${GENOME} -i ${CELLLINE}.bed"
	region_analysis.py -d ensembl -g ${GENOME} -i ${CELLLINE}.bed
	mv ${CELLLINE}.bed.annotated ${GENOME}.${FEATURE}.${CELLLINE}.ensembl.txt
	echo -e "region_analysis.py -d refseq -g ${GENOME} -i ${CELLLINE}.bed"
	region_analysis.py -d refseq -g ${GENOME} -i ${CELLLINE}.bed
	mv ${CELLLINE}.bed.annotated ${GENOME}.${FEATURE}.${CELLLINE}.refseq.txt

	echo -e "genDB_region_analysis_feature.R ${GENOME} ${FEATURE}.${CELLLINE} ensembl ${GENOME}.ensembl.txt"
	Rscript ../genDB_region_analysis_feature.R ${GENOME} ${FEATURE}.${CELLLINE} ensembl ${GENOME}.ensembl.biotype.txt
	echo -e "genDB_region_analysis_feature.R ${GENOME} ${FEATURE}.${CELLLINE} refseq ${GENOME}.refseq.txt"
	Rscript ../genDB_region_analysis_feature.R ${GENOME} ${FEATURE}.${CELLLINE} refseq ${GENOME}.refseq.biotype.txt
}
export -f Bed2DB

echo -e "${CELLLINES[@]}" | xargs -n 1 -I {} -P $CORES bash -c "Bed2DB {} ${GENOME} ${FEATURE}"


# mkdir -p database
# mv *.RData ./database