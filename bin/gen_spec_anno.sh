#! /bin/bash

GENOME=${1}
NPVER=${2}
ENSVER=${3}
SPECIES=${4}
ASSEMBLY=${5}
PRO_BIN_PATH=${6}
## To get data folder under project, Mac OSX doesn't support readlink
## PRO_DATA_PATH=$(readlink -m ${PRO_BIN_PATH}/../data)
PRO_DATA_PATH=$(cd ${PRO_BIN_PATH}/../data; pwd)
NP_DB_PATH=${7}

GENOME_ANNO_DIR="${NP_DB_PATH}/tmp"
RA_PATH="${NP_DB_PATH}/RA/${GENOME}_${NPVER}_En${ENSVER}"

mkdir -p ${RA_PATH}

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
	"select chrom, size from "${GENOME}".chromInfo"  > ${RA_PATH}/${GENOME}.genome

## if no ready info, try to fetch centromere and telomere info from UCSC
if [ ! -d "${PRO_DATA_PATH}/spec_anno/${GENOME}" ]; then
	echo "Querying centromere and telomere information from UCSC..."
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D ${GENOME} -e 'select chrom, chromStart, chromEnd from gap where type="centromere"' \
	| grep -v chromStart > ${RA_PATH}/${GENOME}_centromere.bed
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D ${GENOME} -e 'select chrom, chromStart, chromEnd from gap where type="telomere"' \
	| grep -v chromStart > ${RA_PATH}/${GENOME}_telomere.bed
	if [ ! -s "${RA_PATH}/${GENOME}_centromere.bed" ]; then
		echo "No centromere info in ${GENOME} database, dumb entry will be used!"
		rm ${RA_PATH}/${GENOME}_centromere.bed
		cp -f ${PRO_DATA_PATH}/spec_anno/dumb/dumb_pericentromere.bed ${RA_PATH}/${GENOME}_centromere.bed
		cp -f ${PRO_DATA_PATH}/spec_anno/dumb/dumb_pericentromere.bed ${RA_PATH}/${GENOME}_pericentromere.bed
	fi
	if [ ! -s "${RA_PATH}/${GENOME}_telomere.bed" ]; then
		echo "No telomere info in ${GENOME} database, dumb entry will be used!"
		rm ${RA_PATH}/${GENOME}_telomere.bed
		cp -f ${PRO_DATA_PATH}/spec_anno/dumb/dumb_subtelomere.bed ${RA_PATH}/${GENOME}_telomere.bed
		cp -f ${PRO_DATA_PATH}/spec_anno/dumb/dumb_subtelomere.bed ${RA_PATH}/${GENOME}_subtelomere.bed
	fi
else
	echo "Using predefined pericentromere and subtelomere information..."
	cp -f ${PRO_DATA_PATH}/spec_anno/${GENOME}/${GENOME}_centromere.bed ${RA_PATH}/${GENOME}_centromere.bed
	cp -f ${PRO_DATA_PATH}/spec_anno/${GENOME}/${GENOME}_telomere.bed ${RA_PATH}/${GENOME}_telomere.bed
	cp -f ${PRO_DATA_PATH}/spec_anno/${GENOME}/${GENOME}_pericentromere.bed ${RA_PATH}/${GENOME}_pericentromere.bed
	cp -f ${PRO_DATA_PATH}/spec_anno/${GENOME}/${GENOME}_subtelomere.bed ${RA_PATH}/${GENOME}_subtelomere.bed
fi

# gene gaps: for gene desert
grep genebody ${GENOME_ANNO_DIR}/${GENOME}.ensembl.biotype.txt | grep protein_coding | cut -f 1,2,3 | \
	complementBed -i - -g ${RA_PATH}/${GENOME}.genome \
	> ${RA_PATH}/${GENOME}_gap.bed

echo "Generating SubTelomere annotation..."
if [ ! -f "${RA_PATH}/${GENOME}_subtelomere.bed" ]; then
	# subtelomere: extend 500K bp of telomere, then subtract telomere regions, 
	# and remove the regions of genes (and 10K bp of flank regions of genes)
	grep genebody ${GENOME_ANNO_DIR}/${GENOME}.ensembl.biotype.txt | grep protein_coding | cut -f 1,2,3 | \
	slopBed -i - -g ${RA_PATH}/${GENOME}.genome -b 1e4 \
		> ${RA_PATH}/${GENOME}_gene_ext10k.bed
	slopBed -i ${RA_PATH}/${GENOME}_telomere.bed -g ${RA_PATH}/${GENOME}.genome -b 5e5 | \
		subtractBed -a - -b ${RA_PATH}/${GENOME}_telomere.bed | \
		subtractBed -a - -b ${RA_PATH}/${GENOME}_gene_ext10k.bed > ${RA_PATH}/${GENOME}_subtelomere.bed
fi

echo "Generating PeriCentromere annotation..."
if [ ! -f "${RA_PATH}/${GENOME}_pericentromere.bed" ]; then
	# pericentromere: gaps overlapped with centromere, remove the centromere regions, and only keep the regions longer than 10K bp.
	intersectBed -a ${RA_PATH}/${GENOME}_gap.bed -b ${RA_PATH}/${GENOME}_centromere.bed -wa | \
		subtractBed -a - -b ${RA_PATH}/${GENOME}_centromere.bed | \
		subtractBed -a - -b ${RA_PATH}/${GENOME}_telomere.bed | \
		awk '{if($3-$2>1e4){print $0}}' > ${RA_PATH}/${GENOME}_pericentromere.bed
fi

echo "Generating GeneDesert annotation..."
# # gene desert: gaps, remove the regions of telomere, subtelomere, centromere, and pericentromere, 
# # and only keep the regions longer than 1M bp, and subtract 10K from the flank regions of the nearest genes.
awk 'BEGIN{FS="\t"};{if(($3-$2)>=1e6){print $0}}' ${RA_PATH}/${GENOME}_gap.bed | \
	subtractBed -a - -b ${RA_PATH}/${GENOME}_pericentromere.bed | \
	subtractBed -a - -b ${RA_PATH}/${GENOME}_subtelomere.bed | \
	subtractBed -a - -b ${RA_PATH}/${GENOME}_telomere.bed | \
	subtractBed -a - -b ${RA_PATH}/${GENOME}_centromere.bed | \
	awk 'BEGIN{OFS="\t"}{if(($3-$2)>=1e6){print $1, $2+1e4, $3-1e4}}' > ${RA_PATH}/${GENOME}_geneDesert.bed # 10k from the nearest genes

echo "Generating genome annotation..."
for i in ${NP_DB_PATH}/tmp/${GENOME}*.biotype.txt; do
	## fulfill the empty gene symbol column with "NA"
	grep genebody ${i} | sed 's/^\t/NA\t/;:a;s/\t\t/\tNA\t/g;ta;s/\t$/\tNA/' | \
		egrep "^chr" | grep -v "random" | egrep -v "^NT" > ${i/txt/bed}
	cut -f2,3 ${i/txt/bed} | paste ${i/txt/bed} - > ${NP_DB_PATH}/tmp/${GENOME}.temp.bed
	slopBed -l 3e3 -r 1e3 -g ${RA_PATH}/${GENOME}.genome -s -i ${NP_DB_PATH}/tmp/${GENOME}.temp.bed | \
		awk '{if(($7=="+")&&($3>$2)&&($3>0))print $0}' > ${NP_DB_PATH}/tmp/${GENOME}.temp_ext.bed
	slopBed -l 1e3 -r 3e3 -g ${RA_PATH}/${GENOME}.genome -s -i ${NP_DB_PATH}/tmp/${GENOME}.temp.bed | \
		awk '{if(($7=="-")&&($3>$2)&&($3>0))print $0}' >> ${NP_DB_PATH}/tmp/${GENOME}.temp_ext.bed
	mv ${NP_DB_PATH}/tmp/${GENOME}.temp_ext.bed ${i/biotype.txt/biotype_region_ext.bed}
	rm ${NP_DB_PATH}/tmp/${GENOME}.temp.bed
done

rm ${RA_PATH}/${GENOME}_gap.bed
if [ -f "${RA_PATH}/${GENOME}_gene_ext10k.bed" ]; then
	rm ${RA_PATH}/${GENOME}_gene_ext10k.bed
fi
mv ${NP_DB_PATH}/tmp/${GENOME}*biotype_region_ext.bed ${RA_PATH}/
# mv ${GENOME}*.bed ${RA_PATH}/

${PRO_BIN_PATH}/genRA_json.py ${GENOME} ${NPVER} ${SPECIES} ${ASSEMBLY} \
	${ENSVER} ${RA_PATH} ${PRO_DATA_PATH}
CUR_PATH=`pwd`
cd ${RA_PATH}/..
echo "Compressing annotation package for Region_Analysis..."
tar czvf ${GENOME}_${NPVER}_En${ENSVER}.tar.gz ${GENOME}_${NPVER}_En${ENSVER}
echo "Installing package for Region_Analysis..."
region_analysis_db.py install ${GENOME}_${NPVER}_En${ENSVER}.tar.gz
cd ${CUR_PATH}
