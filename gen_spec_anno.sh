#! /bin/bash

GENOME=${1}
GENOME_ANNO_DIR="./tmp"

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
	"select chrom, size from "${GENOME}".chromInfo"  > ${GENOME}.genome

## if no ready info, try to fetch centromere and telomere info from UCSC
if [ ! -d "./spec_anno/${GENOME}" ]; then
	mkdir "./spec_anno/${GENOME}"

	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D ${GENOME} -e 'select chrom, chromStart, chromEnd from gap where type="centromere"' \
	| grep -v chromStart > ${GENOME}_centromere.bed
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D ${GENOME} -e 'select chrom, chromStart, chromEnd from gap where type="telomere"' \
	| grep -v chromStart > ${GENOME}_telomere.bed
	if [ ! -s "${GENOME}_centromere.bed" ]; then
		echo "No centromere info in ${GENOME} database, dumb entry will be used!"
		rm ${GENOME}_centromere.bed
		cp -f ./spec_anno/dumb/dumb_pericentromeres.bed ./spec_anno/${GENOME}/${GENOME}_centromere.bed
		cp -f ./spec_anno/dumb/dumb_pericentromeres.bed ./spec_anno/${GENOME}/${GENOME}_pericentromere.bed
	fi
	if [ ! -s "${GENOME}_telomere.bed" ]; then
		echo "No telomere info in ${GENOME} database, dumb entry will be used!"
		rm ${GENOME}_telomere.bed
		cp -f ./spec_anno/dumb/dumb_subtelomeres.bed ./spec_anno/${GENOME}/${GENOME}_telomere.bed
		cp -f ./spec_anno/dumb/dumb_subtelomeres.bed ./spec_anno/${GENOME}/${GENOME}_subtelomere.bed
	fi
fi

# gene gaps: for gene desert
grep genebody ${GENOME_ANNO_DIR}/${GENOME}.ensembl.biotype.txt | grep protein_coding | cut -f 1,2,3 | \
	complementBed -i - -g ${GENOME}.genome \
	> ${GENOME}_gap.bed

if [ ! -f "./spec_anno/${GENOME}/${GENOME}_subtelomere.bed" ]; then
	# subtelomere: extend 500K bp of telomere, then subtract telomere regions, 
	# and remove the regions of genes (and 10K bp of flank regions of genes)
	grep genebody ${GENOME_ANNO_DIR}/${GENOME}.ensembl.biotype.txt | grep protein_coding | cut -f 1,2,3 | \
	slopBed -i - -g ${GENOME}.genome -b 1e4 \
	> ${GENOME}_gene_ext10k.bed
	slopBed -i ${GENOME}_telomere.bed -g ${GENOME}.genome -b 5e5 | subtractBed -a - -b ${GENOME}_telomere.bed | \
	subtractBed -a - -b ${GENOME}_gene_ext10k.bed > ./spec_anno/${GENOME}/${GENOME}_subtelomere.bed
fi

if [ ! -f "./spec_anno/${GENOME}/${GENOME}_pericentromere.bed" ]; then
	# pericentromere: gaps overlapped with centromere, remove the centromere regions, and only keep the regions longer than 10K bp.
	intersectBed -a ${GENOME}_gap.bed -b ${GENOME}_centromere.bed -wa | subtractBed -a - -b ${GENOME}_centromere.bed | \
	subtractBed -a - -b ${GENOME}_telomere.bed | awk '{if($3-$2>1e4){print $0}}' > ./spec_anno/${GENOME}/${GENOME}_pericentromere.bed
fi

# # gene desert: gaps, remove the regions of telomere, subtelomere, centromere, and pericentromere, 
# # and only keep the regions longer than 1M bp, and subtract 10K from the flank regions of the nearest genes.
awk 'BEGIN{FS="\t"};{if(($3-$2)>=1e6){print $0}}' ${GENOME}_gap.bed | \
	subtractBed -a - -b ./spec_anno/${GENOME}/${GENOME}_pericentromere.bed | \
	subtractBed -a - -b ./spec_anno/${GENOME}/${GENOME}_subtelomere.bed | \
	subtractBed -a - -b ./spec_anno/${GENOME}/${GENOME}_telomere.bed | \
	subtractBed -a - -b ./spec_anno/${GENOME}/${GENOME}_centromere.bed | \
	awk 'BEGIN{OFS="\t"}{if(($3-$2)>=1e6){print $1, $2+1e4, $3-1e4}}' > ${GENOME}_geneDesert.bed # 10k from the nearest genes

for i in ./tmp/${GENOME}*.biotype.txt; do
	## fulfill the empty gene symbol column with "NA"
	grep genebody ${i} | sed 's/^\t/NA\t/;:a;s/\t\t/\tNA\t/g;ta;s/\t$/\tNA/' |egrep "^chr" | grep -v "random" | egrep -v "^NT" > ${i/txt/bed}
	cut -f2,3 ${i/txt/bed} | paste ${i/txt/bed} - > temp.bed
	slopBed -l 3e3 -r 1e3 -g ${GENOME}.genome -s -i temp.bed | awk '{if(($3>$2)&&($3>0))print $0}' > temp_ext.bed
	mv temp_ext.bed ${i/biotype.txt/biotype_region_ext.bed}
	rm temp.bed
done

rm ${GENOME}_gap.bed
if [ -f "${GENOME}_gene_ext10k.bed" ]; then
	rm ${GENOME}_gene_ext10k.bed
fi
mv ./tmp/${GENOME}*biotype_region_ext.bed ./spec_anno/${GENOME}/
mv ${GENOME}*.bed ./spec_anno/${GENOME}/
mv ${GENOME}.genome ./spec_anno/${GENOME}/

./RA_install_db.py ${GENOME}