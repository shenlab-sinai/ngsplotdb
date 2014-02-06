It is the pipeline for annotations of ensembl and UCSC.
Copied and modified from our ensembl database pipeline.

Usage: 
Add ./bin to $PATH, then makedir a new folder, then run:
> genDB.sh db_list.txt [animal|plant|bacteria].json

Workflow:
1. Read pipeline configuration from [animal|plant|bacteria].json under ./json, and parsing the information of genomes from db_list.txt.
2. Download gene annotations from ensembl and UCSC, then generate the annotations in plain text and RData. If no annotation in UCSC, then only ensembl annotation was generated.
3. Generate meta-info data of the annotations.
4. Process the annotation files for region_analysis (if no UCSC Refseq annotation, this step will be skipped):
	A. Prepare the annotation fiels for pericentromere and subtelomere:
		a. If there are ready info of pericentromere and subtelomere, then use the ready files. If not, then:
		b. If UCSC has standard gap table of genome gap, then download centromere and telomere annotations and calculate pericentromere and subtelomeres. If not, then:
		c. Use the dumb files directly.
	B. Calculate the gene desert regions.
	C. Extend the annotations of genes to the regions needed by region_analysis: upstream 3k, downstream 1k.
	D. Install the new databases to region_analysis.
5. Download the CpG islands annotations from UCSC. If the result is not empty, then:
	A. Annotate CGI by region_analysis.
	B. Generate RData of CGI
6. Pack all annotations.

Attention:
1. Pipeline templates are under ./json folder. Now animal, plant, and bacteria are supported.
2. When re-generate the databases of "animal", the deletion of installed annotations of region_analysis is needed. Generally removing all files under ~/.config/regionsanalysis is enough.
3. Some genomes have wired nominations and take trouble for parsing. For now what I know:
	A. "#" in yeast sacCer3 gene names.
	B. ";" in rice IRGSP-1 transcript ids.
	C. "'" in Arabidopsis Tair10 gene names.
4. No UCSC Refseq annotation in zebrafish genome.
5. CORES, defined at the beginning of genDB.sh, is the threads used in the pipeline. Now I set it as 4.
6. NPVer is defined at the beginning of genDB.sh.


ENCODE_cellline_dhs:
For the annotation of ENCODE DHS regions.
Multiple threads are supported, default is 4.

TODO:
	* Now it needs hg19.ensembl.biotype.txt to get gene name, it is annoying.