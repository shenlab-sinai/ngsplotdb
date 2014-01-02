It is the pipeline for animal annotations of ensembl and UCSC.
Copied and modified from ensembl.

Workflow:
1. Download gene annotations from ensembl and UCSC, then generate the annotations in plain text and RData.
2. Generate meta-info data of the annotations.
3. Process the annotation files for region_analysis:
	A. Prepare the annotation fiels for pericentromere and subtelomere:
	a. If there are ready info of pericentromere and subtelomere, then use the ready files. If not, then:
	b. If UCSC has standard gap table of genome gap, then download centromere and telomere annotations and calculate pericentromere and subtelomeres. If not, then:
	c. Use the dumb files directly.
	B. Calculate the gene desert regions.
	C. Extend the annotations of genes to the regions needed by region_analysis: upstream 3k, downstream 1k.
	C. Install the new databases to region_analysis.
4. Download the CpG islands annotations from UCSC. If the result is not empty, then:
	A. Annotate CGI by region_analysis.
	B. Generate RData of CGI
5. Pack all annotations.

Attention:
1. If region_analysis installed under root, then "sudo" is needed.