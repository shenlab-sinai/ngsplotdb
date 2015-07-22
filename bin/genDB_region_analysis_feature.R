#! /bin/env Rscript

make_feature_list <- function(feature.anno, database, gene.anno=NULL){
	stopifnot(database == 'refseq' || (database == 'ensembl' && !is.null(gene.anno)))
	# Read features.
	feature.cols <- c('chrom','start','end',
		'gname','tid','tstrand','tstart','tend','feature','dis', 'biotype', 'gsymbol')
	if(database == 'ensembl'){
		feature.cols[4] <- 'gid'
	}
	feature <- read.table(feature.anno, sep="\t", col.names=feature.cols, row.names=NULL, stringsAsFactors=FALSE, comment.char="")
	feature[feature$tstrand == ".", ]$tstrand <- '+'	# set transcript strand to "+" for intergenic feature.
	# Format feature table.
	if(database == 'ensembl'){
		gid_name.tbl <- read.table(gene.anno, sep="\t", stringsAsFactors=FALSE)
		gid_name.tbl <- unique(gid_name.tbl[, c(4,5)])	# col 4: gid; 5: gname.
		colnames(gid_name.tbl) <- c('gid','gname')
		feature <- merge(feature, gid_name.tbl, all.x=TRUE, sort=FALSE)
		feature.tbl <- data.frame(chrom=feature$chrom, start=feature$start, end=feature$end,
			gid=feature$gid, gname=feature$gname, tid=feature$tid, strand=feature$tstrand, biotype=feature$biotype,
			byname.uniq=TRUE, bygid.uniq=TRUE)
	}else{
		feature.tbl <- data.frame(chrom=feature$chrom, start=feature$start, end=feature$end,
			gid=NA, gname=feature$gname, tid=feature$tid, strand=feature$tstrand, biotype=feature$biotype,
			byname.uniq=TRUE, bygid.uniq=NA)
	}
	# Split feature table based on feature annotation.
	feature.class <- split(feature.tbl, feature$feature)
	feature.class.biotype <- lapply(feature.class, function(x) split(x, x$biotype))
	feature.class.biotype
}

SplitDB <- function(d, s, g, f){
# Split the original ngs.plot databases into smaller ones based on functions.
# This aims to improve the efficiency in loading genomic coordinates.
# Args:
#   d: database object.
#   s: species name.
#   g: database name, such as refseq or ensembl.
	require(foreach)

	# features annotated by region_analysis, like cgi, dhs
	foreach(n=names(d$feature)) %do% {
		foreach(biotype=names(d$feature[[n]])) %do% {
			d$feature[[n]][[biotype]] -> genome.coord
			if(dim(genome.coord)[1] != 0){
				genome.coord[order(!genome.coord$byname.uniq, !genome.coord$bygid.uniq, genome.coord$gname), ] -> genome.coord
				out.f <- sprintf("%s.%s.%s.%s.%s.RData", s, g, f, n, biotype)
				save(genome.coord, file=out.f)
			}
		}
	}
}

args <- commandArgs(TRUE)
species <- args[1]
feature <- args[2]
db.name <- args[3]

anno.file <- args[4]
feature.anno <- args[5]

anno <- list()
print(feature.anno)
if(file.exists(feature.anno)){
	if(db.name=='ensembl'){
		anno$feature <- make_feature_list(feature.anno, db.name, anno.file)
	}
	if(db.name=='refseq'){
		anno$feature <- make_feature_list(feature.anno, db.name)
	}
}
SplitDB(anno, species, db.name, feature)
