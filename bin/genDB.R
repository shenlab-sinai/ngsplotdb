make_gene_model <- function(txtfile){
	stopifnot(is.character(txtfile))
	# there are gene names with quote in Tair10!
	gene.tbl <- try(read.table(txtfile, comment.char='', sep="\t", as.is=T,  quote="",
		col.names=c('chrom', 'start', 'end', 'gid', 'gname', 'tid', 'strand', 'region', 'class', 'biotype')))
	if(class(gene.tbl) == "try-error"){
		gene.tbl <- read.table(txtfile, comment.char='', sep="\t", as.is=T, 
			col.names=c('chrom', 'start', 'end', 'gid', 'gname', 'tid', 'strand', 'region', 'class'))
		gene.tbl['biotype'] <- 'all'
	}


	# genebody.
	genebody <- subset(gene.tbl, region=='genebody')
	# choose longest transcript for each gene name.
	gb.byname.sorted <- genebody[order(genebody$gname, (genebody$end - genebody$start), decreasing=T), ]
	byname.uniq.tid <- gb.byname.sorted[!duplicated(gb.byname.sorted$gname), ]$tid
	genebody <- cbind(genebody[, c(1:7, 10)], byname.uniq=(genebody$tid %in% byname.uniq.tid))
	# choose longest transcript for each gene id. 
	if(all(is.na(genebody$gid))){
		bygid.uniq.tid <- NULL
		genebody <- cbind(genebody, bygid.uniq=NA)
	}else{
		gb.bygid.sorted <- genebody[order(genebody$gid, (genebody$end - genebody$start), decreasing=T), ]
		bygid.uniq.tid <- gb.bygid.sorted[!duplicated(gb.bygid.sorted$gid), ]$tid
		genebody <- cbind(genebody, bygid.uniq=(genebody$tid %in% bygid.uniq.tid))
	}
	genebody.biotype <- split(genebody, genebody$biotype)
	
	# exon.
	exon <- subset(gene.tbl, region=='exon')
	exon <- cbind(exon[, -8], byname.uniq=(exon$tid %in% byname.uniq.tid))
	if(all(is.na(exon$gid))){
		exon <- cbind(exon, bygid.uniq=NA)
	}else{
		exon <- cbind(exon, bygid.uniq=(exon$tid %in% bygid.uniq.tid))
	}
	exon.class <- split(exon, exon$class)
	exon.class.biotype <- lapply(exon.class, function(x) split(x, x$biotype))

	# exonmodel.
	exon.tid <- split(exon, exon$tid)
	library(IRanges)
	exonmodel <- lapply(exon.tid, function(tid) {
		tid <- tid[order(tid$start), ]
		list(ranges=IRanges(start=tid$start, end=tid$end))
	})
	exonmodel <- exonmodel[match(genebody$tid, names(exonmodel))]
	
	# return a nested list.
	list(genebody=genebody.biotype,
		exon=exon.class.biotype,
		exonmodel=exonmodel)
}

make_cgi_list <- function(cgi.anno, database, gene.anno=NULL){
	stopifnot(database == 'refseq' || (database == 'ensembl' && !is.null(gene.anno)))
	# Read cpg islands.
	cgi.cols <- c('chrom','start','end',
		'gname','tid','tstrand','tstart','tend','feature','dis', 'biotype')
	if(database == 'ensembl'){
		cgi.cols[4] <- 'gid'
	}
	cgi <- read.table(cgi.anno, sep="\t", col.names=cgi.cols, row.names=NULL, stringsAsFactors=FALSE, comment.char='')
	cgi[cgi$tstrand == "", ]$tstrand <- '+'	# set transcript strand to "+" for intergenic cgi.
	# Format cgi table.
	if(database == 'ensembl'){
		gid_name.tbl <- read.table(gene.anno, sep="\t", stringsAsFactors=FALSE, comment.char='')
		gid_name.tbl <- unique(gid_name.tbl[, c(4,5)])	# col 4: gid; 5: gname.
		colnames(gid_name.tbl) <- c('gid','gname')
		cgi <- merge(cgi, gid_name.tbl, all.x=TRUE, sort=FALSE)
		cgi.tbl <- data.frame(chrom=cgi$chrom, start=cgi$start, end=cgi$end,
			gid=cgi$gid, gname=cgi$gname, tid=cgi$tid, strand=cgi$tstrand, biotype=cgi$biotype,
			byname.uniq=TRUE, bygid.uniq=TRUE)
	}else{
		cgi.tbl <- data.frame(chrom=cgi$chrom, start=cgi$start, end=cgi$end,
			gid=NA, gname=cgi$gname, tid=cgi$tid, strand=cgi$tstrand, biotype=cgi$biotype,
			byname.uniq=TRUE, bygid.uniq=NA)
	}
	# Split cgi table based on feature annotation.
	cgi.class <- split(cgi.tbl, cgi$feature)
	cgi.class.biotype <- lapply(cgi.class, function(x) split(x, x$biotype))
	cgi.class.biotype
}

SplitDB <- function(d, s, g){
# Split the original ngs.plot databases into smaller ones based on functions.
# This aims to improve the efficiency in loading genomic coordinates.
# Args:
#   d: database object.
#   s: species name.
#   g: database name, such as refseq or ensembl.
	require(foreach)

	# genebody
	foreach(biotype=names(d$genebody)) %do% {
    	d$genebody[[biotype]] -> genome.coord
    	genome.coord[order(!genome.coord$byname.uniq, !genome.coord$bygid.uniq, genome.coord$gname), ] -> genome.coord
    	out.f <- sprintf("%s.%s.genebody.%s.RData", s, g, biotype)
    	save(genome.coord, file=out.f)
	}
   
    
	# exon
    foreach(n=names(d$exon)) %do% {
    	foreach(biotype=names(d$exon[[n]])) %do% {
        	d$exon[[n]][[biotype]] -> genome.coord
        	genome.coord[order(!genome.coord$byname.uniq, !genome.coord$bygid.uniq, genome.coord$gname), ] -> genome.coord
        	out.f <- sprintf("%s.%s.exon.%s.%s.RData", s, g, n, biotype)
        	save(genome.coord, file=out.f)
    	}
	}

	# cgi
    foreach(n=names(d$cgi)) %do% {
    	foreach(biotype=names(d$cgi[[n]])) %do% {
        	d$cgi[[n]][[biotype]] -> genome.coord
        	genome.coord[order(!genome.coord$byname.uniq, !genome.coord$bygid.uniq, genome.coord$gname), ] -> genome.coord
        	out.f <- sprintf("%s.%s.cgi.%s.%s.RData", s, g, n, biotype)
        	save(genome.coord, file=out.f)
    	}
    }
    
	# exonmodel
    d$exonmodel -> exonmodel
    out.f <- sprintf("%s.%s.exonmodel.RData", s, g)
    save(exonmodel, file=out.f)
}

args <- commandArgs(TRUE)
species <- args[1]
db.name <- args[2]
anno.file <- args[3]

anno <- make_gene_model(anno.file)
cgi.anno <- paste('./tmp/',species,'.cgi.',db.name,'.biotype.txt', sep='')
if(file.exists(cgi.anno)){
	if(db.name=='ensembl'){
		anno$cgi <- make_cgi_list(cgi.anno, db.name, anno.file)
	}
	if(db.name=='refseq'){
		anno$cgi <- make_cgi_list(cgi.anno, db.name)
	}
}
SplitDB(anno, species, db.name)