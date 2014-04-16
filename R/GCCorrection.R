# TODO: Add comment
# 

if (0) {
	
	testGCscores <- function( files=dir(pattern="RawRanges.+Rdata") ) {
		source('/Users/przemol/Documents/workspace/RBeads/rBeads/rB_enriched_regions_one_strand_v1.0.R')
		out = list()
		for(file in files) {
			cat("-> GC TESTING:", file, "\n")
			varname <- load(file)
			d <- sprintf("%s", unlist(strsplit(file, "\\."))[1])
			ER <- rB.EnrichedRegions.OS(get(varname), desc=d )
			GC <- rB.sumGCscores(get(varname), enriched_regions=ER, desc=d, smoothing_spline=FALSE)
			out <- c(out, GC)
		}
		names(out) <- files
		return(out)
	}
	
	
}
# Author: przemol
##############################################################################
require(GenomicRanges)
require(rtracklayer)
require(BSgenome)

GCCorrection <- function(ranges.raw, enriched_regions, nonMappableFilter, genome_fasta=NULL, genome_package=NULL, resize_length=200L desc='', smoothing_spline=FALSE, cutoff=c(35, 140)) {
	
  if( !is.null(genome_fasta) ) {
    REF <- readDNAStringSet( genome_fasta )
  } else {
    package <- grep(genome_package, installed.genomes(), value=TRUE, ignore.case=TRUE)
    if( !length(package) ) { stop('Genome ', genome_package, ' is not installed. Please run "available.genomes()" to get genomes supported by BioConductor or provider reference FASTE file.', call.=FALSE) }
    library(package, character.only = TRUE)
    REF <- getSeq( get(package) ); names(REF) <- seqnames(get(package))
  }
  
  
	#Mask out reads in enriched regions
	if (!is.null(enriched_regions)) {
		catTime("Mask out reads in enriched regions", e={
			ERoverlaps <- findOverlaps(ranges.raw, enriched_regions, select="first") 	
		})
		#INFO: percentage of reads in enriched regions
		cat("\tINFO: percentage of reads in non-enriched regions: ", round(100 * sum(is.na(ERoverlaps))/length(ranges.raw), 2), "%\n", sep='')
	}
	
	#Calculate input GC content from 200bp extended reads (ranges.raw)
	catTime("Calculate input GCcontent", e={
		GCcontent <- as.integer(letterFrequency(REF[ranges.raw], "GC"))
	})
	
	#Sample genome for GCcontent
	catTime("Sample genome for GCcontent", e={
		if (!is.null(enriched_regions)) {
			#Calculate ogical vectors of non-enriched regions
			nonEnrichedRegionsLogi <- !coverage(GRanges(space(enriched_regions), unlist(ranges(enriched_regions)), "*", seqlengths= seqlengths( REF ) ))
			#Perform logical sum of non-enriched regions and mappable regions
			nonEnrichedMappableRegionsLogi <- nonEnrichedRegionsLogi & nonMappableFilter
		} else {
			nonEnrichedMappableRegionsLogi <- nonMappableFilter
		}
		#Calculate GC pecrentage among the chromosomes
		GC <- RleList( lapply(REF, function(x) Rle(letterFrequencyInSlidingView(x, 200, "GC") )) )
    
		#Select only on-enriched regions and mappable regions
		GCgenome <- GC[nonEnrichedMappableRegionsLogi]	
	})
	
	#INFO: percentage of genome to be sampled
	cat("\tINFO: percentage of genome to be sampled: ", sum(elementLengths(GCgenome)) / sum(seqlengths( REF )), "\n") 

	#Calculate histograms for genomic (a) nad and sample (b) GC content
	catTime("Calculate histograms for genomic (a) nad and sample (b) GC content", e={
		a <- hist(as.integer(unlist(GCgenome, use.names=FALSE)), 0:200, plot=FALSE)
		if (!is.null(enriched_regions)) {
			b <- hist(GCcontent[is.na(ERoverlaps)], 0:200, plot=F)
		} else {
			b <- hist(GCcontent, 0:200, plot=FALSE)
		}
	
		pdf(sprintf("IMG - GCdistribution - %s.pdf", desc), width = 12.0, height = 7.5, onefile = FALSE, paper = "special", encoding = "TeXtext.enc")
		plot(b$mids, b$density, col="red", type="l", main=sprintf("GCdistribution - %s", desc))
		lines(b$mids, a$density, type="l", col="blue")
		legend("topleft", c("Non enriched GENOMIC GC content distribution", "Non enriched SAMPLE GC content distribution"), fill=c("blue", "red") )
		dev.off()
	})
	
	cat('\tINFO: scales <- c(', paste(a$density/b$density, collapse=','), ')\n', sep='')
	
	#Calculate GC weighting vector
	catTime("Calculate GC weighting vector", e={												
		scales <- a$density/b$density
	
		pdf(file=sprintf("IMG - GCoutlier - %s.pdf", desc), width = 12.0, height = 7.5, onefile = FALSE, paper = "special", encoding = "TeXtext.enc")
		plot(scales, main=sprintf("GCoutlier - %s", desc))
		
		#Remove outliers from the vecor
		#Grubbs test
		#while (grubbs.test(scales[,2])$p < 0.1 ) {
		#	scales <- scales[-which.max(scales[,2]), ]
		#}
		#Alternative ways of doing that
		#1) Quantile
		#scales <- scales[which(scales[,2] < quantile(scales[,2], .9)), ]
		#2) Boxplot
		scales[ scales %in% boxplot(scales, plot=F)$out ] <- NA
		#3) Fixed
		scales[c( 1:(cutoff[1]-1), (cutoff[2]+1):200 )] <- NA
		abline(h=0, v=c(cutoff[1], cutoff[2]), col = "gray60")
		
		
		#Roboust prediction of scaling ceficient by fitting a cubic smoothing spline to the supplied data.
		SSpoints <- cbind(1:200, scales)[!is.na(scales) & !(scales %in% boxplot(scales, plot=F)$out), ]	
		SS <- smooth.spline(SSpoints, df=10)
		scales.fit <- cbind(predict(SS, 0:200)$x, predict(SS, 0:200)$y)
		lines(scales.fit, col="green")
		
		if (smoothing_spline) {
			scales.f <- scales.fit
		}
		points(scales, col="red", pch="*")
		points(SSpoints, col="green", pch="o")
		dev.off()
	})
	
	#Assign GC score to every read
	catTime("Assign GC score to every read", e={
		GCcontent[GCcontent == 0] <- 1
		scores <- scales[GCcontent]
	})
	
	#Prepare GCscore-enriched GRanges
	catTime("Prepare GCscore-enriched GRanges", e={				
		values(ranges.raw) <- scores
		ranges.f1 <- ranges.raw[!is.na(values(ranges.raw)),]
	})
	
	#calculate GC weighted coverage with given accuracy
  # TODO; simplify it and chack on real life example
	catTime("calculate GC weighted coverage with given accuracy", e={
		acc = 1000
		w <- sapply( seqlevels(ranges.f1), function(x) values( ranges.f1[seqnames(ranges.f1) == x] )$value * acc )
		#Output coverage
		cov <- coverage(ranges.f1, weight=w) / acc
		cov.r <- round(cov)
	})

	#Mask nonGC correctable regions
	catTime("Masking non-GCcorrectable regions", e={
		notGCcorrectableReads <- IntegerList( sapply( names( REF ), function(x) {  
			GCchr <- letterFrequencyInSlidingView(REF[[x]], 200, "GC")
			which( ! (GCchr >= cutoff[1] & GCchr <= cutoff[2]) )
		}) )
		#notGCcorrectableRegions <- GRanges(seqnames=as.data.frame(notGCcorrectableReads)$space, ranges=IRanges(as.data.frame(notGCcorrectableReads)$value, width=200))
		notGCcorrectableRegions <- GRanges(seqnames=Rle(names(notGCcorrectableReads), sapply(notGCcorrectableReads, length)), ranges=IRanges(unlist(notGCcorrectableReads), width=200))
		#seqlengths(notGCcorrectableRegions) <- seqlengths( REF )
		notGCcorrectableRegions <- c(notGCcorrectableRegions, GRanges(seqnames=seqlevels( REF ), ranges=IRanges( seqlengths( REF )-199, width=200) ))
		cov.r[ coverage(notGCcorrectableRegions)[names(cov.r)] > 0 ] <- NA
	})

	return(cov.r)
}

catTime <- function(..., e=NULL, file="", gc=FALSE) {
	cat(..., "...", sep="", file=file, append=TRUE)
	cat("\t<", system.time(e, gcFirst=gc)[3], "s>\n", sep="", file=file, append=TRUE)	
}
