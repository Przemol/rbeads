# TODO: Add comment
# 
# Author: przemol
###############################################################################


# TODO: Add comment
# 
# Author: przemol
###############################################################################


# TODO: Add comment
# 
# Author: przemol
if(0){
	files <- sapply( dir(pattern="Ranges.+Rdata"), function(x) {cat("Loading:", x, "\n"); return(load(x, envir=.GlobalEnv))} )
	gsub("(\\W|get)", "_", deparse(substitute(get("a"))), perl=T)
	
	ERandBin <- function(files) {
		for(file in files) {
			desc <- sprintf("%s_OS", unlist(strsplit(file, "\\."))[1])
			varname <- desc
			assign(varname, importBAM(file, desc=desc, resize_length=200, quality_cutoff=10, export_bin=TRUE, export_track=TRUE)) 
				
			#varname <- load(file)
			rB.EnrichedRegions.OS(get(varname), desc=desc )
			#binTrack(coverage(get(varname)), n=25, smooth=FALSE, out=sprintf("%s.wig", desc), type="WIG")
		}
	} 
	ERandBin( dir(pattern="RawRanges.+Rdata") )
	ERandBin( dir(pattern="bam") )
	
}
#	
#
###############################################################################
require(GenomicRanges)
require(rtracklayer)
#library(multicore)

EnrichedRegions <- function(ranges.raw, desc="EnrichedRegions1" ) {
	
	
	cat("GC correction - peak calling", "\n")
	
	#INTEGRATED PEAK CALLER
	
	#Calculate coverage
	catTime("Calculate coverage", e={
				combExtCoverRep1 <- coverage(ranges.raw)
			})	
	
	#Do peak calling
	catTime("Do peak calling", e={
				a = mean(quantile(combExtCoverRep1, .75))
				enriched_regions <- slice(combExtCoverRep1, lower = a)
				peakSumsRep1 <-viewSums(enriched_regions)
				enriched_regions <- RangedData(as( enriched_regions[peakSumsRep1 >= quantile(peakSumsRep1[[3]], .90)] , "IRangesList"))
			})
	
	#INFO: prepare the peak calling bed file to view in IGB/IGV
	#catTime("INFO: prepare the peak calling bed file to view in IGB/IGV", e={															
	#			export.bed(enriched_regions, "ERpeakCall.bed")	
	#		})
	return(enriched_regions)
}


