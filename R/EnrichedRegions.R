EnrichedRegions <-
function(ranges.raw, desc="EnrichedRegions1" ) {
	
	
	cat("GC correction - peak calling", "\n")
	
	#INTEGRATED PEAK CALLER
	
	#Calculate coverage
	catTime("Calculate coverage", e={
				combExtCoverRep1 <- coverage(ranges.raw)
			})	
	
	#Do peak calling
	catTime("Do peak calling", e={
				a = quantile(suppressWarnings(as.vector(combExtCoverRep1)), 0.75)
				enriched_regions <- slice(combExtCoverRep1, lower = a)
				peakSumsRep1 <-viewSums(enriched_regions)
				enriched_regions <- RangedData(as( enriched_regions[peakSumsRep1 >= quantile(peakSumsRep1[[3]], .90)] , "IRangesList"))
			})
	
	#INFO: prepare the peak calling bed file to view in IGB/IGV
	catTime("INFO: prepare the peak calling bed file to view in IGB/IGV", e={															
				export.bed(enriched_regions, sprintf("%s.ERpeakCall.bed", desc))	
			})
	return(enriched_regions)
}
