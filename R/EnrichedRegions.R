# Przemyslaw Stempor, 2014
##############################################################################
EnrichedRegions <- function(ranges.raw, REF=NULL ) {
	
	#INTEGRATED "PEAK CALLER"
	
	#Calculate coverage
	catTime("Calculate coverage", e={
    
				combExtCoverRep1 <- coverage(ranges.raw)
				sorted <-  lapply(combExtCoverRep1, sort)
				tabulation <- lapply(sorted, function(y) {
				  cbind(value=as.numeric(runValue(y)), length=as.numeric(runLength(y)))
				})
				tabulation <- do.call(rbind, tabulation)
				tabulation <- tapply(tabulation[,'length'], tabulation[,'value'], sum)
				a <- as.integer(names(tabulation))[ which.max( cumsum(tabulation) / sum(tabulation) >= .75 ) ]
        
			})	
	message('INFO: a = ', a)
  
	#Do peak calling
	catTime("Calling enriched regions", e={
				enriched_regions <- slice(combExtCoverRep1, lower = a )
				peakSumsRep1 <- viewSums(enriched_regions)
				enriched_regions <- RangedData(as( enriched_regions[peakSumsRep1 >= quantile(unlist(peakSumsRep1), .90)] , "IRangesList"))
				enriched_regions <- as( enriched_regions , 'GRanges' )
        seqinfo(enriched_regions) <- seqinfo(ranges.raw)[seqlevels(enriched_regions)]
        #if( !is.null(REF) ) seqlengths(enriched_regions) <- seqlengths(REF)[seqlevels(enriched_regions)]
			})
	
	return(enriched_regions)
}


