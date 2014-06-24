# Przemyslaw Stempor, 2014
##############################################################################

MappabilityCorrection <- function(GCnormTrack, mappabilityTrack=NULL, cutoff=0.25) {
	
  # 1) masking regions of mappability lower than cuttof value [mappability score=0.25] from precalculated mappability track
	catTime("Mappability correction", e={			
	  GCnormTrack[mappabilityTrack[names(GCnormTrack)] < cutoff] <- NA	
	})
  
  message("\tINFO: Total ", round(sum(as.numeric(sum(is.na( GCnormTrack )))) / sum(as.numeric(elementLengths( GCnormTrack )))*100, 2), "% of genome in masked.")
	return(GCnormTrack)
  
}