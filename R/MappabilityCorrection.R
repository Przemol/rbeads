# Przemyslaw Stempor, 2014
##############################################################################

MappabilityCorrection <- function(GCnormTrack, mappabilityTrack=NULL, cutoff=0.25) {
	
  # 1) masking regions of mappability lower than 100 (mappability score=0.25) from precalculated mappability track
	catTime("Mappability correction", e={			
	  GCnormTrack[mappabilityTrack[names(GCnormTrack)] < 100] <- NA	
	})
	return(GCnormTrack)
  
}