# Przemyslaw Stempor, 2014
##############################################################################

DivStep <- function(sample, control) {

	catTime("DivStep: Divideing", e={	
        commonChr <- intersect(names(sample), names(control))
				divstep <- round(sample[commonChr] / (control[commonChr]+1), 3)
			})
	
	catTime("DivStep [median ]: Scaling", e={
	  
	  sorted <-  lapply(divstep[!is.na(divstep)], sort)
	  tabulation <- lapply(sorted, function(y) {
	    cbind(value=as.numeric(runValue(y)), length=as.numeric(runLength(y)))
	  })
	  tabulation <- do.call(rbind, tabulation)
	  tabulation <- tapply(tabulation[,'length'], tabulation[,'value'], sum)
	  #tabulation <- tabulation[names(tabulation)!=0]
	  aa <- as.numeric(names(tabulation))[ which.max( cumsum(tabulation) / sum(tabulation) >= .5 ) ]
    
    divstep <- round(divstep / aa, 3)
	})
	message('INFO: BEADS score scaling coefficient = ', aa, '\n', sep='')
	
	return(divstep)
	
}


