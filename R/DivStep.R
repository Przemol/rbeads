# Przemyslaw Stempor, 2014
##############################################################################

DivStep <- function(sample, control) {
	
	catTime("DivStep: Divideing", e={	
				divstep <- sample / control[names(sample)]
			})
	
	catTime("DivStep [median ]: Scaling", e={
				aa <- median( as.numeric(unlist(divstep, use.names=FALSE)), na.rm=T)
				divstep <- divstep / aa
			})
	cat('INFO: Scaling coefficient = ', aa, '\n', sep='')
	
	return(round(divstep, 3))
	
}




