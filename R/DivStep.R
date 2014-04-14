# TODO: Add comment
# 
# Author: przemol
###############################################################################
rlel2int <- function(v) {
	eval(parse( text = paste('c(as.integer(', deparse(substitute(v)), '[["',  paste(names(v), collapse=paste('"]]), as.integer(', deparse(substitute(v)), '[["', sep='')), '"]]))', sep='') ))
}

DivStep <- function(sample, control, enriched_regions=NULL) {
	
	## if(class(control) == 'chracter') {
	##     if (control == "formaldehyde") {
	##         catTime("Loading default formaldehyde track", e={			
	##             varname <- load("/Users/przemol/Documents/workspace/RBeads/rBeads/precalculated/summed_formaldehyde_input.Rdata")
	##         })
	## 
	##         catTime("Performing DIV step", e={
	##             formaldehyde_median <- 189 		
	##             aa <- formaldehyde_median / median(as.numeric(sample), na.rm=T)
	##             divstep <- sample / (summed_formaldehyde_input / aa)
	##         })		
	##     } else if (control == "EGS") {
	##         catTime("Loading default EGS track", e={			
	##             varname <- load("/Users/przemol/Documents/workspace/RBeads/rBeads/precalculated/summed_EGS_input.Rdata")
	##         })
	## 
	##         catTime("Performing DIV step", e={
	##             EGS_median <- 76 		
	##             aa <- EGS_median / median(as.numeric(sample), na.rm=T)
	##             divstep <- sample / (summed_EGS_input / aa)
	##         })
	##     }
	## } else {		
		## if (!is.null(enriched_regions)) {
		##     catTime("Using custom input track and performing DIV step with enreached region", e={
		##         #Calculate ogical vectors of non-enriched regions
		##         nonEnrichedRegionsLogi <- !coverage(GRanges(space(enriched_regions), unlist(ranges(enriched_regions)), "*", seqlengths=seqlengths(Celegans)))
		##         aa <- median(as.numeric(control), na.rm=T) / median(as.numeric(sample[nonEnrichedRegionsLogi[c(1,2,3,4,7,5,6)]]), na.rm=T)
		##     })
		## } else {
		##     catTime("Using custom input track amd performing DIV step - NO enreached region found!", e={
		##         aa <- median(as.numeric(control), na.rm=T) / median(as.numeric(sample), na.rm=T)
		##     })
		## }
		## catTime("Divideing and scaling", e={	
		##     divstep <- sample / (control / aa)		
		## })
	## }
	
	catTime("DivStep: Divideing", e={	
				divstep <- sample / control[names(sample)]
			})
	
	catTime("DivStep [median2x]: Scaling", e={
				aa <- median(median(divstep, na.rm=T))
				divstep <- divstep / aa
			})
	cat('INFO: Scaling coefficient = ', aa, '\n', sep='')
	
	return(round(divstep, 3))
	
}

catTime <- function(..., e=NULL, file="", gc=FALSE) {
	cat(..., "...", sep="", file=file, append=TRUE)
	cat("\t<", system.time(e, gcFirst=gc)[3], "s>\n", sep="", file=file, append=TRUE)	
}




