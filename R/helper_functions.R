catTime <-
function(..., e=NULL, file="", gc=FALSE) {
	cat(..., "...", sep="", file=file, append=TRUE)
	cat("\t<", system.time(e, gcFirst=gc)[3], "s>\n", sep="", file=file, append=TRUE)	
}


toBW <- function(grm, proc, name) {
  message('Exporting to BigWiggle: ', name)
  
  input <- get(grm)
  if(class(input) == 'GRanges') {
    message('Calculating coverage from ranges')
    input <- coverage(input)
  }
  
  export.bw(input, reName(name, proc, 'linear', '1bp', '.bw') )
}
