catTime <-
function(..., e=NULL, file="", gc=FALSE) {
	cat(..., "...", sep="", file=file, append=TRUE)
	cat("\t<", system.time(e, gcFirst=gc)[3], "s>\n", sep="", file=file, append=TRUE)	
}

toBW <- function(input, proc, name, missing2zero=FALSE) {
  message('Exporting to BigWiggle: ', reName(name, proc, 'linear', '1bp', '.bw'))
  if(class(input) == 'GRanges') {
    message('Calculating coverage from ranges')
    input <- coverage(input)
  }
  if(missing2zero) {
    input[is.na(input)] <- 0
  }
  export.bw(input, reName(name, proc, 'linear', '1bp', '.bw') )
}

toBW_missing <- function(input, proc, name, resolution=1L) {
  message('Exporting to BigWiggle: ', reName(name, proc, 'linear', paste0(resolution, 'bp'), '.bw'))
  if(class(input) == 'GRanges') {
    message('Calculating coverage from ranges')
    input <- coverage(input)
  }
  gr <- as(input, 'GRanges')
  gr <- gr[!is.na(gr$score)]
  gr <- gr[!is.infinite(gr$score)]
  export.bw(gr, reName(name, proc, 'linear', paste0(resolution, 'bp'), '.bw') )
}

getREF <- function(genome) {
  if( file.exists(genome) ) {
    REF <- readDNAStringSet( genome )
  } else {
    package <- grep(genome, installed.genomes(), value=TRUE, ignore.case=TRUE)
    if( !length(package) ) { stop('Genome ', genome, ' is not installed. Please run "available.genomes()" to get genomes supported by BioConductor or provider reference FASTE file.', call.=FALSE) }
    library(package, character.only = TRUE)
    REF <- getSeq( get(package) ); names(REF) <- seqnames(get(package))
  }
  return(REF)
}

binRle2Rle <- function(input, bin=25L, precision=0L) {
  out <- RleList(sapply(names(input), function(chr) {
    br <- breakInChunks(length(input[[chr]]), bin)
    sc <- round(viewMeans(Views(input[[chr]], br), na.rm = TRUE), precision)
    Rle(sc, width(br))
  }))
  return(out)
}