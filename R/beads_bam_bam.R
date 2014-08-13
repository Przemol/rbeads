#' BEADS internal
#' 
#' The function runs BEADS for signature BAM, BW, BW, ANY: \cr
#' \code{signature(c(experiment='BamFile', control='BigWigFile', mappability='BigWigFile', genome='ANY'))}
#' 
#' This function is package internal and should not be executed directly
#' by users.
#' 
#' @keywords internal
#'    
beads_bam_bam <- function(bam.file, bam.control, bw.mappability, genome, uniq=TRUE, insert=200L, mapq_cutoff=10L, export='BEADS', rdata=FALSE, export_er=TRUE, quickMap=TRUE) {
  
  #Importing refference genome
  message('Importing refference genome...')
  REF <- getREF(genome)
  
  #Setting up mappability values
  if( quickMap ) {
    message('Importing mappability [binning]...')
    if( class(bw.mappability) == 'BigWigFile' ) { bwf <- bw.mappability } else { bwf <- BigWigFile( bw.mappability ) }
    MAP  <- IRanges::summary(bwf, size=GenomicRanges::seqlengths(bwf)/insert, asRle=TRUE)
    MAP[is.na(MAP)] <- 0
    MAPF <- ( MAP > 0.5 )
  } else {
    message('Importing mappability [running median smoothing]...')
    gem  <- import.bw( bw.mappability, as = "RleList" )
    MAP  <- round( runmean(gem, ifelse(insert %% 2 == 0, insert+1, insert), endrule = "constant"), 2 )
    MAP[is.na(MAP)] <- 0
    MAPF <- ( gem > 0.5 )
  }

  #Loading input from BigWig file
  message('Importing input...')
  control.re <-  ImportBAM(bam.control, REF=REF, uniq=uniq, resize_length=insert, quality_cutoff=mapq_cutoff)
  control.gc <-  GCCorrection(control.re, enriched_regions=NULL, REF=REF, nonMappableFilter=MAPF, RL=insert, desc=gsub('.bam$', '', basename(bam.control)), smoothing_spline=FALSE)
  #message('Smoothing...')
  #control.gc <-  round( runmean(control.gc, ifelse(insert %% 2 == 0, insert+1, insert), endrule = "constant"), 2 )
  control.map <- MappabilityCorrection(GCnormTrack=control.gc, mappabilityTrack=MAP)
  
  #Full beads
  sample.re <-   ImportBAM(bam.file, REF=REF, uniq=uniq, resize_length=insert, quality_cutoff=mapq_cutoff)
  sample.er <-   EnrichedRegions(sample.re, REF=REF)
  sample.gc <-   GCCorrection(sample.re, enriched_regions=sample.er, REF=REF, nonMappableFilter=MAPF, RL=insert, desc=gsub('.bam$', '', basename(bam.file)), smoothing_spline=FALSE)
  sample.map <-  MappabilityCorrection(GCnormTrack=sample.gc, mappabilityTrack=MAP)
  
  sample.norm <- DivStep(sample.map, control.map)
  
  if(export_er) {
    message('Exporing ER...'); er_con <- file(gsub('.bam$', '_EnrichedRegions.bed', basename(bam.file)))
    export.bed(sample.er, er_con); close(er_con)
  }
  if(rdata) { 
    message('Exporing Rdata binaries...')
    save( list = ls(), file = gsub('.bam$', '.Rdata', basename(bam.file)), envir = environment() ) 
  }
  
  message('Exporing BigWig tracks...')
  exp <- list('control.re', 'control.gc', 'control.map', 'sample.re', 'sample.gc', 'sample.map', 'sample.norm')
  names(exp) <- c('control_readsCoverage', 'control_GCcorected', 'control_GCandMap', 'readsCoverage', 'GCcorected', 'GCandMap', 'BEADS')
  bw_lst <- lapply( names(exp[export]), function(x) toBW_missing(get(exp[[x]]), x, basename(bam.file)) )
  
  names(bw_lst) <- names(exp[export])
  return(invisible( new("BigWigFileList", listData=bw_lst) ))
  
}
