beads_bamFiles <- function(bam.file, bam.control, mappability, uniq=TRUE, sample.d='', insert=200L, mapq_cutoff=10L, export='BEADS') {
  
  #Setting up mappability values
  gem  <- import.bw( mappability, as = "RleList" )
  MAP  <- round( runmean(gem, ifelse(insert %% 2 == 0, insert+1, insert), endrule = "constant"), 2 )
  MAPF <- ( gem > 0.5 )

  #Custom input processing
  control.re <-  ImportBAM(bam.control, uniq=uniq, resize_length=insert, quality_cutoff=mapq_cutoff)
  control.gc <-  GCCorrection(control.re, enriched_regions=NULL, nonMappableFilter=MAPF, desc=sample.d, smoothing_spline=FALSE)
  control.map <- MappabilityCorrection(GCnormTrack=control.gc, mappabilityTrack=MAP)
  
  #Full beads
  sample.re <-   ImportBAM(bam.file,    uniq=uniq, resize_length=insert, quality_cutoff==mapq_cutoff)
  sample.er <-   EnrichedRegions(sample.re, desc=sample.d)
  sample.gc <-   GCCorrection(sample.re, enriched_regions=sample.er, nonMappableFilter=MAPF, desc=sample.d, smoothing_spline=FALSE)
  sample.map <-  MappabilityCorrection(GCnormTrack=sample.gc, mappabilityTrack=MAP)
  
  sample.norm <- DivStep(sample.map, control.map, sample.er)
  
  exp <- list('control.re', 'control.gc', 'control.map', 'sample.coverage', 'sample.gc', 'sample.map', 'sample.norm')
  names(exp) <- c('control_readsCoverage', 'control_GCcorected', 'control.GCandMap', 'readsCoverage', 'GCcorected', 'GCandMap', 'BEADS')
  lapply( names(exp[export]), function(x) toBW(exp[[x]], x) )
}
