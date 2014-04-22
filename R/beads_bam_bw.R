beads_bam_bw <- function(bam.file, bw.control, bw.mappability, genome_fasta=NULL, genome_package=NULL, uniq=TRUE, sample.d='', insert=200L, mapq_cutoff=10L, export='BEADS') {
  
  #Importing refference genome
  message('Importing refference genome...')
  if( !is.null(genome_fasta) ) {
    REF <- readDNAStringSet( genome_fasta )
  } else {
    package <- grep(genome_package, installed.genomes(), value=TRUE, ignore.case=TRUE)
    if( !length(package) ) { stop('Genome ', genome_package, ' is not installed. Please run "available.genomes()" to get genomes supported by BioConductor or provider reference FASTE file.', call.=FALSE) }
    library(package, character.only = TRUE)
    REF <- getSeq( get(package) ); names(REF) <- seqnames(get(package))
  }
  
  #Setting up mappability values
  message('Importing mappability...')
  gem  <- import.bw( bw.mappability, as = "RleList" )
  MAP  <- round( runmean(gem, ifelse(insert %% 2 == 0, insert+1, insert), endrule = "constant"), 2 )
  MAPF <- ( gem > 0.5 )
  
  #Loading input from BigWig file
  message('Importing input...')
  control.map <- import.bw(bw.control, as="RleList")
  
  #Full beads
  sample.re <-   ImportBAM(bam.file, REF=REF, uniq=uniq, resize_length=insert, quality_cutoff=mapq_cutoff)
  sample.er <-   EnrichedRegions(sample.re, desc=sample.d)
  sample.gc <-   GCCorrection(sample.re, enriched_regions=sample.er, REF=REF, nonMappableFilter=MAPF, RL=insert, desc=sample.d, smoothing_spline=FALSE)
  sample.map <-  MappabilityCorrection(GCnormTrack=sample.gc, mappabilityTrack=MAP)
  
  sample.norm <- DivStep(sample.map, control.map, sample.er)
  
  exp <- list('control.re', 'control.gc', 'control.map', 'sample.coverage', 'sample.gc', 'sample.map', 'sample.norm')
  names(exp) <- c('control_readsCoverage', 'control_GCcorected', 'control.GCandMap', 'readsCoverage', 'GCcorected', 'GCandMap', 'BEADS')
  lapply( names(exp[export]), function(x) toBW(exp[[x]], x) )
}
