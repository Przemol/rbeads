#' Creates summed input track from multiple ChIP-seq input experiments 
#' 
#' The function creates the \href{http://genome.ucsc.edu/goldenPath/help/bigWig.html}{BigWig formatted}
#' summed input from aligned input experiments in \href{http://genome.ucsc.edu/goldenPath/help/bam.html}{BAM format}.
#' The first two step of BEADS normalization, i.e. GC normalization and mappability correction are performed 
#' independently on all input experiments. Further, the obtain signals are summed. The tracks are normalized for
#' tag count, so each experiment contribute in equally to summed input. It is important to select only good quality run 
#' for summed input preparation. \href{http://www.bioinformatics.babraham.ac.uk/projects/fastqc/}{FastQC}
#' is a good program for sequencing experiment quality assessment. It is important to create the summed input with the 
#' same parameters as will be used for future BEADS runs.
#'
#' @param bam.controls The character vector containing paths to control (input) alignment files in BAM format
#' @param bw.mappability The path to mappability track in BigWiggle format (accepts \code{\link[rtracklayer]{BigWigFile}} cless as well)
#' @param genome The path reffrence genome FASTA or UCSC identifier fo installed \code{\link[BSgenome]{BSgenome}} packages e.g. "hg19" for human
#' @param out_name The prefix for exported BigWiggle, full name \code{out_name}_SummedInput_linear_\code{bin}bp.bw
#' @param uniq If TRUE the alignemnt will be uniqued, i.e. only one of non-unique reads will be used
#' @param insert The expected insert size in base pairs.
#' @param mapq_cutoff The cutoff parameter used to filter BAM alignments for low mapping quality reads.
#' @param quickMap If TRUE the quick mappability processing be used, otherwise the mappability track will be processed by running mean smoothing 
#' @param bin The desired binning window. 1L disables binning. It is recommended 
#' to set it >1L, as it greatly reduces the size of the output BW files and speeds 
#' up further BEADS runs.
#' 
#' @return \code{\link[rtracklayer]{BigWigFile}} class containing connection to summed input BigWig file.
#' 
#' @references \url{http://beads.sourceforge.net/} \cr \url{http://www.ncbi.nlm.nih.gov/pubmed/21646344}
#' 
#' @author Przemyslaw Stempor
#' @export
#' 
#' @examples
#' # Get the paths of example files
#' input_bam <- system.file("extdata", "Input_fE3_AA169.bam", package="rbeads")
#' map_bw <- system.file("extdata", "ce10_mappability_chrI_100Kb_sample.bw", package="rbeads")
#' ref_fa <- system.file("extdata", "ce10_chrI_100Kb_sample.fa", package="rbeads")
#' 
#' # Set the directory where the output files will be crated
#' setwd(tempdir())
#' 
#' # Calculate summed input; in this example we will pretend that we have 3 test input files,
#' # in real applications these should be different experiments
#' out1 <- sumBAMinputs(c(input_bam, input_bam, input_bam), map_bw, ref_fa)
#'
#' # Sanity check: The mean signal ratio between previous example (using 3 identical inputs) 
#' # and single normalized input should be very close to 3 
#' # (limited by binning function and BigWig summary numerical precision) 
#' out2 <- sumBAMinputs(input_bam, map_bw, ref_fa)
#' rtracklayer::summary(out1)[[1]]$score / rtracklayer::summary(out2)[[1]]$score
#'
sumBAMinputs <- function(bam.controls=dir(pattern='bam$'), bw.mappability, 
                         genome, out_name=paste0('Summed_', length(bam.controls), '_experiments'), uniq=TRUE, insert=200L, 
                         mapq_cutoff=10L, quickMap=TRUE, bin=25L) {

  #Importing refference genome
  message('Importing refference genome...')
  REF <- getREF(genome)
  
  #Setting up mappability values
  if( quickMap ) {
    message('Importing mappability [binning]...')
    if( class(bw.mappability) == 'BigWigFile' ) { bwf <- bw.mappability } else { bwf <- BigWigFile( bw.mappability ) }
    MAP  <- IRanges::summary(bwf, size=seqlengths(bwf)/insert, asRle=TRUE)
    MAP[is.na(MAP)] <- 0
    MAPF <- ( MAP > 0.5 )
  } else {
    message('Importing mappability [running median smoothing]...')
    gem  <- import.bw( bw.mappability, as = "RleList" )
    MAP  <- round( runmean(gem, ifelse(insert %% 2 == 0, insert+1, insert), endrule = "constant"), 2 )
    MAP[is.na(MAP)] <- 0
    MAPF <- ( gem > 0.5 )
  }
	
	inputs.list <- list()
	nreads <- NULL
	
	for(i in 1:length(bam.controls))  {
		
	  bam.control <- bam.controls[i]
	  
    message(i, ") Processing: ", as.character(bam.control))
	  #Loading input from BigWig file
	  message('Importing input...')
	  control.re <-  ImportBAM(bam.control, REF=REF, uniq=uniq, resize_length=insert, quality_cutoff=mapq_cutoff)
	  control.gc <-  GCCorrection(control.re, enriched_regions=NULL, REF=REF, nonMappableFilter=MAPF, RL=insert, desc=gsub('.bam$', '', basename(bam.control)), smoothing_spline=FALSE)
	  #message('Smoothing...')
	  #control.gc <-  round( runmean(control.gc, ifelse(insert %% 2 == 0, insert+1, insert), endrule = "constant"), 2 )
	  control.map <- MappabilityCorrection(GCnormTrack=control.gc, mappabilityTrack=MAP)
	
		inputs.list[[i]] <- control.map
		nreads[i] <- length(control.re)	
	}
	
	# Sum the input tracks (scaled by number of reads) 
	for(i in seq(along = nreads)) {
	  message(sprintf('Summing %s', bam.controls[i]), ' [weight=', (mean(nreads) / nreads[i]), ']...')
    coverage.i <- inputs.list[[i]] * (mean(nreads) / nreads[i])
		if (exists('summed.input')) {
			summed.input <- summed.input + coverage.i
		} else {
			summed.input <- coverage.i
		}	
	}
	summed.input <- binRle2Rle(summed.input, bin)
  
  bw <- toBW_missing(summed.input, 'SummedInput', out_name, resolution=bin)
  message('\nSummed inpud exported to BigWiggle file: ',  path(bw) )
  
  return(invisible(bw))

}
