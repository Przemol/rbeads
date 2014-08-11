sumBAMinputs <- function(bam.controls=dir(pattern='bam$'), bw.mappability, 
                         genome, out_name='SummedInput', uniq=TRUE, insert=200L, 
                         mapq_cutoff=10L, export='BEADS', rdata=FALSE, 
                         export_er=TRUE, quickMap=TRUE, bin=25L) {

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
	  message(sprintf('Summing %s...', bam.controls[i]))
    coverage.i <- inputs.list[[i]] * (mean(nreads) / nreads[i])
		if (exists('summed.input')) {
			summed.input <- summed.input + coverage.i
		} else {
			summed.input <- coverage.i
		}	
	}
	summed.input <- binRle2Rle(summed.input, bin)
	export.bw(summed.input, paste0(out_name, '_binned25bp.wig') )

}
