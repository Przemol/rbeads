# Przemyslaw Stempor, 2014
##############################################################################
ImportBAM <- function(bam.file=dir(pattern="\\.bam$")[1], REF, uniq=FALSE, resize_length=200, quality_cutoff=10) {
	
	#Read sam allignment file
	catTime("Reading alignment file [", bam.file, "]",  e={
		what <- c("rname", "strand", "pos",  "qwidth", "mapq")
		flag <- scanBamFlag(isUnmappedQuery = FALSE)
		param <- ScanBamParam(flag = flag, simpleCigar = FALSE, what = what)
		
		aln2 <- scanBam(bam.file, param=param)
	})
	
	catTime("Processing alignment file [uniq=", uniq,"]:", e={	
				
		#Filer ouut reads of quality score lower than "quality_cutoff" [10]
		lg <- ( aln2[[1]]$mapq >= quality_cutoff ) 
		 
		#Construct GRanges object
		ranges.raw <- GRanges(seqnames = aln2[[1]]$rname[lg], ranges = IRanges(aln2[[1]]$pos[lg], width=aln2[[1]]$qwidth[lg]), strand = aln2[[1]]$strand[lg]);
		suppressWarnings({
		  if( class(seqinfo(BamFile(bam.file))) == "Seqinfo") {
		    seqinfo(ranges.raw) <- seqinfo(BamFile(bam.file))[seqlevels(ranges.raw)]
      } else {
        seqlengths(ranges.raw) <- seqlengths(REF)[seqlevels(ranges.raw)]
      }
    })
    
		#Determinig if the ranges should be unified
		if (uniq == TRUE) { 
		  suppressWarnings( ranges.raw <- unique(ranges.raw) )
		}
			
		#Resize sequences to 200bp //This resize method can be better with smooth end resizeing (as in peak calling)
		if( !is.null(resize_length) ) { 
		  suppressWarnings( ranges.raw <- resize(ranges.raw, resize_length) )
		}
		
		#Calculate oryginal BAM file statistics
		num <- countBam(bam.file)
  
	})
	
  if(uniq) {
    message(sprintf("INFO: %0.0f out of %0.0f (%0.2f%%) reads mapped uniquely with %0.0f quality cutoff: [All=%0.0f, Mapped=%0.0f, Qmapped=%0.0f, Quniq=%0.0f].", 
                    length(ranges.raw), num$records, 100*length(ranges.raw)/num$records, quality_cutoff, num$records, length(lg), sum(lg), length(ranges.raw) ))
  } else {
    message(sprintf("INFO: %0.0f out of %0.0f [%0.2f%%] reads mapped with %0.0f quality cutoff: [All=%0.0f, Mapped=%0.0f, Qmapped=%0.0f].", 
                    length(ranges.raw), num$records, 100*length(ranges.raw)/num$records, quality_cutoff, num$records, length(lg), sum(lg) ))
  }
	
  return(ranges.raw)

}