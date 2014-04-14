
# Version log:
#
# 1.3.1: 	'uniq' parameter determining if the ranges should be unified prior output
# 						
# 
#frame_files <- lapply(sys.frames(), function(x) x$ofile)
#frame_files <- Filter(Negate(is.null), frame_files)
#PATH <- dirname(frame_files[[length(frame_files)]])

# Author: przemol
###############################################################################

require(Rsamtools)
require(BSgenome.Celegans.UCSC.ce10)
#require(BSgenome.Celegans.UCSC.ce6)

ImportBAM <- function(bam.file=dir(pattern="\\.bam$")[1], desc=unlist(strsplit(bam.file, "\\."))[1], uniq=FALSE, resize_length=200, quality_cutoff=10, export_bin=TRUE, export_track=TRUE) {
	
	
	#Read sam allignment file
	catTime("Reading alignment file [", bam.file, "]",  e={
		what <- c("rname", "strand", "pos",  "qwidth", "mapq")
		flag <- scanBamFlag(isUnmappedQuery = F)
		param <- ScanBamParam(flag = flag, simpleCigar = F, what = what)
		
		aln2 <- scanBam(bam.file, param=param)
	})
	
	catTime("Processing alignment file [uniq=", uniq,"]:", e={	
				
			#Filer ouut reads of quality score lower than "quality_cutoff" [10]
		lg <- (aln2[[1]]$mapq >= quality_cutoff) & (!grepl('M', aln2[[1]]$rname ))
		
			#Construct GRanges object
		ranges.raw <- GRanges(seqnames = aln2[[1]]$rname[lg], ranges = IRanges(aln2[[1]]$pos[lg], width=aln2[[1]]$qwidth[lg]), strand = aln2[[1]]$strand[lg]);
		
			#Determinig if the ranges should be unified
		if (uniq == TRUE) { 
			#unique.gr <- getMethod('unique', signature='GenomicRanges', where='GenomicRanges')
			ranges.raw <- unique(ranges.raw)
		}
		
			#Sort out chromosome naming and sequence lengths (ce6)
		seqlevels (ranges.raw) <- seqlevels(ranges.raw)[!grepl('M', seqlevels(ranges.raw))]
		seqlevels (ranges.raw) <- seqlevels (Celegans)[c(1,2,3,4,5,6)]
		seqlengths(ranges.raw) <- seqlengths(Celegans)[c(1,2,3,4,5,6)]
			
			#Resize sequences to 200bp //This resize method can be better with smooth end resizeing (as in peak calling)
		if( !is.null(resize_length) ) { ranges.raw <- resize(ranges.raw, resize_length) }
		
			#Calculate oryginal BAM file statistics
		num <- countBam(bam.file)
	})
	
	cat(sprintf("\tINFO: %0.0f out of %0.0f [%0.2f%%] reads mapped with %0.0f quality cutoff (%0.0f out of %0.0f [%0.2f%%] nucleotides). \n", length(ranges.raw), num$records, 100*length(ranges.raw)/num$records, quality_cutoff, mean(aln2[[1]]$qwidth[lg])*length(ranges.raw), num$nucleotides, 100*(mean(aln2[[1]]$qwidth[lg])*length(ranges.raw))/num$nucleotides))
	
	rm(aln2)
	
	if (export_bin | export_track) { 
			catTime("Exporting Rdata and/or Wiggle file:\n", e={ 
			if (export_bin) {
				assign(sprintf("RawRanges.%s", desc), ranges.raw)
				save(list=sprintf("RawRanges.%s", desc), file=sprintf("RawRanges_%s.Rdata",  desc))
				rm(list=sprintf("RawRanges.%s", desc))
			}
			if (export_track) { 
				BinTrack(coverage(ranges.raw), n=25, smooth=FALSE, out=sprintf("RAW_%s.wig", desc), type="WIG")
			}					
		})
	}
	return(ranges.raw)
}


catTime <- function(..., e=NULL, file="", gc=FALSE) {
	cat(..., "...", sep="", file=file, append=TRUE)
	cat("\t<", system.time(e, gcFirst=gc)[3], "s>\n", sep="", file=file, append=TRUE)	
}