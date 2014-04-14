ImportSolexa <-
function(file=dir(pattern="export.fq$"), data.dir=getwd(), desc=unlist(strsplit(file, "\\."))[1], resize_length=200, quality_cutoff=10, export_bin=FALSE, export_track=FALSE) {
	
	
	catTime("Reading Eland file:", e={	
				aln0 <- readAligned(data.dir, paste('^', file, '$', sep=''), "SolexaExport", filter=alignQualityFilter(threshold=10L))
				#sQ <- quality(alignQuality(aln0))
				#Q = 10 * log(1 + 10 ** (sQ / 10.0)) / log(10);
			})
	catTime("Processing alignment file:", e={	
				#Filer ouut reads of quality score lower than "quality_cutoff" [10]
				#aln2 <- (data.frame(aln2)[aln2[[1]]$mapq >= quality_cutoff,]);
				
				#Construct GRanges object
				#ranges.raw <- GRanges(seqnames = aln2$rname, ranges = IRanges(aln2$pos, width=aln2$qwidth), strand = aln2$strand);
				ranges.raw <- as(aln0, 'GRanges')
				ranges.raw <- ranges.raw[alignData(aln0)$filtering=='Y']
				rm(aln0)
				values(ranges.raw) <- NULL
				
				#Sort out chromosome naming and sequence lengths (ce6)
				seqlevels (ranges.raw) <- seqlevels (Celegans)[c(1,2,3,4,7,5,6)]
					#ovr <- findOverlaps(ranges.raw, GRanges(seqnames=seqnames(Celegans), IRanges(start=rep(1, length(seqlengths(Celegans))), end=seqlengths(Celegans)), seqlengths=seqlengths(Celegans)), type="within")
					overflow <- sapply(c("chrI", "chrII", "chrIII", "chrIV", "chrM", "chrV", "chrX"), function(x) which( end(ranges.raw) >= seqlengths(Celegans)[x]  & as.logical(seqnames(ranges.raw)==x) )  )
					if(length(unlist(overflow)) != 0) {	
						warning( sprintf("%i sequences found oudside reference Ce6 genome at chrs: %s!", length(unlist(overflow)), paste(names(sapply(overflow, length)[sapply(overflow, length)>0]), collapse=', ') )) 
						ranges.raw <- ranges.raw[-unlist(overflow)]					
					}
					seqlengths(ranges.raw) <- seqlengths(Celegans)[c(1,2,3,4,7,5,6)]
					#sapply(c("chrI", "chrII", "chrIII", "chrIV", "chrM", "chrV", "chrX"), function(x) sum( end(ranges.raw)[as.logical(seqnames(ranges.raw)==x)] >= seqlengths(Celegans)[x])  )
				
				
				
				#Resize sequences to 200bp //This resize method can be better with smooth end resizeing (as in peak calling)
				if( !is.null(resize_length) ) { ranges.raw <- resize(ranges.raw, resize_length) }
			})
	
	if (export_bin | export_track) { 
		catTime("Exporting Rdata and/or Wiggle file:\n", e={ 
					if (export_bin) {
						assign(sprintf("RawRanges.%s", desc), ranges.raw)
						save(list=sprintf("RawRanges.%s", desc), file=sprintf("RawRanges_%s.Rdata",  desc))
						rm(list=sprintf("RawRanges.%s", desc))
					}
					if (export_track) { 
						binTrack(coverage(ranges.raw), n=25, smooth=FALSE, out=sprintf("RAW_%s.wig", desc), type="WIG")
					}					
				})
	}
	return(ranges.raw)
}
