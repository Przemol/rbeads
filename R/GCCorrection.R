# Przemyslaw Stempor, 2014
##############################################################################
GCCorrection <- function(ranges.raw, enriched_regions, nonMappableFilter, REF, RL=200L, desc='', smoothing_spline=FALSE, maskLowHighGC=FALSE) {

	#Mask out reads in enriched regions
	if (!is.null(enriched_regions)) {
		catTime("Mask out reads in enriched regions", e={
			ERoverlaps <- findOverlaps(ranges.raw, enriched_regions, select="first") 	
		})
		#INFO: percentage of reads in enriched regions
		message(" INFO: percentage of reads in non-enriched regions: ", round(100 * sum(is.na(ERoverlaps))/length(ranges.raw), 2), "%" )
	}
 
	#Calculate input GC content from RL [200bp] extended reads (ranges.raw)
	catTime("Calculate input GCcontent", e={
      GCcontent <- unlist( lapply( split(ranges.raw, seqnames(ranges.raw)), function(x) {
        cat('.'); letterFrequency(REF[x], 'GC') 
      } ), use.names=FALSE)
		#GCcontent <- as.integer(letterFrequency(REF[ranges.raw], "GC"))
	})
	
	#Sample genome for GCcontent
	catTime("Sample genome for GCcontent", e={
		if (!is.null(enriched_regions)) {
			#Calculate logical vectors of non-enriched regions
			nonEnrichedRegionsLogi <- !coverage( enriched_regions )
			#Perform logical sum of non-enriched regions and mappable regions
			nonEnrichedMappableRegionsLogi <- nonEnrichedRegionsLogi & nonMappableFilter[names(nonEnrichedRegionsLogi)]
		} else {
			nonEnrichedMappableRegionsLogi <- nonMappableFilter
		}
		#Calculate GC pecrentage among the chromosomes
		if( any(sum(is.na(nonEnrichedMappableRegionsLogi))) ) nonEnrichedMappableRegionsLogi[ is.na(nonEnrichedMappableRegionsLogi) ] <- FALSE

		GCgenome <- Map(function(x, y) {
      cat('.'); 
      tabulate( letterFrequencyInSlidingView(x, RL, 'GC')[ as.logical(y) ], RL )
    }, REF[names(nonEnrichedMappableRegionsLogi)], nonEnrichedMappableRegionsLogi)
 
	})
	
	#INFO: percentage of genome to be sampled
	#message( "INFO: percentage of genome to be sampled: ", sum(as.numeric(elementLengths(GCgenome))) / sum(as.numeric(seqlengths( REF ))) ) 

	#Calculate histograms for genomic (a) nad and sample (b) GC content
	catTime("Calculate histograms for genomic (a) nad and sample (b) GC content", e={
		a <- rowSums(do.call(cbind, GCgenome))
		if (!is.null(enriched_regions)) {
			b <- tabulate(GCcontent[is.na(ERoverlaps)], RL)
		} else {
			b <- tabulate(GCcontent, RL)
		}
	
		pdf(sprintf("IMG - GCdistribution - %s.pdf", desc), width = 12.0, height = 7.5, onefile = FALSE, paper = "special", encoding = "TeXtext.enc")
    
		plot(  a/sum(a), col="blue", type="l", xlab='GCcontent', ylab='probablility', main=sprintf("GCdistribution - %s", desc), lwd=2)
		lines( b/sum(b), col="red",  type="l", lwd=2)
		legend("topleft", c("Non enriched GENOMIC GC content distribution", "Non enriched SAMPLE GC content distribution"), fill=c("blue", "red") )
		abline(h=0, v=range(which( b/sum(b) > 0.0001 )), col = "gray60", lwd=1.5)
    dev.off()
	})
	
	#Calculate GC weighting vector
	catTime("Calculate GC weighting vector", e={												
		scales <- (a/sum(a)) / (b/sum(b))
		scales[is.infinite(scales)] <- NA
	
		pdf(file=sprintf("IMG - GCoutlier - %s.pdf", desc), width = 12.0, height = 7.5, onefile = FALSE, paper = "special", encoding = "TeXtext.enc")
		
    plot(scales, ylim=c(0,5), main=sprintf("GCoutlier - %s", desc), pch=20, col='grey', cex=1)
		abline(h=0, v=range(which( b/sum(b) > 0.0001 )), col = "gray60", lwd=2)
    
		scales[ b/sum(b) <= 0.0001 ] <- NA
    
		SSpoints <- cbind(1:RL, scales)[!is.na(scales) & !(scales %in% boxplot(scales, plot=F)$out), ]	
		SS <- smooth.spline(SSpoints, df=10)
		scales.fit <- cbind(predict(SS, 0:RL)$x, predict(SS, 0:RL)$y)
		lines(scales.fit, col="darkgreen")
		points(SSpoints, col="darkgreen", pch="o")
		if (smoothing_spline) { scales <- scales.fit }
    
		points(scales, pch=20, col='black')
		
		dev.off()
	})
	
	#Assign GC score to every read
	catTime("Assign GC score to every read", e={
		GCcontent[GCcontent == 0] <- 1
		scores <- scales[GCcontent]
	})
	
	#Prepare GCscore-enriched GRanges
	catTime("Prepare GCscore-enriched GRanges", e={				
		ranges.raw$value <- scores
		ranges.f1 <- ranges.raw[ !is.na(ranges.raw$value) ]
	})

	#calculate GC weighted coverage with given accuracy
	catTime("Calculate GC weighted coverage", e={
		#Output coverage
		cov <- coverage(ranges.f1, weight='value')
		cov <- round(cov, 3)
	})
 
	#Mask nonGC correctable regions
  if( maskLowHighGC ) {
	  catTime("Masking non-GCcorrectable regions", e={
      rang <- as(seqinfo(ranges.f1), 'GRanges')
      cutoff <- range(which( b/sum(b) > 0.0001 ))
		  notGCcorrectableRegions2 <- sapply( split(rang, seqnames(rang)) , function(x) {
			  GCchr <- letterFrequencyInSlidingView(REF[x][[1]], RL, "GC"); cat('.')
        ind <- which( ! (GCchr >= cutoff[1] & GCchr <= cutoff[2]) ); cat('.')
        if( length(ind) ) {
			  GRanges( seqnames=seqnames(x), ranges=reduce(IRanges( ind , width=RL)) )
        } else { NULL }
		  })
		  notGCcorrectableRegions <- GRangesList(notGCcorrectableRegions2[elementLengths(notGCcorrectableRegions2)>0 ])
      cov[ coverage(notGCcorrectableRegions)[names(cov)] > 0 ] <- NA
	  })
	  message("INFO: masking ", round(sum(as.numeric(sum( coverage(notGCcorrectableRegions)[names(cov)] > 0 ))) / sum(as.numeric(seqlengths(rang)))*100, 2), "% of genome in gc-correctable regions")
  }

	return(cov)
}
