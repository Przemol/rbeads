BinTrack <-
function(divstep, n=25, smooth=FALSE, out="TEST.gff", type="WIG", zscore=FALSE, log2=FALSE, name=out, col='darkblue') {
	
	if(n > 1) {

		catTime("BIN: Preprocessing input [with smoothing=", smooth, "]", e={			
			if(smooth) {
				div.m <- divstep
				div.m[is.na(div.m)] <- 0
				div.m[div.m == Inf] <- 0
				div.m <- runmean(div.m, smooth, endrule="constant")
				div.m[is.na(divstep)] <- NA
				divstep <- div.m
				rm(div.m)			
			}
		})
		
		catTime("BIN: Binning @ ", n, "bp", e={				
			suppressWarnings( tiled <- genomeBlocks(Celegans, sort(seqnames(Celegans)[-7]), n) )
			if(names(divstep)[1] == "I") {
				seqlevels(tiled) <- c("I", "II", "III", "IV", "V", "X")
			}
			tiled.split <- split(start(tiled), as.character(seqnames(tiled)))
			values(tiled)$score <- unlist(lapply(names(divstep), function(x) as.numeric(divstep[[x]][tiled.split[[x]]])))		
		})

	} else {
		catTime("BIN: 1bp resolution, filtering only", e={		
					tiled <- as(as(divstep, "RangedData"), "GRanges")
		})
	}
	
	catTime("BIN: Filtering output", e={		
		tiled <- tiled[!is.na(values(tiled)$score)]
		tiled <- tiled[!is.nan(values(tiled)$score)]
		tiled <- tiled[!is.infinite(values(tiled)$score)]
		tiled <- tiled[!values(tiled)$score==0]
	})
	
	if(zscore) {
		catTime("BIN: Z-scoreing", e={		
			values(tiled)$score <- ( values(tiled)$score - mean( values(tiled)$score ) ) / sd( values(tiled)$score )	
		})
	} else if(log2) {
		catTime("BIN: log2 transforming", e={		
			values(tiled)$score <- log2(values(tiled)$score)	
		})
	}


	#seqlengths(Celegans)
	#names(seql) <- c("I", "II", "III", "IV", "V", "X", "MtDNA")
	#seqlengths(tiled) <- seqlengths(Celegans)
	catTime("BIN: Exporting output to file [", out, "]",  e={
		values(tiled)$score <- round(values(tiled)$score, 2)
		#seqlevels(tiled) <- c("I", "II", "III", "IV", "MtDNA", "V", "X")
		if(type == "GFF") { export.gff(as(tiled, "RangedData")[,2], out) }
		if(type == "BGR") { export.bed(as(tiled, "RangedData")[,2], out, variant = "bedGraph", color = NULL) }
		if(type == "WIG") { export.wig(as(tiled, "RangedData")[,2], out, dataFormat="variableStep", name=name, color=as.integer(col2rgb(col))) }
	})	
}
