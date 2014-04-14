SumBAMinputs <-
function(files=dir(pattern='bam$'), uniq=FALSE, export.all=FALSE, qc=FALSE, out.name='TestSummedInput_v01') {

	MAP <- data('MappabilityCe6', package='rBEADS')
	MAPF <- data('nonMappableFilterCe6', package='rBEADS')
	
	inputs.list <- list()
	nreads <- NULL
	
	for(i in 1:length(files))  {
		
		##Set up naming
		control.d <- paste(unlist(strsplit(files[i], "\\."))[-length(unlist(strsplit(files[i], "\\.")))], collapse='_')
		
		if(export.all) {
			##Log output
			dir.create(control.d); setwd(control.d)
			sink(file = 'LOG.txt', type = "output", split=TRUE)
			cat("/* File created on", date(), "*/\n")
			fp <- file.path('..', files[i])
		} else {
			fp <- file.path('.', files[i])
		}
		
		cat(i, ") Processing: ", as.character(files[i]), "\n", sep='')
		
		##Input BAM file, process GC correction and mappability mapping
		control.re <- ImportBAM(fp, uniq=uniq, desc=control.d, resize_length=200, quality_cutoff=10, export_bin=FALSE, export_track=FALSE)
		control.gc <- GCCorrection(control.re, enriched_regions=NULL, nonMappableFilter=get(MAPF), desc=control.d, smoothing_spline=FALSE)	
		control.map <- MappabilityCorrection(GCnormTrack=control.gc, mappabilityTrack=get(MAP))
		inputs.list[control.d] <- control.map
		
		nreads[i] 			<- length(control.re)
		
		if(export.all) {
			#EXPORT debug tracks
			BinTrack(coverage(control.re), n=25, smooth=FALSE, out=reName(control.d, 'Raw', 'linear', '25bp', '.wig'), type="WIG", name='RAW alignment', col='black')
			BinTrack(control.map, n=25, smooth=FALSE, out=reName(control.d, 'GCmap', 'linear','25bp', '.wig'), type="WIG", name='Mappability corrected', col='blue')
			
			assign(sprintf("INPUT_%s", control.d), control.map)
			save(   list=sprintf("INPUT_%s", control.d), file=reName(control.d, 'GCmap', 'linear', '01bp', '.Rdata'))
			
			cat("/* Calculation complited on", date(), "*/\n")
			sink(); setwd('..');
		}
		gc()
	}
	
	# Sum the input tracks (scaled by number of reads) 
	suppressWarnings(rm(summed.input))
	for(i in seq(along = nreads)) {
		cat(sprintf('Processing %s...\n', names(inputs.list[i])))
		coverage.i <- inputs.list[[i]] * (mean(nreads) / nreads[i])
		
		if (exists('summed.input')) {
			summed.input <- summed.input + coverage.i
		} else {
			summed.input <- coverage.i
		}	
	}
	summed.input <- round(summed.input, 0)
	save(summed.input, file=paste(out.name,'Rdata',sep='.'))
	BinTrack(summed.input, n=25, smooth=FALSE, out=paste(out.name,'@25bp.wig',sep=''), type="WIG", name=out.name, col='gray')
	
	#QC
	if ( qc ) {
		cat('Building quality control report...', '\n')
		library(ShortRead)
		qas <- list()
		for (i in files) {
			cat(i, '\n')
			qas[i] <- qa(getwd(), pattern=gsub('\\^', '.', i), type="BAM")
			gc()
		}
		qa <- do.call(rbind, qas)
		report(qa, dest='QAreport', type="html")		
	}

}
