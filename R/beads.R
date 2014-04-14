#' Main rBEADS wrapper
#'
#' This function should be started ni the BAM/Solexa data directory. 
#' It inputs CSV config file and run BEADS algoritm with specified settings.
#'
#' @param config Configuration file in CSV format.
#' @return LogFile - Read the log file for further run info.
#' 
#' @details
#'  rBEADS configuration file are 4 columns Comma Separated Values (CSV) tables. They can be created using any spreadsheet editor or plain text editor. The structure of sample input file:
#'  \preformatted{
#'    Sample,Control,rep,ER
#'    ChIPFile1.bam,InputFile1.bam,rep1,ERfile.bed
#'    ChIPFile2.bam,FRM,rep2a,none
#'    ChIPFile3.bam,EGS,rep2b,auto
#'  }
#' 
#'  The tabular structure:
#'  \tabular{llll}{
#'    Sample \tab Control \tab rep \tab ER \cr
#'    ChIPFile1.bam \tab InputFile1.bam \tab rep1 \tab ERfile.bed \cr
#'    ChIPFile2.bam \tab FRM \tab rep2a \tab none \cr
#'    ChIPFile3.bam \tab EGS \tab rep2b \tab auto \cr 
#'   }
#' 
#'  The \code{Sample} column contains ChIP experiment aligned short reads files.
#' 
#'  The \code{Control} column contains Input file corresponding to ChIP experiment. This filed may hold 2 special values - \code{EGS} forces BEADS wrapper to use pre-calculate EGS summed input track shipped with R package. Similarly \code{FRM} value indicates using summed formaldehyde input (also available in rBEADS package).
#' 
#'  \code{rep} column is not currently used by and have been incorporated to easier organize coning table. This filed can hold any value (excluding commas).
#' 
#'  The last column - \code{ER} - indicates which enriched regions should be used in BEADS run. The \code{auto} value indicates build in rBEADS algorithm for finding enriched regions, while \code{none} switches off enriched regions masking from the algorithm (useful for testing input files). Finally, if user decides to use external enriched regions finder (e. g. by peak calling) this field may take BED file name.
#' 
#' @references http://beads.sourceforge.net/
#' @author Przemyslaw Stempor
#' @keywords beads
#' @export
#' 
#' @examples
#' source("http://bioconductor.org/biocLite.R")
#' biocLite("BSgenome.Celegans.UCSC.ce10")
#' 
#' require(BSgenome.Celegans.UCSC.ce10)
#' require(rbeads)
#' beads(config='BeadsConfig.csv')


beads <-
function(config='BeadsConfig.csv') {
	
	data.dir <- getwd()
	setwd(data.dir)
	files <- read.csv(config)
	
	MAP <- data('MappabilityCe6', package='rBEADS')
	MAPF <- data('nonMappableFilterCe6', package='rBEADS')
	#FRMI <- load(file.path(source.dir, 'precalculated/SummedFormaldehydeInput_v2.Rdata'))
	
	for(i in 1:nrow(files))  {
		
		setwd(data.dir)
		
		##Set up naming
		sample.d <- reName(as.character(files$Sample[i]), 'dir', 'na', 'na', '', 'D', 0)
		
		##Set up 'uniq' parameter with backward compatibility to versions >=1.2.5
		if ( is.null(files$uniq[i]) ) {
			uniq <- FALSE
		} else {
			uniq <- as.logical(files$uniq[i])
		}
		
		##Log output
		dir.create(sample.d); setwd(sample.d)
		sink(file = 'LOG.txt', type = "output", split=TRUE)
		##zz <- file("LOG.ERROR.txt", open="wt"); sink(zz, split=TRUE); sink(zz, type="message")
		
		cat("/* File created on", date(), "*/\n")
		cat(i, ") Processing: ", as.character(files$Sample[i]), ' & ', as.character(files$Control[i]), "\n", sep='')
		
		if(files$Control[i] == 'FRM') {
			cat('INFO: Using summed formaldehyde input!\n')
			if(exists('summed.frm')) { control.map <- summed.frm } else { control.map <- get(data('summed.frm', package='rBEADS')) }
		} else if(files$Control[i] == 'EGS') {
			cat('INFO: Using summed EGS input!\n')
			if(exists('summed.egs')) { control.map <- summed.egs } else { control.map <- get(data('summed.egs', package='rBEADS')) }
		}  else if(grepl('Rdata$', as.character(files$Control[1]))) {
			cat('INFO: Using summed', as.character(files$Control[1]), 'input!\n')
			control.map <- get(load( file.path(data.dir, as.character(files$Control[1])) ))
		} else {
			##Control [INPUT]
			control.d <- reName(as.character(files$Control[i]), 'na', 'na', 'na', '', 'D', 0)
			if (grepl('bam$', as.character(files$Control[1]))) {
				control.re <- ImportBAM(bam.file=file.path(data.dir, as.character(files$Control[i])), uniq=uniq, desc=sample.d, resize_length=200, quality_cutoff=10, export_bin=F, export_track=F)
			} else if (grepl('export.fq', as.character(files$Control[1]))) {
				control.re <- ImportSolexa(file=as.character(files$Control[i]), data.dir=data.dir, desc=sample.d, resize_length=200, quality_cutoff=10, export_bin=FALSE, export_track=FALSE)
			} else {
				stop('Unknown input file type for sample!')
			}
			control.gc <- GCCorrection(control.re, enriched_regions=NULL, nonMappableFilter=get(MAPF), desc=control.d, smoothing_spline=FALSE)	
			control.map <- MappabilityCorrection(GCnormTrack=control.gc, mappabilityTrack=get(MAP))
		} 
		
		##Sample [ChIP]
		if (grepl('bam$', as.character(files$Sample[1]))) {
			sample.re <- ImportBAM(bam.file=file.path(data.dir, as.character(files$Sample[i])), uniq=uniq, desc=sample.d, resize_length=200, quality_cutoff=10, export_bin=F, export_track=F)
		} else if (grepl('export.fq', as.character(files$Sample[1]))) {
			sample.re <- ImportSolexa(file=as.character(files$Sample[i]), data.dir=data.dir, desc=sample.d, resize_length=200, quality_cutoff=10, export_bin=FALSE, export_track=FALSE)
		} else if (grepl('Rdata$', as.character(files$Sample[1]))) {
			sample.re <- get(load( file.path(data.dir, as.character(files$Sample[i])) ))
		} else {
			stop('Unknown input file type for sample!')
		}
		if(is.null(files$ER[i]) |  files$ER[i]=='auto' ){ 
			sample.er <- EnrichedRegions(sample.re, desc=sample.d)
		} else if (files$ER[i]=='no') { 
			cat('INFO: Skipping ER auto finder!\n')
			sample.er <- NULL
		} else {
			sample.er <- import.bed(files$ER[i], genome = "hg18")
		}
		sample.gc <- GCCorrection(sample.re, enriched_regions=sample.er, nonMappableFilter=get(MAPF), desc=sample.d, smoothing_spline=FALSE)
		sample.map <- MappabilityCorrection(GCnormTrack=sample.gc, mappabilityTrack=get(MAP))
		
		##Normalization
		sample.norm <- DivStep(sample.map, control.map, sample.er)
		
		
		#EXPORT	
		dir.create('Rbinaries'); setwd('Rbinaries')
		sample.re.cov <- coverage(sample.re)
		save(sample.re.cov, file=reName(sample.d, 'Raw', 'linar', '01bp', ext='.Rdata', 'B', 0))
		save(sample.gc, 	file=reName(sample.d, 'GCCorrected', 'linar', '01bp', ext='.Rdata', 'B', 0))
		save(sample.map, 	file=reName(sample.d, 'MapCorrected', 'linar', '01bp', ext='.Rdata', 'B', 0))
		save(sample.norm, 	file=reName(sample.d, 'NORM', 'linar', '01bp', ext='.Rdata', 'B', 0))
		setwd('..')
		
		if (exists('sample.re')) {
			dir.create('ChIP_QC'); setwd('ChIP_QC')
			BinTrack(coverage(sample.re), n=25, smooth=FALSE, out=reName(sample.d, 'Raw', 'linear', '25bp', ext='.wig'), type="WIG", name='RAW alignment', col='darkred')
			BinTrack(sample.gc, n=25, smooth=FALSE, out=reName(sample.d, 'GCCorrected', 'linear', '25bp', ext='.wig'), type="WIG", name='GC corrected', col='darkblue')
			BinTrack(sample.map, n=25, smooth=FALSE, out=reName(sample.d, 'MapCorrected', 'linear', '25bp', ext= '.wig'), type="WIG", name='Mappability corrected', col='blue')
			setwd('..')
		}
		if (exists('control.re')) {
			dir.create('Input_QC'); setwd('Input_QC')
			BinTrack(coverage(control.re), n=25, smooth=FALSE, out=reName(control.d, 'Raw', 'linear', '25bp', ext='.wig'), type="WIG", name='RAW alignment', col='black')
			BinTrack(control.gc, n=25, smooth=FALSE, out=reName(control.d, 'GCCorrected', 'linear', '25bp', ext='.wig'), type="WIG", name='GC corrected', col='darkblue')
			BinTrack(control.map, n=25, smooth=FALSE, out=reName(control.d, 'MapCorrected', 'linear', '25bp', ext='.wig'), type="WIG", name='Mappability corrected', col='blue')
			setwd('..')
		}
		dir.create('NORMALIZED'); setwd('NORMALIZED')
		BinTrack(sample.norm, n=25, smooth=FALSE, out=reName(sample.d, 'NORM', 'linear', '25bp', ext='.wig', 'N', 0), type="WIG", name=paste(sample.d, 'BEADS Normalised'), col='darkgreen')
		BinTrack(sample.norm, n=25, smooth=FALSE, out=reName(sample.d, 'NORM', 'Log2', '25bp', ext='.wig', 'N', 0), type="WIG", zscore=FALSE, log2=TRUE, name=paste(sample.d, 'BEADS Normalised & log2'), col='darkgreen')
		BinTrack(sample.norm, n=25, smooth=FALSE, out=reName(sample.d, 'NORM', 'Zscore', '25bp', ext='.wig', 'N', 0), type="WIG", zscore=TRUE, log2=FALSE, name=paste(sample.d, 'BEADS Normalised & z-scored'), col='darkgreen')
		BinTrack(sample.norm, n=25, smooth=FALSE, out=reName(sample.d, 'NORM', 'linear', '25bp', ext='.gff', 'N', 0), type="GFF")
		setwd('..')
		
		cat("/* Calculation complited on", date(), "*/\n")
		##sink(type="message"); sink(); close(zz); 
		sink(); setwd('..');
	}
}
