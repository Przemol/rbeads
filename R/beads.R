#' BEADS normaliztion algorithm
#'
#' This function should be started ni the BAM/Solexa data directory. 
#' It inputs CSV config file and run BEADS algoritm with specified settings.
#' http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability
#' http://www.plosone.org/article/fetchObject.action?uri=info%3Adoi%2F10.1371%2Fjournal.pone.0030377&representation=PDF
#' http://wiki.bits.vib.be/index.php/Create_a_mappability_track
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



setGeneric("beads",
  function(experiment, control, mappability, genome, uniq=TRUE, insert=200L, mapq_cutoff=10L, export='BEADS', rdata=FALSE, quickMap=TRUE, ...) 
  standardGeneric("beads")
)

setMethod("beads", signature(experiment='BamFile', control='BamFile', mappability='BigWigFile', genome='ANY'),
  function(experiment, control, mappability, genome, uniq=TRUE, insert=200L, mapq_cutoff=10L, export='BEADS', rdata=FALSE, quickMap=TRUE,...) {
    beads_bam_bam(path(experiment), path(control), mappability, genome, uniq, insert, mapq_cutoff, export, rdata, quickMap, ...)
  }
)

setMethod("beads", c(experiment='BamFile', control='BigWigFile', mappability='BigWigFile', genome='ANY'),
  function(experiment, control, mappability, genome, uniq=TRUE, insert=200L, mapq_cutoff=10L, export='BEADS', rdata=FALSE, quickMap=TRUE, ...) {
    beads_bam_bw(path(experiment), control, mappability, genome, uniq, insert, mapq_cutoff, export, rdata, quickMap, ...)
  }
)

setMethod("beads", signature(experiment='character', control='character', mappability='character', genome='character'),
  function(experiment, control, mappability, genome, uniq=TRUE, insert=200L, mapq_cutoff=10L, export='BEADS', rdata=FALSE, quickMap=TRUE, ...) {
    beads(FileForFormat(experiment), FileForFormat(control), FileForFormat(mappability), genome, uniq, insert, mapq_cutoff, export, rdata, quickMap, ...)
  }
)






