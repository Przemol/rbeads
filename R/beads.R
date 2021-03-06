#' BEADS normalization algorithm for ChIP-seq experiments 
#' 
#' The function runs BEAD algorithm on aligned files. It requires short read alignment 
#' in \href{http://genome.ucsc.edu/goldenPath/help/bam.html}{BAM format}, 
#' input track (either single experiment \href{http://genome.ucsc.edu/goldenPath/help/bam.html}{BAM}
#' or summed input created using \code{\link{sumBAMinputs}} function), mappability track in 
#' \href{http://genome.ucsc.edu/goldenPath/help/bigWig.html}{BigWig format} (see details) 
#' and reference genome in \href{http://en.wikipedia.org/wiki/FASTA_format}{FASTA fomat} 
#' or as \code{\link[BSgenome]{BSgenome}} genomic package.
#'
#' @param experiment The path to experiment (sample) alignment file in BAM format or \code{\link[Rsamtools]{BamFile}} class.
#' @param control The path to control (input) alignment file in BAM format or summed input file in BigWiggle format (accepts \code{\link[Rsamtools]{BamFile}} and \code{\link[rtracklayer]{BigWigFile}} classes as well).
#' @param mappability The path to mappability track in BigWiggle format (accepts \code{\link[rtracklayer]{BigWigFile}} class as well).
#' @param genome The path reference genome FASTA or UCSC identifier for installed \code{\link[BSgenome]{BSgenome}} packages e.g. "hg19" for human
#' @param uniq If TRUE the alignment will be uniqued, i.e. only one of non-unique reads will be used.
#' @param insert The expected insert size in base pairs.
#' @param mapq_cutoff The cutoff parameter used to filter BAM alignments for low mapping quality reads.
#' @param export The character vector of BW tracks to be exported.
#' @param rdata If TRUE all data will be exported as R binaries in addition to BigWiggle tracks.
#' @param export_er If TRUE the enriched regions will be exported to BED format.
#' @param quickMap If TRUE the quick mappability processing be used, otherwise the mappability track will be processed by running mean smoothing.
#' @param ... parameters passed by reference to rbeads internal functions.
#'  
#' @return Named \code{\link[rtracklayer]{BigWigFileList}} containing exported BW file connections.
#' 
#' @details
#' Mappability/alignability tracks gives numeric score for level of reference sequence uniqueness.
#' Short reads cannot be confidently aligned to non-unique sequences, so BEADS masks the out.
#' The pre-calculated tracks for many species can be found in genome databases, e.g. following link 
#' gives the collection of human tracks for reference genome GRCh37/hg19:
#' \url{http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=340327143&g=wgEncodeMapability}.
#' 
#' For other species mappability track can be easily calculated from reference FASTA file
#' using GEM-mappability software (\url{http://www.ncbi.nlm.nih.gov/pubmed/22276185}). The 
#' GEM library binaries are available on SourceForge: (\url{http://sourceforge.net/projects/gemlibrary/files/gem-library/}),
#' while the reference genome FASTA file can be obtained from UCSC: (\url{http://hgdownload.soe.ucsc.edu/downloads.html}).
#' The following example illustrates the procedure for \emph{C. elegans} reference genome
#' (36bp read length and 8 parallel threads set in options):
#'
#' \code{gem-indexer -i ce10.fa -o ce10 -T 8}\cr
#' \code{gem-mappability -I ce10.gem -l 36 -o ce10 -T 8}\cr
#' \code{gem-2-wig -I ce10.gem -i ce10.map -o ce10}\cr
#' \code{wigToBigWig ce10.wig ce10.chrom.sizes ce10.bw}
#' 
#' By default only BEADS normalized track is exported. The export files can be 
#' control by \code{export} parameter,  which is the character vector containing 
#' following values: 
#' \describe{
#'  \item{\code{'BEADS'}}{fully normalized files, i.e. GC correction, mappability correction, division by input and scaling by median}
#'  \item{\code{'GCandMap'}}{GC correction and mappability correction}
#'  \item{\code{'GCcorected'}}{GC correction}
#'  \item{\code{'readsCoverage'}}{raw reads coverage, after extending to \code{insert} 
#'    length, filtering /code{mapq_cutoff} and uniquing if \code{uniq}}
#'  \item{\code{'control_GCandMap'}}{GC normalized and mappability corrected input, only works if input is a BAM file}
#'  \item{\code{'control_GCcorected'}}{GC normalized, only works if input is a BAM file}
#'  \item{\code{'control_readsCoverage'}}{input raw reads coverage, after extending to
#'    \code{insert} length, filtering \code{mapq_cutoff} and uniquing if 
#'    \code{uniq}, only works if input is a BAM file}
#' }
#' For example: \code{  beads(sample_bam, input_bam, map, fa, export=c('control_readsCoverage', 'control_GCcorected', 'control_GCandMap', 'readsCoverage', 'GCcorected', 'GCandMap', 'BEADS'))}
#' 
#' @references \url{http://beads.sourceforge.net/} \cr \url{http://www.ncbi.nlm.nih.gov/pubmed/21646344}
#' 
#' @author Przemyslaw Stempor
#' @export
#' 
#' @examples
#' # Get the paths of example files
#' sample_bam <- system.file("extdata", "GSM1208360_chrI_100Kb_q5_sample.bam", package="rbeads")
#' input_bam <- system.file("extdata", "Input_fE3_AA169.bam", package="rbeads")
#' SummedInput_bw <- system.file("extdata", "Ce10_HiSeqFRMInput_UNIQ_bin25bp_chrI_100Kb_sample.bw", package="rbeads")
#' map_bw <- system.file("extdata", "ce10_mappability_chrI_100Kb_sample.bw", package="rbeads")
#' ref_fa <- system.file("extdata", "ce10_chrI_100Kb_sample.fa", package="rbeads")
#' 
#' # Set the directory where the output files will be crated
#' setwd(tempdir())
#' 
#' # Run BEADS for BAM input file
#' beads(sample_bam, input_bam, map_bw, ref_fa)
#' 
#' # Run BEADS for SummedInput (BigWig) input file
#' beads(sample_bam, SummedInput_bw, map_bw, ref_fa)
#' 
#' \dontrun{
#' ## Run BEADS for BSgenome package, the reference genome package have to be installed prior to running this example
#' # source("http://bioconductor.org/biocLite.R")
#' # biocLite("BSgenome.Celegans.UCSC.ce10")
#' # library(BSgenome.Celegans.UCSC.ce10)
#' beads(sample_bam, SummedInput_bw, map_bw, genome='ce10')
#' 
#' ## Run BEADS for all BAM files in the directory
#' #lapply(dir(pattern='bam$'), beads, control=input, mappability=map_bw, genome=ref_fa)
#' }
#' 
#' \dontrun{
#' ## In parallel:
#' # library(parallel)
#' # mclapply(dir(pattern='.bam$'), beads, control=input, mappability=map_bw, genome='ce10', mc.cores=parallel::detectCores())
#' }
setGeneric("beads",
  function(experiment, control, mappability, genome, uniq=TRUE, insert=200L, mapq_cutoff=10L, export='BEADS', rdata=FALSE, export_er=TRUE, quickMap=TRUE, subsample=0L, ...) 
  standardGeneric("beads")
)

#' @describeIn beads Method for signature experiment='BamFile', control='BamFile', mappability='BigWigFile', genome='ANY'
setMethod("beads", signature(experiment='BamFile', control='BamFile', mappability='ANY', genome='ANY'),
  function(experiment, control, mappability, genome, uniq=TRUE, insert=200L, mapq_cutoff=10L, export='BEADS', rdata=FALSE, quickMap=TRUE, subsample=0L, ...) {
    beads_bam_bam(path(experiment), path(control), mappability, genome, uniq, insert, mapq_cutoff, export, rdata, quickMap, subsample, ...)
  }
)

#' @describeIn beads Method for signature experiment='BamFile', control='BigWigFile', mappability='BigWigFile', genome='ANY'
setMethod("beads", c(experiment='BamFile', control='BigWigFile', mappability='ANY', genome='ANY'),
  function(experiment, control, mappability, genome, uniq=TRUE, insert=200L, mapq_cutoff=10L, export='BEADS', rdata=FALSE, quickMap=TRUE, subsample=0L, ...) {
    beads_bam_bw(path(experiment), control, mappability, genome, uniq, insert, mapq_cutoff, export, rdata, quickMap, subsample, ...)
  }
)

#' @describeIn beads Method for signature experiment='character', control='character', mappability='character', genome='character'
setMethod("beads", signature(experiment='character', control='character', mappability='ANY', genome='character'),
  function(experiment, control, mappability, genome, uniq=TRUE, insert=200L, mapq_cutoff=10L, export='BEADS', rdata=FALSE, quickMap=TRUE, subsample=0L, ...) {
    beads(FileForFormat(experiment), FileForFormat(control), mappability, genome, uniq, insert, mapq_cutoff, export, rdata, quickMap, subsample, ...)
  }
)

.onAttach <- function(libname, pkgname) {
  packageStartupMessage('The rbeads package ready, type "?beads" for help.')
}






