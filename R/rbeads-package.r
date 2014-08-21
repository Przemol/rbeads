#' rbeads - R implementation of Bias Elimination Algorithm for Deep Sequencing
#'
#' BEADS is a normalization scheme that corrects nucleotide composition bias, mappability variations and differential local DNA structural effects in deep sequencing data. In high-throughput sequencing data, the recovery of sequenced DNA fragments is not uniform along the genome. In particular, GC-rich sequences are often over-represented and AT-rich sequences under-represented in sequencing data. In addition, the read mapping procedure also generates regional bias. Sequence reads that can be mapped to multiple sites in the genome are usually discarded. Genomic regions with high degeneracy therefore show lower mapped read coverage than unique portions of the genome. Mappability varies along the genome and thus creates systematic bias. Furthermore, local DNA or chromatin structural effects can lead to coverage inhomogeneity of sequencing data.
#' 
#' For detailed instruction how to run BEADS see \code{\link{beads}} function documentation.
#' 
#' @references 
#' \url{http://beads.sourceforge.net/} \cr 
#' \url{http://www.ncbi.nlm.nih.gov/pubmed/21646344}
#' 
#' @seealso 
#' \code{\link{beads}} \cr
#' \code{\link{sumBAMinputs}}
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
#' @author Przemyslaw Stempor
#' @import Rsamtools digest rtracklayer GenomicRanges BSgenome Biostrings IRanges methods
#' @keywords package
#' 
#' @docType package
#' @name rbeads-package
#' @aliases rbeads
NULL