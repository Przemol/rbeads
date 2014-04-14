#' @docType package
#' @name rbeads
#' 
#' @title rbeads - R implementation of Bias Elimination Algorithm for Deep Sequencing
#' @description
#' BEADS is a normalization scheme that corrects nucleotide composition bias, mappability variations and differential local DNA structural effects in deep sequencing data. In high-throughput sequencing data, the recovery of sequenced DNA fragments is not uniform along the genome. In particular, GC-rich sequences are often over-represented and AT-rich sequences under-represented in sequencing data. In addition, the read mapping procedure also generates regional bias. Sequence reads that can be mapped to multiple sites in the genome are usually discarded. Genomic regions with high degeneracy therefore show lower mapped read coverage than unique portions of the genome. Mappability varies along the genome and thus creates systematic bias. Futhermore, local DNA or chromatin structural effects can lead to coverage inhomogeneity of sequencing data.
#' 
#' @references http://beads.sourceforge.net/manual.php
#' @author Przemyslaw Stempor
#' @details
#' \tabular{ll}{
#'  Package: \tab rBEADS\cr
#'  Type: \tab Package\cr
#'  Version: \tab 1.0\cr
#'  Date: \tab 2011-07-04\cr
#'  License: \tab LGPL\cr
#'  LazyLoad: \tab yes\cr
#' }
#' <<echo=FALSE>>=
#' <<double>>
#'   rnorm(1)
#' 
#' 
#' @keywords CHiP-seq beads
#' @examples example(beads)
#' 
#' @import IRanges Rsamtools digest rtracklayer GenomicRanges
#' @exportPattern ^[[:alpha:]]+
NULL