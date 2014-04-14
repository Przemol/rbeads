library(devtools)
my_description <- list(
  Package='rBEADS',
  Type='Package',
  Title='R implementation of Bias Elimination Algorithm for Deep Sequencing',
  Version='0.3.0',
  Date=Sys.Date(),
  Author='Przemyslaw Stempor',
  Maintainer='Przemyslaw Stempor <ps562@cam.ac.uk>',
  Depends='R (>= 2.12.1), Rsamtools, BSgenome.Celegans.UCSC.ce10, digest, rtracklayer, GenomicRanges',
  Description='BEADS is a normalization scheme that corrects nucleotide composition bias, mappability variations and differential local DNA structural effects in deep sequencing data. In high-throughput sequencing data, the recovery of sequenced DNA fragments is not uniform along the genome. In particular, GC-rich sequences are often over-represented and AT-rich sequences under-represented in sequencing data. In addition, the read mapping procedure also generates regional bias. Sequence reads that can be mapped to multiple sites in the genome are usually discarded. Genomic regions with high degeneracy therefore show lower mapped read coverage than unique portions of the genome. Mappability varies along the genome and thus creates systematic bias. Futhermore, local DNA or chromatin structural effects can lead to coverage inhomogeneity of sequencing data.',
  License='LGPL',
  LazyLoad='yes'
  )
system('rm -Rf rbeads')
create('./rbeads', my_description, check=TRUE)
