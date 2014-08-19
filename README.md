rBEADS
======

[![Build Status](https://travis-ci.org/Przemol/rbeads.svg?branch=master)](https://travis-ci.org/Przemol/rbeads)

rBEADS is the R implementation of Bias Elimination Algorithm for Deep Sequencing.

BEADS is a normalization scheme that corrects nucleotide composition bias, mappability variations and differential local DNA structural effects in deep sequencing data. In high-throughput sequencing data, the recovery of sequenced DNA fragments is not uniform along the genome. In particular, GC-rich sequences are often over-represented and AT-rich sequences under-represented in sequencing data. In addition, the read mapping procedure also generates regional bias. Sequence reads that can be mapped to multiple sites in the genome are usually discarded. Genomic regions with high degeneracy therefore show lower mapped read coverage than unique portions of the genome. Mappability varies along the genome and thus creates systematic bias. Furthermore, local DNA or chromatin structural effects can lead to coverage inhomogeneity of sequencing data.

