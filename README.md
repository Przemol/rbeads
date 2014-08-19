rBEADS [![Build Status](https://travis-ci.org/Przemol/rbeads.svg?branch=master)](https://travis-ci.org/Przemol/rbeads)
======

>  **rBEADS** is the R implementation of Bias Elimination Algorithm for Deep Sequencing.

<<<<<<< HEAD
>  **rBEADS** is the R implementation of Bias Elimination Algorithm for Deep Sequencing.

**BEADS** is a normalization scheme that corrects nucleotide composition bias, mappability variations and differential local DNA structural effects in deep sequencing data. In high-throughput sequencing data, the recovery of sequenced DNA fragments is not uniform along the genome. In particular, GC-rich sequences are often over-represented and AT-rich sequences under-represented in sequencing data. In addition, the read mapping procedure also generates regional bias. Sequence reads that can be mapped to multiple sites in the genome are usually discarded. Genomic regions with high degeneracy therefore show lower mapped read coverage than unique portions of the genome. Mappability varies along the genome and thus creates systematic bias. Furthermore, local DNA or chromatin structural effects can lead to coverage inhomogeneity of sequencing data.

## Installation

=======
**BEADS** is a normalization scheme that corrects nucleotide composition bias, mappability variations and differential local DNA structural effects in deep sequencing data. In high-throughput sequencing data, the recovery of sequenced DNA fragments is not uniform along the genome. In particular, GC-rich sequences are often over-represented and AT-rich sequences under-represented in sequencing data. In addition, the read mapping procedure also generates regional bias. Sequence reads that can be mapped to multiple sites in the genome are usually discarded. Genomic regions with high degeneracy therefore show lower mapped read coverage than unique portions of the genome. Mappability varies along the genome and thus creates systematic bias. Furthermore, local DNA or chromatin structural effects can lead to coverage inhomogeneity of sequencing data.

## Installation

>>>>>>> e79f665fb04aa8f17898a1080ee784fb7c81dfd6
First, install required BioConductor packages, by running in R:
```{r}
source("http://bioconductor.org/biocLite.R")
biocLite(c('methods','IRanges','BSgenome','digest','Rsamtools','rtracklayer','GenomicRanges','Biostrings'))
```

To install the latest development version directly from GitHub, run in R:
```{r}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("przemol/rbeads"")
```

## Getting Started

Run following in R to load the library and see package help:
```{r}
library(rbeads)
help(rbeads)
```
