rBEADS
======

> The R implementation of <strong>B</strong>ias <strong>E</strong>limination <strong>A</strong>lgorithm for <strong>D</strong>eep <strong>S</strong>equencing.

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.11427.png)](http://dx.doi.org/10.5281/zenodo.11427)
[![Build Status](https://travis-ci.org/Przemol/rbeads.svg?branch=master)](https://travis-ci.org/Przemol/rbeads)

:exclamation:  **RELEASE NOTE**  :exclamation:

rBEADS is in pre-release (alpha) stage. The software is provided for testing purposes. Please report the problem, bugs, unexpected behaviors and missing features [here](../../issues/new).

BEADS algorithm requires deep inputs (high reads coverage) to work properly. This means >50 million reads for worm and fly experiments and proportionally higher number for mammalian experiments. It is suggested to pool multiple input experiments using ```sumBAMinputs``` function from rBEADS package.

## Introduction

**BEADS** is a normalization scheme that corrects nucleotide composition bias, mappability variations and differential local DNA structural effects in deep sequencing data. In high-throughput sequencing data, the recovery of sequenced DNA fragments is not uniform along the genome. In particular, GC-rich sequences are often over-represented and AT-rich sequences under-represented in sequencing data. In addition, the read mapping procedure also generates regional bias. Sequence reads that can be mapped to multiple sites in the genome are usually discarded. Genomic regions with high degeneracy therefore show lower mapped read coverage than unique portions of the genome. Mappability varies along the genome and thus creates systematic bias. Furthermore, local DNA or chromatin structural effects can lead to coverage inhomogeneity of sequencing data.

## Installation

First, install required BioConductor packages, by running in R:
```{r}
source("http://bioconductor.org/biocLite.R")
biocLite(c('methods','IRanges','BSgenome','digest','Rsamtools','rtracklayer','GenomicRanges','Biostrings'))
```

To install the latest development version directly from GitHub, run in R:
```{r}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("przemol/rbeads")
```

## Getting Started

Run following in R to load the library and see package help:
```{r}
library(rbeads)
help(rbeads)
```
R style reference manual (PDF) can be found [**here**](https://github.com/Przemol/rbeads/releases/download/v0.3.1-alpha/rbeads.pdf).

## Mappability tracks

Following pre-calculated mappabiliti tracks (BigWig files) are avilable ta the moment:
* [```ce10_gem-mappability_36bp.bw```](https://github.com/Przemol/rbeads/releases/download/v0.3.1-alpha/ce10_gem-mappability_36bp.bw) - *C. elegans* mappability track for 36bp reads
* [```dm3_gem-mappability_36bp.bw```](https://github.com/Przemol/rbeads/releases/download/v0.3.1-alpha/dm3_gem-mappability_36bp.bw) - *D. melanogaster* mappability track for 36bp reads

Human tracks from [UCSC](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/):
* [```wgEncodeCrgMapabilityAlign24mer.bigWig```](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign24mer.bigWig) - *H. sapiens* mappability track for 24bp reads
* [```wgEncodeCrgMapabilityAlign36mer.bigWig```](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig) - *H. sapiens* mappability track for 36bp reads
* [```wgEncodeCrgMapabilityAlign40mer.bigWig```](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign40mer.bigWig) - *H. sapiens* mappability track for 40bp reads
* [```wgEncodeCrgMapabilityAlign50mer.bigWig ```](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig) - *H. sapiens* mappability track for 50bp reads
* [```wgEncodeCrgMapabilityAlign75mer.bigWig```](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig) - *H. sapiens* mappability track for 75bp reads
* [```wgEncodeCrgMapabilityAlign100mer.bigWig ```](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig) - *H. sapiens* mappability track for 100bp reads
