if (0) {
  
  setwd('/Volumes/raid0/rBeadsNew/PanGenomicTests/HUMAN/')
  options(error=recover)
  ex=c('control_readsCoverage', 'control_GCcorected', 'control_GCandMap', 'readsCoverage', 'GCcorected', 'GCandMap', 'BEADS')
  bigtest <- function() {
    beads(
      'GSM733760_hg19_wgEncodeBroadHistoneHsmmH2azStdAlnRep1.bam',
      'GSM733648_hg19_wgEncodeBroadHistoneHsmmtControlStdAlnRep1.bam',
      'wgEncodeCrgMapabilityAlign36mer.bigWig', 
      '/Volumes/raid0/_ref_genomes_/hg19.fa',
      export=ex, rdata=TRUE, quickMap=TRUE
    )   
    gc()
    beads(
      'GSM733760_hg19_wgEncodeBroadHistoneHsmmH2azStdAlnRep2.bam',
      'GSM733648_hg19_wgEncodeBroadHistoneHsmmtControlStdAlnRep2.bam',
      'wgEncodeCrgMapabilityAlign36mer.bigWig', 
      '/Volumes/raid0/_ref_genomes_/hg19.fa',
      export=ex,rdata=TRUE, quickMap=TRUE
    ) 
  }
  
  require(rbeads)
  setwd('/Volumes/raid0/rBeadsNew/PanGenomicTests/HUMAN/HUMAN_K562/')
  options(error=recover)
  ex=c('control_readsCoverage', 'control_GCcorected', 'control_GCandMap', 'readsCoverage', 'GCcorected', 'GCandMap', 'BEADS')
  bigtest <- function() {
    beads(
      'GSM733680_hg19_wgEncodeBroadHistoneK562H3k4me3StdAlnRep1.bam',
      'GSM733780_hg19_wgEncodeBroadHistoneK562ControlStdAlnRep1.bam',
      '/Volumes/raid0/rBeadsNew/PanGenomicTests/HUMAN/wgEncodeCrgMapabilityAlign36mer.bigWig', 
      '/Volumes/raid0/_ref_genomes_/hg19.fa',
      export=ex, rdata=TRUE, quickMap=TRUE
    )   
    gc()
    beads(
      'GSM733680_hg19_wgEncodeBroadHistoneK562H3k4me3StdAlnRep2.bam',
      'GSM733780_hg19_wgEncodeBroadHistoneK562ControlStdAlnRep1.bam',
      '/Volumes/raid0/rBeadsNew/PanGenomicTests/HUMAN/wgEncodeCrgMapabilityAlign36mer.bigWig', 
      '/Volumes/raid0/_ref_genomes_/hg19.fa',
      export=ex,rdata=TRUE, quickMap=TRUE
    ) 
  }
  
  #weighted.mean(sum(is.na(sample.norm))/elementLengths(sample.norm), elementLengths(sample.norm))
  toBW(sample.norm, 'BEADS_withNAN', 'rep1_test')
  toBW(sample.norm, 'BEADS_NAto0', 'rep1_test', missing2zero=TRUE)
  toBW_missing(sample.norm, 'BEADS_missing', 'rep1_test')
  
  beads_bam_bam(
    system.file("extdata", "GSM1208360_chrI_100Kb_q5_sample.bam", package="rbeads"),
    system.file("extdata", "GSM1208360_chrI_100Kb_q5_sample.bam", package="rbeads"),
    system.file("extdata", "ce10_mappability_chrI_100Kb_sample.bw", package="rbeads"), 
    system.file("extdata", "ce10_chrI_100Kb_sample.fa", package="rbeads")
  )  
  
  setwd('/Volumes/raid0/rBeadsNew/PanGenomicTests/Worm_Cfp1/SAMPLE_test/')
  beads_bam_bw(
    'GSM1208360_chrI_100Kb_q5_sample.bam',
    '/Volumes/raid0/rBeadsNew/PanGenomicTests/ce10_commons/Ce10_HiSeqFRMInput_UNIQ_bin25bp.bw',
    '/Volumes/raid0/rBeadsNew/PanGenomicTests/ce10_commons/ce10_mappability_36.bw',
    '/Volumes/raid0/rBeadsNew/PanGenomicTests/ce10_commons/ce10.fa'
  )
  
  setwd('/Volumes/raid0/rBeadsNew/PanGenomicTests/Worm_H3K4me3/')
  beads(
    'H3K4me3^ab8580_rc03^F^Cfp1-GFP^L3_aligned^NA^NA_RC007^F3b06602.bam',
    '/Volumes/raid0/rBeadsNew/PanGenomicTests/ce10_commons/Ce10_HiSeqFRMInput_UNIQ_bin25bp.bw',
    '/Volumes/raid0/rBeadsNew/PanGenomicTests/ce10_commons/ce10_mappability_36.bw',
    '/Volumes/raid0/rBeadsNew/PanGenomicTests/ce10_commons/ce10.fa',
    rdata=TRUE, quickMap=FALSE
  )
  
  setwd('/Volumes/raid0/rBeadsNew/PanGenomicTests/Worm_Cfp1')
  beads(
    'antiGFP^ab_rc04^F^Cfp1-GFP^L3_aligned^NA^NA_RC010^Ffd06602.bam', 
    'Ce10_HiSeqFRMInput_UNIQ_bin25bp.bw', 
    'ce10_mappability_36.bw', 
    'ce10.fa',rdata=TRUE, quickMap=FALSE
  )
  
  oldwd <- getwd()
  setwd('/Volumes/raid0/rBeadsNew/tests2/SummedInputs/')
  toBW(load('Ce10_GAIIxEGSInput_UNIQ_bin25bp.Rdata'), 'summedInput','Ce10_GAIIxEGSInput_UNIQ_bin25bp')
  toBW(load('Ce10_GAIIxFRMInput_UNIQ_bin25bp.Rdata'), 'summedInput','Ce10_GAIIxFRMInput_UNIQ_bin25bp')
  toBW(load('Ce10_HiSeqEGSInput_UNIQ_bin25bp.Rdata'), 'summedInput','Ce10_HiSeqEGSInput_UNIQ_bin25bp')
  toBW(load('Ce10_HiSeqFRMInput_UNIQ_bin25bp.Rdata'), 'summedInput','Ce10_HiSeqFRMInput_UNIQ_bin25bp')
  test(oldwd)
  
  
  
  testGCscores <- function( files=dir(pattern="RawRanges.+Rdata") ) {
    source('/Users/przemol/Documents/workspace/RBeads/rBeads/rB_enriched_regions_one_strand_v1.0.R')
    out = list()
    for(file in files) {
      cat("-> GC TESTING:", file, "\n")
      varname <- load(file)
      d <- sprintf("%s", unlist(strsplit(file, "\\."))[1])
      ER <- rB.EnrichedRegions.OS(get(varname), desc=d )
      GC <- rB.sumGCscores(get(varname), enriched_regions=ER, desc=d, smoothing_spline=FALSE)
      out <- c(out, GC)
    }
    names(out) <- files
    return(out)
  }
  
  
}