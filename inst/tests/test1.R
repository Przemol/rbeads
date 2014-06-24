
test_that("Test beads with ref genome in fasta file", {
  
  tmp <- file.path(tempdir(), 'BAM_BW_BW_FASTA')
  
  context("BEADS tests for [BAM, BW, BW, FASTA]")
  cat(tmp,'\n')
  expect_that({  
    dir.create(tmp, showWarnings = FALSE); setwd(tmp)
    beads(
      system.file("extdata", "GSM1208360_chrI_100Kb_q5_sample.bam", package="rbeads"), 
      system.file("extdata", "Ce10_HiSeqFRMInput_UNIQ_bin25bp_chrI_100Kb_sample.bw", package="rbeads"), 
      system.file("extdata", "ce10_mappability_chrI_100Kb_sample.bw", package="rbeads"), 
      system.file("extdata", "ce10_chrI_100Kb_sample.fa", package="rbeads")
    ) 
  }, is_true())
  
  expect_true(file.exists( file.path(tmp, 'GSM1208360_chrI_100Kb_q5_sample_BEADS_linear_1bp.bw') ))
  
  cc <- cor(
    as.numeric(IRanges::unlist( import.bw(file.path(tmp, 'GSM1208360_chrI_100Kb_q5_sample_BEADS_linear_1bp.bw'), as="RleList"), use.names=FALSE)), 
    as.numeric(IRanges::unlist( get(load(system.file("extdata", "GSM1208360_expected.Rdata", package="rbeads"))), use.names=FALSE)),
  , use='comp')
  expect_more_than(cc, .99, info=cc)

})

test_that("Test beads with ref genome in BSGenome package", {
  
  tmp <- file.path(tempdir(), 'BAM_BW_BW_BSGenome')
  
  context("BEADS tests for [BAM, BW, BW, BSGenome]")
  cat(tmp,'\n')
  expect_that({  
    dir.create(tmp, showWarnings = FALSE); setwd(tmp)
    beads(
      system.file("extdata", "GSM1208360_chrI_100Kb_q5_sample.bam", package="rbeads"), 
      system.file("extdata", "Ce10_HiSeqFRMInput_UNIQ_bin25bp_chrI_100Kb_sample.bw", package="rbeads"), 
      system.file("extdata", "ce10_mappability_chrI_100Kb_sample.bw", package="rbeads"), 
      genome='ce10'
    ) 
  }, is_true())
  
  expect_true(file.exists( file.path(tmp, 'GSM1208360_chrI_100Kb_q5_sample_BEADS_linear_1bp.bw') ))
  
  cc <- cor(
    as.numeric(IRanges::unlist( import.bw(file.path(tmp, 'GSM1208360_chrI_100Kb_q5_sample_BEADS_linear_1bp.bw'), as="RleList"), use.names=FALSE)), 
    as.numeric(IRanges::unlist( get(load(system.file("extdata", "GSM1208360_expected.Rdata", package="rbeads"))), use.names=FALSE)),
    , use='comp')
  expect_more_than(cc, .99, info=cc)
  
})

test_that("Test beads for two BAM", {
  
  tmp <- file.path(tempdir(), 'BAM_BAM_BW_FATA')
  
  context("BEADS tests for [BAM, BAM, BW, FATA]")
  cat(tmp,'\n')
  expect_that({  
    dir.create(tmp, showWarnings = FALSE); setwd(tmp)
    beads(
      system.file("extdata", "GSM1208360_chrI_100Kb_q5_sample.bam", package="rbeads"),
      system.file("extdata", "GSM1208360_chrI_100Kb_q5_sample.bam", package="rbeads"),
      system.file("extdata", "ce10_mappability_chrI_100Kb_sample.bw", package="rbeads"), 
      system.file("extdata", "ce10_chrI_100Kb_sample.fa", package="rbeads")
    )   
  }, is_true())
  
  expect_true(file.exists( file.path(tmp, 'GSM1208360_chrI_100Kb_q5_sample_BEADS_linear_1bp.bw') ))
  
  cc <- cor(
    as.numeric(IRanges::unlist( import.bw(file.path(tmp, 'GSM1208360_chrI_100Kb_q5_sample_BEADS_linear_1bp.bw'), as="RleList"), use.names=FALSE)), 
    as.numeric(IRanges::unlist( get(load(system.file("extdata", "GSM1208360_expected.Rdata", package="rbeads"))), use.names=FALSE)),
    , use='comp')
  #expect_more_than(cc, .99, info=cc)
  
})