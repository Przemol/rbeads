#Load system files
sample_bam <- system.file("extdata", "GSM1208360_chrI_100Kb_q5_sample.bam", package="rbeads")
input_bam <- system.file("extdata", "Input_fE3_AA169.bam", package="rbeads")
SummedInput_bw <- system.file("extdata", "Ce10_HiSeqFRMInput_UNIQ_bin25bp_chrI_100Kb_sample.bw", package="rbeads")
map_bw <- system.file("extdata", "ce10_mappability_chrI_100Kb_sample.bw", package="rbeads")
ref_fa <- system.file("extdata", "ce10_chrI_100Kb_sample.fa", package="rbeads")

test_that("Test beads with ref genome in fasta file", {
  
  tmp <- file.path(tempdir(), 'BAM_BW_BW_FASTA')
  
  context("BEADS tests for [BAM, BW, BW, FASTA]")
  dir.create(tmp, showWarnings = FALSE); setwd(tmp)
  
  out <- beads(sample_bam, SummedInput_bw, map_bw, ref_fa) 
  expect_is(out, "BigWigFileList")
  expect_true(file.exists( file.path(tmp, 'GSM1208360_chrI_100Kb_q5_sample_BEADS_linear_1bp.bw') ))
  
  cc <- cor(
    as.numeric(unlist( import.bw(file.path(tmp, 'GSM1208360_chrI_100Kb_q5_sample_BEADS_linear_1bp.bw'), as="RleList"), use.names=FALSE)), 
    as.numeric(unlist( get(load(system.file("extdata", "GSM1208360_expected.Rdata", package="rbeads"))), use.names=FALSE)),
  use='comp')
  expect_gt(cc, .99)

})

test_that("Test beads with ref genome in BSGenome package", {
  
  tmp <- file.path(tempdir(), 'BAM_BW_BW_BSGenome')
  
  context("BEADS tests for [BAM, BW, BW, BSGenome]")
  dir.create(tmp, showWarnings = FALSE); setwd(tmp)
  out <- beads(sample_bam, SummedInput_bw, map_bw, genome='ce10') 
  expect_is(out, "BigWigFileList")
  expect_true(file.exists( file.path(tmp, 'GSM1208360_chrI_100Kb_q5_sample_BEADS_linear_1bp.bw') ))
  
  cc <- cor(
    as.numeric(unlist( import.bw(file.path(tmp, 'GSM1208360_chrI_100Kb_q5_sample_BEADS_linear_1bp.bw'), as="RleList"), use.names=FALSE)), 
    as.numeric(unlist( get(load(system.file("extdata", "GSM1208360_expected.Rdata", package="rbeads"))), use.names=FALSE)),
    use='comp')
  expect_gt(cc, .99)
  
})

test_that("Test beads for two BAM", {
  
  tmp <- file.path(tempdir(), 'BAM_BAM_BW_FATA')
  
  context("BEADS tests for [BAM, BAM, BW, FATA]")
  dir.create(tmp, showWarnings = FALSE); setwd(tmp)
  out <- beads(sample_bam, input_bam, map_bw, ref_fa) 
  expect_is(out, "BigWigFileList")
  expect_true(file.exists( file.path(tmp, 'GSM1208360_chrI_100Kb_q5_sample_BEADS_linear_1bp.bw') ))
  
  cc <- cor(
    as.numeric(unlist( import.bw(file.path(tmp, 'GSM1208360_chrI_100Kb_q5_sample_BEADS_linear_1bp.bw'), as="RleList"), use.names=FALSE)), 
    as.numeric(unlist( get(load(system.file("extdata", "GSM1208360_expected.Rdata", package="rbeads"))), use.names=FALSE)),
    use='comp')
  expect_gt(cc, .94)
  
})


test_that("Test summed input", {
  
  tmp <- file.path(tempdir(), 'SummedInputTest')
  
  context("Binning helper function")
  track <- import.bw(system.file("extdata", "Ce10_HiSeqFRMInput_UNIQ_bin25bp_chrI_100Kb_sample.bw", package="rbeads"), as="RleList")
  expect_equal(as.integer(runValue(binRle2Rle(track, 100000L))), as.integer(mean(as.numeric(unlist(track, use.names=FALSE)), na.rm=TRUE)))
  
  context("Summed input creation")
  a <- rtracklayer::summary(sumBAMinputs(c(input_bam, input_bam), map_bw, ref_fa, out_name = tmp))
  b <- rtracklayer::summary(sumBAMinputs(c(input_bam, input_bam, input_bam), map_bw, ref_fa, out_name = tmp))
  expect_equal(a[[1]]$score*(3/2), b[[1]]$score, tolerance=0.01)
  
})