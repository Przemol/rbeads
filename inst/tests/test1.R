context("Basic tests")


test_that("Test number 1", {
  expect_that(5 * 2, equals(10))
  expect_that(sqrt(2) ^ 2, equals(2))
  
  expect_equal(
    import.bw('/Volumes/raid0/rBeadsNew/tests2/SummedInputs/Ce10_GAIIxEGSInput_UNIQ_bin25bp_summedInput_linear_1bp.bw', as="RleList"), 
    get(load('/Volumes/raid0/rBeadsNew/tests2/SummedInputs/Ce10_GAIIxEGSInput_UNIQ_bin25bp.Rdata')), 
    tolerance=1e-6
  )
  
  ## Not run: 
  ## expect_that(sqrt(2) ^ 2, is_identical_to(2))
})


context("Basic tests2")
test_that("trigonometric functions match identities", {
  expect_that(sin(pi / 4), equals(1 / sqrt(2)))
  expect_that(cos(pi / 4), equals(1 / sqrt(2)))
  expect_that(tan(pi / 4), equals(1))
})



if (0) {
  setwd('/Volumes/raid0/rBeadsNew/PanGenomicTests')
  beads_bam_bw('H3K79me1^pAB082050_am01^F^let418ts^Embryo_aligned^NA^NA_AM003^F4208555.bam', 
               'Ce10_HiSeqFRMInput_UNIQ_bin25bp.bw', 
               'ce10_mappability_36.bw', 
               'ce10.fa')
  
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