library(RUnit)
library(ADvariantExplorer)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   studyIDs <- "GTEx_ge_brain_frontal_cortex"
   chrom <- "8"
   start <- 27610984
   end <- 27610987

   fetcher <- EQTLCatalogFetcher$new(studyIDs, chrom, start, end, simplify=TRUE)
   checkTrue(all(c("EQTLCatalogFetcher", "R6") %in% class(fetcher)))

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_fetchFrom_1_study <- function()
{
   message(sprintf("--- test_fetchFrom_1_study"))

   studyIDs <- "GTEx_ge_brain_frontal_cortex"
   chrom <- "8"
   start <- 27610984
   end <- 27610987

   fetcher <- EQTLCatalogFetcher$new(studyIDs, chrom, start, end, simplify=TRUE)
   tbl <- fetcher$fetch()

} # test_fetchFrom_1_study
#----------------------------------------------------------------------------------------------------
