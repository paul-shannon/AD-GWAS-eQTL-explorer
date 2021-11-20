library(RUnit)
library(ADvariantExplorer)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_gwasCatalog()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    targetGene <- "CLU"
    loc.chrom <- "chr8"
    loc.start <- 27447528
    loc.end   <- 27764088
    avx <- ADvariantExplorer$new(targetGene, loc.chrom, loc.start, loc.end)
    checkTrue(all(c("R6", "ADvariantExplorer") %in% class(avx)))
    checkEquals(avx$getTargetGene(), targetGene)

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_gwasCatalog <- function()
{
    message(sprintf("--- test_ctor"))

    targetGene <- "CLU"
    loc.chrom <- "chr8"
    loc.start <- 27447528
    loc.end   <- 27764088
    avx <- ADvariantExplorer$new(targetGene, loc.chrom, loc.start, loc.end)
    tbl.f <- avx$getFilteredGwasTable()
    checkTrue(nrow(tbl.f) < 30)
    tbl.all <- avx$getFullGwasTable(trim.columns=TRUE)
    dim(tbl.all)
    checkTrue(nrow(tbl.all) > 110000)

} # test_gwasCatalog
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
