library(RUnit)
library(ADvariantExplorer)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_gwasCatalog()
    test_eqtlCatalogVariants()
    test_eqtCatalogVariants_combineSlightlyDiscordantStudies()

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
test_tablixToEBIeQTLCatalog <- function()
{
  tbl <- fetch_restAPI(unique_id="Alasoo_2018.macrophage_naive",
                                chrom="8", bp_lower=27603335, bp_upper=27608281)

   cmd <- sprintf("%s %s/%s %s",
                   "tabix",
                   "ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv",
                   "Alasoo_2018/ge/Alasoo_2018_ge_macrophage_naive.all.tsv.gz",
                   "8:27603335-27608281")
   system(cmd)

} # test_tabixToEBIdQTLCatalog
#----------------------------------------------------------------------------------------------------
test_gwasCatalog <- function()
{
    message(sprintf("--- test_gwasCatalog"))

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
test_eqtlCatalogSummary <- function()
{
    message(sprintf("--- test_eqtlCatalog"))

    targetGene <- "CLU"
    loc.chrom <- "chr8"
    loc.start <- 27447528
    loc.end   <- 27764088

    avx <- ADvariantExplorer$new(targetGene, loc.chrom, loc.start, loc.end)
    tbl.cat <- avx$geteQTLSummary()
    study.ids <- unique(tbl.cat$unique_id)
    checkEquals(grep("brain_naive", study.ids, value=TRUE), "ROSMAP.brain_naive")

} # test_eqtlCatalogSummary
#----------------------------------------------------------------------------------------------------
test_eqtlCatalogStudySearch <- function()
{
    message(sprintf("--- test_eqtlCatalogStudySearch"))

    targetGene <- "CLU"
    loc.chrom <- "chr8"
    loc.start <- 27447528
    loc.end   <- 27764088

    avx <- ADvariantExplorer$new(targetGene, loc.chrom, loc.start, loc.end)

    studies <- avx$geteqtlStudyNamesForGroup("macrophage_naive")
    checkEquals(studies, c("Alasoo_2018.macrophage_naive", "Nedelec_2016.macrophage_naive"))

    studies.expected.to.fail <- avx$geteqtlStudyNamesForGroup("bogus")
    checkTrue(is.na(studies.expected.to.fail))

    studies <- avx$geteqtlStudyNamesForGroup("brain")
    checkEquals(length(studies), 15)

} # test_eqtlCatalogSummarySearch
#----------------------------------------------------------------------------------------------------
test_eqtlCatalogVariants <- function()
{
    message(sprintf("--- test_eqtlCatalog"))

    targetGene <- "CLU"
    loc.chrom <- "chr8"
    loc.start <- 27447528
    loc.end   <- 27764088
    avx <- ADvariantExplorer$new(targetGene, loc.chrom, loc.start, loc.end)
    tbl.cat <- avx$geteQTLSummary()
    dim(tbl.cat) # 397 12
    checkTrue(nrow(tbl.cat) > 390)
    checkEquals(ncol(tbl.cat), 12)

       # find study ids like this
    sort(unique(tbl.cat$unique_id))   # 112 of them (24 nov 2021)
    sort(unique(grep("macrophage_naive", tbl.cat$unique_id, value=TRUE)))

    study.1 <- "Alasoo_2018.macrophage_naive"
    study.2 <- "Nedelec_2016.macrophage_naive"
    checkTrue(all(c(study.1, study.2) %in% grep("macrophage_naive", tbl.cat$unique_id, value=TRUE)))
    chrom <- "8"
    start <- 27603335
    end   <- 27608281
    1 + end - start
    tbl.1 <- avx$geteQTLsByLocationAndStudyID(chrom, start, end, study.1, simplify=TRUE)

    tbl.2 <- avx$geteQTLsByLocationAndStudyID(chrom, start, end, study.2, simplify=TRUE)
    tbl.12 <- avx$geteQTLsByLocationAndStudyID(chrom, start, end, c(study.1, study.2), simplify=TRUE)
    checkEquals(nrow(tbl.1) + nrow(tbl.2), nrow(tbl.12))

    tbl.rosmap.brain <- avx$geteQTLsByLocationAndStudyID(chrom, start, end,
                                                          "ROSMAP.brain_naive", simplify=TRUE)
    dim(tbl.rosmap.brain)
    checkTrue(nrow(tbl.rosmap.brain) > 160)

    # this fails with a 404
    #tbl.gtex.bc <- avx$geteQTLsByLocationAndStudyID(chrom, start, end, "GTEx.brain_cortex", simplify=TRUE)
    #dim(tbl.gtex.bc)
    #checkTrue(nrow(tbl.gtex.bc) > 430)

} # test_eqtlCatalogVariants
#----------------------------------------------------------------------------------------------------
# "BrainSeq.brain"     "ROSMAP.brain_naive" have 23 and 25 columns, requring rbind with fill
test_eqtCatalogVariants_combineSlightlyDiscordantStudies <- function()
{
    message(sprintf("--- test_eqtCatalogVariants_combineSlightlyDiscordantStudies"))

    targetGene <- "CLU"
    chrom <- rep("8", 3)
    start <- c(27354393, 27337603,27362469)
    end   <- c(27354394, 27337604, 27362470)
    avx <- ADvariantExplorer$new(targetGene, chrom[2], start[2], end[2])

    studies <- c("BrainSeq.brain",  "ROSMAP.brain_naive")
    i <- 2
    tbl <- avx$geteQTLsByLocationAndStudyID(chrom[i], start[i], end[i], studies, method="REST", simplify=TRUE)
    checkTrue(nrow(tbl) > 35)
    checkEquals(ncol(tbl), 6)

} # test_eqtCatalogVariants_combineSlightlyDiscordantStudies
#----------------------------------------------------------------------------------------------------
# rs867230  8:27610986 (GRCh38)
viz <- function()
{
    igv <- start.igv("CLU")
    tbl.track <- data.frame(chrom="chr8", start=27610985, end=27610986, stringsAsFactors=FALSE)
    track <- DataFrameAnnotationTrack("rs867230", tbl.track, color="brown")
    displayTrack(igv, track)

} # viz
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
