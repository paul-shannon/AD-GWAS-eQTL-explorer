library(RUnit)
library(ADvariantExplorer)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_simpleFetch()
    test_fullGwasCatalog()
    test_filteredGwasCatalog_CLU()

    test_eqtlCatalogVariants()
    test_eqtCatalogVariants_combineSlightlyDiscordantStudies()
    test_getCellTypes()

} # runTests
#----------------------------------------------------------------------------------------------------
# defaults to method=REST, depends upon dev branch as of (27 mar 2022)
test_simpleFetch <- function()
{

  require(catalogueR)

        # 10 bases on either side of rs4575098
   loc.chrom <- "1"
   loc.start <- 161185592
   loc.end   <- 161185612
   study <- "GTEx_V8.Brain_Hippocampus"

   tbl <- eQTL_Catalogue.fetch(unique_id=study,
                               quant_method="ge",
                               #use_tabix=FALSE,
                               chrom = loc.chrom,
                               bp_lower=loc.start,
                               bp_upper=loc.end,
                               verbose=TRUE)
   checkEquals(ncol(tbl), 24)
   checkTrue(nrow(tbl) > 50)   # 57 on

   subset(tbl, pvalue.QTL < 1e-2)[, c("gene.QTL", "pvalue.QTL", "beta.QTL", "rsid.QTL")]

   targetGene <- "NDUFS2"
   require(ADvariantExplorer)

   avx <- ADvariantExplorer$new(targetGene, loc.chrom, loc.start, loc.end)
   tbl.2 <- avx$geteQTLsByLocationAndStudyID(loc.chrom, loc.start, loc.end, study, simplify=TRUE)
   tbl.top <- subset(tbl.2, pvalue <= 0.01)
   checkTrue(nrow(tbl.top) >= 1)
   #        rsid      pvalue   gene total.alleles      beta                        id
   # 1 rs4575098 0.000259162 NDUFS2           330 -0.222035 GTEx_V8.Brain_Hippocampus
   # 2 rs4575098 0.007248770  TSTD1           330  0.205420 GTEx_V8.Brain_Hippocampus

} # test_simpleFetch
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
test_fullGwasCatalog <- function()
{
    message(sprintf("--- test_fullGwasCatalog"))

    targetGene <- "CLU"
    loc.chrom <- "chr8"
    loc.start <- 27447528
    loc.end   <- 27764088
    avx <- ADvariantExplorer$new(targetGene, loc.chrom, loc.start, loc.end)

       #------------------------------------------------------------
       # get the entire table, with and without trimming columns
       #------------------------------------------------------------

    tbl.all <- avx$getFullGwasTable(trim.columns=TRUE)
    dim(tbl.all)
    checkTrue(nrow(tbl.all) > 110000)
    expected <- c("SNPS", "P.VALUE", "MAPPED_GENE", "FIRST.AUTHOR", "DATE", "INITIAL.SAMPLE.SIZE", "STUDY")
    checkEquals(colnames(tbl.all), expected)

    tbl.all <- avx$getFullGwasTable(trim.columns=FALSE)
    checkTrue(nrow(tbl.all) > 110000)
    checkTrue(ncol(tbl.all) > 40)

    checkTrue(length(unique(tbl.all$STUDY)) > 3150)

       #----------------------------------------------------
       #  identify the studies with alzheimer in their name
       #----------------------------------------------------

    ad.studies <- unique(grep("alzh", tbl.all$STUDY, ignore.case=TRUE, v=TRUE))
    checkTrue(length(ad.studies) > 65)

} # test_fullGwasCatalog
#----------------------------------------------------------------------------------------------------
# CLU is a long-standing AD locus, primarily associated with rs11136000
#    data.dir <- "~/github/TrenaProjectAD/inst/extdata/gwasLoci"
#    tbl.wms <- get(load(file.path(data.dir, "williams-natureNeuroscience2020.RData")))
#    tbl.posthuma <- get(load(file.path(data.dir, "tbl.posthuma-38-geneAssociations-curated-3828x12.RData")))

test_filteredGwasCatalog_CLU <- function()
{
    message(sprintf("--- test_filteredGwasCatalog"))

    targetGene <- "CLU"
    loc.chrom <- "chr8"
    loc.start <- 27447528
    loc.end   <- 27764088
    avx <- ADvariantExplorer$new(targetGene, loc.chrom, loc.start, loc.end)

    tbl.1 <- avx$getFilteredGwasTable(targetGeneOnly=FALSE)
    checkTrue(length(tbl.1$MAPPED_GENE) > 1)
    checkTrue(nrow(tbl.1) > 20)

    tbl.2 <- avx$getFilteredGwasTable(targetGeneOnly=TRUE)
    mapped.genes <- unique(tbl.2$MAPPED_GENE)
      # e.g.,  "CLU", "CLU - SCARA3"
    checkEquals(length(grep(targetGene, mapped.genes, ignore.case=TRUE)),
                length(mapped.genes))
    checkTrue(nrow(tbl.1) > nrow(tbl.2))

    tbl.3 <- avx$getFilteredGwasTable(targetGeneOnly=FALSE, studyString="")
    checkTrue(nrow(tbl.3) > 50)
    checkTrue(length(unique(tbl.3$STUDY)) > length(unique(tbl.1$STUDY)))
    mapped.genes <- unique(tbl.3$MAPPED_GENE)
    checkTrue(length(grep(targetGene, mapped.genes, ignore.case=TRUE)) < length(mapped.genes))

       #------------------------------------------------------------
       # now get the entire table, with and without trimming columns
       #------------------------------------------------------------

    tbl.all <- avx$getFullGwasTable(trim.columns=TRUE)
    dim(tbl.all)
    checkTrue(nrow(tbl.all) > 110000)

    tbl.all.noTrim <- avx$getFullGwasTable(trim.columns=FALSE)
    checkTrue(nrow(tbl.all.noTrim) > 110000)
    checkTrue(ncol(tbl.all.noTrim) > 40)

} # test_filteredGwasCatalog_CLU
#----------------------------------------------------------------------------------------------------
obtainAndSaveAdGwasVariants <- function()
{
    message(sprintf("--- obtainAndSaveAdGwasVariants"))

       # these mandatory parameters are mere expediences, and are not actually used in what follows
    targetGene <- "CLU"
    loc.chrom <- "chr8"
    loc.start <- 27447528
    loc.end   <- 27764088
    avx <- ADvariantExplorer$new(targetGene, loc.chrom, loc.start, loc.end)

       #------------------------------------------------------------
       # get the entire table, with and without trimming columns
       #------------------------------------------------------------

    tbl.all <- avx$getFullGwasTable(trim.columns=FALSE)
    checkTrue(nrow(tbl.all) > 100000)
    checkEquals(colnames(tbl.all), c("SNPS","P.VALUE","MAPPED_GENE","FIRST.AUTHOR","DATE","INITIAL.SAMPLE.SIZE","STUDY"))
    ad.studies <- sort(unique(grep("alz", tbl.all$STUDY, ignore.case=TRUE, value=TRUE)))
    checkTrue(length(ad.studies) >= 7)
    as.data.frame(sort(table(grep("alz", tbl.all$STUDY, ignore.case=TRUE, value=TRUE))))
    tbl.marioni <- subset(tbl.all, STUDY=="GWAS on family history of Alzh")
    dim(tbl.marioni)
    tbl.posthuma <- subset(tbl.all, STUDY=="Genome-wide meta-analysis identifies new loci and functional pathways influencing Alzheimer's disease risk.")
    dim(tbl.posthuma) # 174 43
    tbl.schwartzentruber <- subset(tbl.all,
                                   STUDY=="Genome-wide meta-analysis, fine-mapping and integrative prioritization implicate new Alzheimer's disease risk genes.")
    dim(tbl.schwartzentruber) # 30 43

    tbl.small <- unique(tbl.posthuma[, c("SNPS", "SNP_GENE_IDS", "P.VALUE")])
    dim(tbl.small)          # 123 2
    library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
    gr.snps <- snpsById(SNPlocs.Hsapiens.dbSNP151.GRCh38, tbl.small$SNPS)
    tbl.snps <- as.data.frame(gr.snps)
    dim(tbl.snps)
    colnames(tbl.snps) <- c("chrom", "hg38", "strand", "rsid", "allele")
    tbl.snps$chrom <- as.character(tbl.snps$chrom)
    tbl.snps <- tbl.snps[, c("chrom", "hg38", "rsid", "allele")]
    tbl.snps <- merge(tbl.snps, tbl.small[, c("SNPS", "P.VALUE")], by.x="rsid", by.y="SNPS", all.x=TRUE)
    dups <- which(duplicated(tbl.snps[, 1:3]))
    tbl.snps <- tbl.snps[-dups,]
    colnames(tbl.snps)[5] <- "pvalue"
    dim(tbl.snps)  # 123 5
    save(tbl.snps, file="~/github/TrenaProjectAD/inst/extdata/gwasLoci/posthuma-2019-with-hg38-locs.RData")

    tbl.small <- unique(tbl.schwartzentruber[, c("SNPS", "SNP_GENE_IDS", "P.VALUE")])
    dim(tbl.small)
    gr.snps <- snpsById(SNPlocs.Hsapiens.dbSNP151.GRCh38, tbl.small$SNPS, ifnotfound="drop")
    tbl.snps <- as.data.frame(gr.snps)
    colnames(tbl.snps) <- c("chrom", "hg38", "strand", "rsid", "allele")
    tbl.snps$chrom <- as.character(tbl.snps$chrom)
    tbl.snps <- tbl.snps[, c("chrom", "hg38", "rsid", "allele")]
    dim(tbl.snps)
    tbl.snps <- merge(tbl.snps, tbl.small[, c("SNPS", "P.VALUE")], by.x="rsid", by.y="SNPS", all.x=TRUE)
    dups <- which(duplicated(tbl.snps[, 1:3]))
    if(length(dups) > 0)
       tbl.snps <- tbl.snps[-dups,]
    colnames(tbl.snps)[5] <- "pvalue"

    lapply(tbl.snps, class)
    save(tbl.snps, file="~/github/TrenaProjectAD/inst/extdata/gwasLoci/schwartzentruber-2021-with-hg38-locs.RData")


} # obtainAndSaveAdGwasVariants
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
    checkTrue(nrow(tbl.cat) > 490)
    checkEquals(ncol(tbl.cat), 12)
      # the unique_id values are used in subsequent fetch of actual datasets
    study.ids <- unique(tbl.cat$unique_id)
    checkEquals(grep("brain_naive", study.ids, value=TRUE), "ROSMAP.brain_naive")
       # GTEx_V8 eqtls added to catalogueR, dev branch, on (1 dec 2021)
    checkTrue(length(grep("GTEx_V8", study.ids)) >= 49)

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
    dim(tbl.cat) # 498 12
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
    checkTrue(nrow(tbl.rosmap.brain) > 8)

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
    tbl <- avx$geteQTLsByLocationAndStudyID(chrom[i], start[i], end[i], studies, simplify=TRUE)
    checkTrue(nrow(tbl) > 1)
    checkEquals(ncol(tbl), 6)

} # test_eqtCatalogVariants_combineSlightlyDiscordantStudies
#----------------------------------------------------------------------------------------------------
test_eqtCatalogVariants_AD_associated_NDUFS2 <- function()
{
   message(sprintf("--- test_eqtCatalogVariants_AD_associated_NDUFS2"))

   targetGene <- "NDUFS2"
   tag.snp <- "rs4575098"
   tag.snp.loc <- 161185602
   chrom <- "1"
   ld.snp  <- "rs11585858"

   shoulder <- 100
   avx <- ADvariantExplorer$new(targetGene, chrom, tag.snp.loc-shoulder, tag.snp.loc+shoulder)
   tbl.gwas <- avx$getFilteredGwasTable()
   dim(tbl.gwas)
   tbl.gwas

   studies <- c("BrainSeq.brain",  "ROSMAP.brain_naive")
   i <- 2
   #tbl <- avx$geteQTLsByLocationAndStudyID(chrom, tag.snp.loc-shoulder, tag.snp.loc+shoulder,
   #                                        studies, simplify=TRUE)

   tbl.eCat <- avx$geteQTLSummary()
   dim(tbl.eCat)
   checkTrue(nrow(tbl.eCat) > 450)
   checkTrue(ncol(tbl.eCat) >= 10)

   brain.geneExpression.studies <- subset(tbl.eCat, grepl("brain", tissue_label, ignore.case=TRUE) &
                                                   quant_method=="ge")$unique_id
   gtex.v8.brain.studies <- grep("GTEx_V8", brain.geneExpression.studies, value=TRUE)
   length(gtex.v8.brain.studies)
   checkTrue(length(gtex.v8.brain.studies) > 10 & length(gtex.v8.brain.studies) < 20)

   tbl.eqtl <- avx$geteQTLsByLocationAndStudyID(chrom, tag.snp.loc-shoulder, tag.snp.loc+shoulder,
                                           gtex.v8.brain.studies, simplify=TRUE)
   tbl.eqtl.ndufs2 <- subset(tbl.eqtl, gene=="NDUFS2")
   checkTrue(nrow(tbl.eqtl.ndufs2) > 10)
   checkTrue(nrow(tbl.eqtl.ndufs2) < 20)

   also.examine.the.variant <- function(){
      library(EndophenotypeExplorer)
      etx <- EndophenotypeExplorer$new("NDUFS2", "hg38", vcf.project="ADNI")
      etx$getAggregatedAlleleFrequencies("rs4575098", quiet=FALSE)
      mtx.geno <- etx$getGenoMatrixByRSID("rs4575098")
      }

} # test_eqtCatalogVariants_AD_associated_NDUFS2
#----------------------------------------------------------------------------------------------------
test_getCellTypes <- function()
{
   message(sprintf("--- test_getCellTypes"))

   targetGene <- "NDUFS2"
   tag.snp <- "rs4575098"
   tag.snp.loc <- 161185602
   chrom <- "1"
   ld.snp  <- "rs11585858"

   shoulder <- 10000
   avx <- ADvariantExplorer$new(targetGene, chrom, tag.snp.loc-shoulder, tag.snp.loc+shoulder)
   tbl.meis2 <- avx$getCellTypes("MEIS2")
   checkTrue(nrow(tbl.meis2) > 500)
   meis2.dist <- as.list(sort(table(tbl.meis2$ct), decreasing=TRUE))
   checkEquals(names(meis2.dist), c("exc", "ast", "oli", "opc", "mic", "inh", "end"))
      # excitatory neurons, astrocytes, oligodendrocytes, oligodendroctye precursor cells,
      # microglia, inhibitory neurons, endothelial, [pericytes]

   tbl.all <- avx$getCellTypes()
   checkTrue(nrow(tbl.all) >  3250000)
      #  "exc" "mic" "ast" "inh" "oli" "end" "opc" "per" "unk" "nk "

} # test_getCellTypes
#----------------------------------------------------------------------------------------------------
explore.failed.fetches <- function()
{
   targetGene <- "NDUFS2"
   tag.snp <- "rs4575098"
   tag.snp.loc <- 161185602
   chrom <- "1"
   ld.snp  <- "rs11585858"

   shoulder <- 10000
   avx <- ADvariantExplorer$new(targetGene, chrom, tag.snp.loc-shoulder, tag.snp.loc+shoulder)

   tbl.cat <- avx$geteQTLSummary()
   studies <- unique(subset(tbl.cat, quant_method=="ge")$unique_id)
   length(studies)  # 157
   old.gtex.notWorking <- grep("GTEx.", studies, fixed=TRUE)
   length(old.gtex.notWorking) # 49
   studies <- studies[-old.gtex.notWorking]
   length(studies)   # 108
   tbl <- avx$geteQTLsByLocationAndStudyID(chrom, tag.snp.loc-shoulder, tag.snp.loc+shoulder,
                                           studies, simplify=TRUE)

       #----------------------------------------
       # 25/108 failures (1 dec 2021)
       #----------------------------------------
    failed.studies <- c("Alasoo_2018.macrophage_IFNg+Salmonella",
                        "Schmiedel_2018.Tfh_memory",
                        "Schmiedel_2018.Th17_memory",
                        "Schmiedel_2018.Th1_memory",
                        "Schmiedel_2018.Th2_memory",
                        "Schmiedel_2018.Th1-17_memory",
                        "Schmiedel_2018.B-cell_naive",
                        "Schmiedel_2018.CD4_T-cell_naive",
                        "Schmiedel_2018.CD8_T-cell_naive",
                        "Schmiedel_2018.monocyte_CD16_naive",
                        "Schmiedel_2018.monocyte_naive",
                        "Schmiedel_2018.NK-cell_naive",
                        "Braineac2.putamen",
                        "Braineac2.substantia_nigra",
                        "CommonMind.DLPFC_naive",
                        "CAP.LCL_naive",
                        "CAP.LCL_statin",
                        "iPSCORE.iPSC",
                        "Peng_2018.placenta_naive",
                        "PhLiPS.iPSC",
                        "PhLiPS.HLC",
                        "Steinberg_2020.synovium_naive",
                        "Steinberg_2020.low_grade_cartilage_naive",
                        "Steinberg_2020.high_grade_cartilage_naive",
                        "Young_2019.microglia_naive")

   tbl <- avx$geteQTLsByLocationAndStudyID(chrom, tag.snp.loc-shoulder, tag.snp.loc+shoulder,
                                           failed.studies, simplify=TRUE)


   dim(tbl)
   save(tbl, file="~/github/ADvariantExplorer/explore/adamts4/tbl.eqtls.25initialFailures.61389x6.RData")

    # error message, in console and in web browser
    #  {"status": 404, "error": "Resource Not Found", "message":
    #  "Not found. Study :Alasoo_2018 with qtl_group: macrophage_IFNg Salmonella and quantification method: ge"}
   subset(tbl.cat, unique_id == study)

} # explore.failed.fetches
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
