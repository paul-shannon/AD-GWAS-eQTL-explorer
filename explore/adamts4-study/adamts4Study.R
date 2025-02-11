# adamts4Study: explore gwas, eQTL, rna-seq and trena around this traditional AD locus
# a few variants in the vicinity of ADAMTS4 are highly associated with AD
# conventional wisdom used to assert that this/these variants, encompassed by ADAMTS4, affected
# that gene.  nowadays the community prefers other genes, or is agnostic.
# we try to shed light on that question here.
#----------------------------------------------------------------------------------------------------
library(EndophenotypeExplorer)
library(ADvariantExplorer)

etx <- EndophenotypeExplorer$new(NA, NA, initialize.snpLocs=TRUE)

targetGene <- "ADAMTS4"
tag.snp <- "rs4575098"
tag.snp.2 <- "rs2070902"

tbl.hap <- read.table("~/github/ADvariantExplorer/explore/adamts4-study/haploreg-rs4575098-0.2.tsv", sep="\t", as.is=TRUE, header=TRUE)
dim(tbl.hap)
tbl.hap.2 <- read.table("~/github/ADvariantExplorer/explore/adamts4-study/haploreg-rs2070902-0.2.tsv", sep="\t", as.is=TRUE, header=TRUE)
dim(tbl.hap.2)

chrom.loc <- "1"
shoulder <- 10000
start.loc <-  min(tbl.hap$hg38) - shoulder
end.loc   <-  max(tbl.hap$hg38) + shoulder
printf("span: %d", 1 + end.loc-start.loc)


# calculated in test_ADvariantExplorer.R, with these parameters
#    tag.snp <- "rs4575098"
#   tag.snp.loc <- 161185602
#   chrom <- "1"
#   ld.snp  <- "rs11585858"
#    shoulder <- 10000

if(!exists("tbl.eQTLs")){
    tbl.1 <- get(load("tbl.eqtls.1420652x6.RData"))
    tbl.2 <- get(load("tbl.eqtls.25initialFailures.61389x6.RData"))
    tbl.eQTLs <- rbind(tbl.1, tbl.2)
    new.order <- order(tbl.eQTLs$pvalue, decreasing=FALSE)
    tbl.eQTLs <- tbl.eQTLs[new.order,]
    rownames(tbl.eQTLs) <- NULL
    soi <- c(grep("brain", unique(tbl.eQTLs$id), ignore.case=TRUE, v=TRUE),
             grep("microglia", unique(tbl.eQTLs$id), ignore.case=TRUE, v=TRUE),
             grep("macrophage", unique(tbl.eQTLs$id), ignore.case=TRUE, v=TRUE))
    tbl.eqtl.choice <- subset(tbl.eQTLs, id %in% soi)
    dim(tbl.eqtl.choice)
    }

#----------------------------------------------------------------------------------------------------

if(!exists("avx")){
   avx <- ADvariantExplorer$new(targetGene, chrom.loc, start.loc, end.loc)
   }
#----------------------------------------------------------------------------------------------------
gwas.survey <- function()
{
    tbl.gwasAll <- avx$getFullGwasTable()
    dim(tbl.gwasAll)  # 111685 43
    tbl.gwas <- avx$getFilteredGwasTable(targetGeneOnly=FALSE, studyString="alzheimer")
    tbl.gwas <- avx$getFilteredGwasTable(targetGeneOnly=FALSE) #, studyString="alzheimer")

    dim(tbl.gwas)

} # gwas.survey
#----------------------------------------------------------------------------------------------------
viz <- function()
{
   if(!exists("igv")){
      igv <<- start.igv(targetGene, "hg38")
      showGenomicRegion(igv, "chr1:161,068,852-161,253,859")
      }

   tbl.tag.snp <- subset(tbl.hap, rsid==tag.snp)[, c("chrom", "hg38", "hg38", "rsid")]
   colnames(tbl.tag.snp) <- c("chrom", "start", "end", "rsid")
   tbl.tag.snp$start <- tbl.tag.snp$start - 1
   track <- DataFrameAnnotationTrack(tag.snp, tbl.tag.snp, color="red")
   displayTrack(igv, track)

   tbl.hap.snp <- tbl.hap[, c("chrom", "hg38", "hg38", "rSquared")]
   colnames(tbl.hap.snp) <- c("chrom", "start", "end", "score")
   tbl.hap.snp$start <- tbl.hap.snp$start - 1
   track <- DataFrameQuantitativeTrack("haplo.snps R^2", tbl.hap.snp, color="brown", autoscale=TRUE)
   displayTrack(igv, track)

   if(FALSE){
      if(!exists("tbl.eqtl.ndufs2")){
         tbl.eqtls <- get(load("tbl.eqtls.240310x6-10.brain.tissues.RData"))
         tbl.eqtl.ndufs2 <- subset(tbl.eqtls, gene=="NDUFS2" & pvalue < 1e-2)
         dim(tbl.eqtl.ndufs2) # 457 6
         tbl.rsid.locs <- etx$rsidToLoc(unique(tbl.eqtl.ndufs2$rsid))
         tbl.rsid.locs <- unique(tbl.rsid.locs)
         dim(tbl.rsid.locs)
         tbl.eqtl.ndufs2 <- merge(tbl.eqtl.ndufs2, tbl.rsid.locs, by="rsid", all.x=TRUE)
         save(tbl.eqtl.ndufs2, file="tbl.eqtl.ndufs2-10.brainTissues.pval.01.RData")
         as.data.frame(table(tbl.eqtl.ndufs2$id))
           # 1                         GTEx_V8.Brain_Amygdala    1
           # 2   GTEx_V8.Brain_Anterior_cingulate_cortex_BA24   40
           # 3            GTEx_V8.Brain_Cerebellar_Hemisphere  103
           # 4                           GTEx_V8.Brain_Cortex   32
           # 5               GTEx_V8.Brain_Frontal_Cortex_BA9   18
           # 6                      GTEx_V8.Brain_Hippocampus   70
           # 7                     GTEx_V8.Brain_Hypothalamus   35
           # 8  GTEx_V8.Brain_Nucleus_accumbens_basal_ganglia   60
           # 9         GTEx_V8.Brain_Spinal_cord_cervical_c-1   40
           # 10                            ROSMAP.brain_naive   58
          } # if exists
      for(tissue in sort(unique(tbl.eqtl.ndufs2$id))){
          tbl.track <- subset(tbl.eqtl.ndufs2, id==tissue)[, c("chrom", "hg38", "hg38", "pvalue", "rsid")]
          colnames(tbl.track) <- c("chrom", "start", "end", "score", "name")
          tbl.track$start <- tbl.track$start - 1
          tbl.track$score <- -log10(tbl.track$score)
          track.name <- sprintf("NDUFS2 eqtl %s", sub("GTEx_V8.", "", tissue, fixed=TRUE))
          track <- DataFrameQuantitativeTrack(track.name, tbl.track,
                                              autoscale=FALSE, min=0, max=10,
                                              color="random",
                                              trackHeight=30)
          displayTrack(igv, track)
          } # for tissue
       } # manhattan plots for NDUFS2 in multiple tissues

} # viz
#----------------------------------------------------------------------------------------------------
eQTLs <- function()
{
     # first, set the igv view to include all the haplotypes of the tag.snp
   roi <- getGenomicRegion(igv)
   avx <- ADvariantExplorer$new(targetGene, roi$chrom, roi$start, roi$end)

   tbl.cat <- avx$geteQTLSummary()
   studies <- unique(subset(tbl.cat, quant_method=="ge")$unique_id)
   length(studies)  # 157
   studies.brain <- grep("gtex_v8.brain", studies, ignore.case=TRUE, v=TRUE)
   studies.brain <- c("ROSMAP.brain_naive", studies.brain)
   #old.gtex.notWorking <- grep("GTEx.", studies, fixed=TRUE)
   #length(old.gtex.notWorking) # 49
   #studies <- studies[-old.gtex.notWorking]
   length(studies.brain)   # 14

   tbl.eqtls <- avx$geteQTLsByLocationAndStudyID(roi$chrom, roi$start, roi$end,
                                                 studies.brain, method="REST", simplify=TRUE)
   save(tbl.eqtls, file="tbl.eqtls.240310x6-10.brain.tissues.RData")

    strong.hap.rsid <- subset(tbl.hap, rSquared >= 0.5)$rsid
    length(strong.hap.rsid)  # 19

    strong.hap.rsid <- subset(tbl.hap, rSquared >= 0.8)$rsid
    length(strong.hap.rsid)  # 2

    pval.threshold <- 10e-4
    beta.threshold <- 0.1
    tbl.eStrong <-
       subset(tbl.eQTLs, rsid %in% strong.hap.rsid & pvalue <= pval.threshold & abs(beta) > beta.threshold)
    dim(tbl.eStrong) # 104 6
    rownames(tbl.eStrong) <- NULL

    tbl.eStrong <-
       subset(tbl.eqtl.choice, rsid %in% strong.hap.rsid & pvalue <= pval.threshold & abs(beta) > beta.threshold)
    rownames(tbl.eStrong) <- NULL
    dim(tbl.eStrong)  # 46 6
    table(tbl.eStrong$gene)   # ADAMTS4 B4GALT3  FCER1G  NDUFS2
                              #      2      12      19      13
    table(tbl.eStrong$id, tbl.eStrong$gene)

} # eQTLs
#----------------------------------------------------------------------------------------------------
examine.breaks <- function()
{
    f <- "motif.breaks.2119.variantsno.pvals.Sun-Sep-26-10:10:16-2021.RData"
    file.exists(f)
    print(load(f))
    dim(tbl.breaks)   # 519k
    head(tbl.breaks)
    colnames(tbl.breaks)[c(1,6)] <- c("chrom", "rsid")
    tbl.breaks$chrom <- as.character(tbl.breaks$chrom)
    tbl.breaks$pctDelta <- with(tbl.breaks, pctAlt - pctRef)

    head(subset(tbl.breaks, rsid==tag.snp)) # 277

} # examine.breaks
#----------------------------------------------------------------------------------------------------
runFimo <- function()
{
  source("~/github/fimoService/batchMode/fimoBatchTools.R")
  meme.file <- "jaspar2018-hocomocoCore.meme"
  library(MotifDb)
  motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "HOCOMOCOv11-core-A"))

  length(motifs)
  export(motifs, con=meme.file, format="meme")

  tbls.fimo <- list()
  for(goi in c("NDUFS2", "FCER1G", "B4GALT3", "ADAMTS4", "ATF6")){
     printf("--- fimoBatch on %s", goi)
     tbl.gh <- ghRegions[[goi]]
     tbl.fimo <- fimoBatch(tbl.gh, matchThreshold=1e-3, genomeName="hg38", pwmFile=meme.file)
     tbls.fimo[[goi]] <- tbl.fimo
  } # for goi

  save(ghRegions, tbls.fimo, file="fimo.1e-3.RData")

} # runFimo
#----------------------------------------------------------------------------------------------------
# of at least these strong pval, high beta genes, eqtl targets of rs4575098 and hap-sib rs11585858
# down-regulated: NDUFS2, FCER1G, B4GALT3
# up-regulated: ATF6
build.models <- function()
{
   igv <- start.igv("NDUFS2", "hg38")
   library(ghdb)
   ghdb <- GeneHancerDB()

   ghRegions <- list()
   for(goi in c("NDUFS2", "FCER1G", "B4GALT3", "ADAMTS4", "ATF6")){
      tbl.gh <- retrieveEnhancersFromDatabase(ghdb, goi, tissues="all")
      tbl.gh$score <- asinh(tbl.gh$combinedscore)
      trackName <- sprintf("GH %s", goi)
      track <- DataFrameQuantitativeTrack(trackName, tbl.gh[, c("chrom", "start", "end", "score")],
                                       autoscale=TRUE, color="brown")
      #displayTrack(igv, track)
      ghRegions[[goi]] <- tbl.gh
      } # for goi



} # build.models
#----------------------------------------------------------------------------------------------------
