library(ADvariantExplorer)
library(EndophenotypeExplorer)

if(!exists("state"))
    state <- new.env(parent=emptyenv())

targetGene <- "CLU"

collab.rsids <- c("rs73223431", "rs28834970", "rs2322599")
#----------------------------------------------------------------------------------------------------
gwas.lookup.collaborator.variants <- function()
{
    loc.chrom <- "chr8"
    loc.start <- 27447528
    loc.end   <- 27764088
    span <- 1 + loc.end - loc.start
    avx <- ADvariantExplorer$new(targetGene, loc.chrom, loc.start, loc.end)
    tbl.f <- avx$getFilteredGwasTable()
    gwas.rsids <- unique(tbl.f$SNPS)
    printf("%d variants in %d kb", length(gwas.rsids), span)
    intersect(collab.rsids, gwas.rsids)

    tbl.all <- avx$getFullGwasTable(trim.columns=TRUE)
    dim(tbl.all)
    tbl.collab <- subset(tbl.all, SNPS %in% collab.rsids)  # 9 7
    new.order <- order(tbl.collab$P.VALUE, decreasing=FALSE)
    tbl.collab <- tbl.collab[new.order,]
    rownames(tbl.collab) <- NULL

} # gwas.lookup.collaborator.variants
#----------------------------------------------------------------------------------------------------
eqtl.lookup.collaborator.variants <- function()
{
    etx <- EndophenotypeExplorer$new(targetGene, "hg38", initialize.snpLocs=TRUE)
    tbl.locs <- etx$rsidToLoc(collab.rsids)
    printf("span: %d", max(tbl.locs$hg38) - min(tbl.locs$hg38))  # 25k
    tbl.cat <- avx$geteQTLSummary()
    dim(tbl.cat)

    studies <- subset(tbl.cat, tissue_label %in% c("brain (DLPFC)", "brain (cortex)") & quant_method=="ge")$unique_id
    gtex <- grep("GTEx", studies)
    studies <- studies[-gtex]
      #  "BrainSeq.brain"    "ROSMAP.brain_naive"   "GTEx.brain_cortex"  "GTEx.brain_frontal_cortex"

    tbls <- list()
    for(r in seq_len(nrow(tbl.locs))){
       chrom <- tbl.locs$chrom[r]
       start <- tbl.locs$hg38[r] - 1
       end <- tbl.locs$hg38[r]
       tbl <- avx$geteQTLsByLocationAndStudyID(chrom, start, end, studies, method="REST", simplify=TRUE)
       tbls[[r]] <- tbl
       }


} # eqtl.lookup.collaborator.variants
#----------------------------------------------------------------------------------------------------
viz <- function()
{
   if(!exists("igv")) {
     igv <- start.igv("CLU", "hg38")
     showGenomicRegion(igv, "chr8:27,274,766-27,723,224")  # from PTK2B to SCARA3
     }


} # viz
#----------------------------------------------------------------------------------------------------
