library(TrenaMultiScore)
library(TrenaProjectAD)
library(EndophenotypeExplorer)
library(ADvariantExplorer)
library(RUnit)
library(trena)
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
#----------------------------------------------------------------------------------------------------
targetGene <- "NDUFS2"
etx <- EndophenotypeExplorer$new(targetGene, "hg38", initialize.snpLocs=TRUE)
tpad <- TrenaProjectAD()
tbl.eqtl <- get(load("../tbl.eqtls.adamts4Study.RData"))
tbl.fimo <- get(load("tbl.fimo.NDUFS2.RData"))
tbl.oc <- get(load("~/github/TrenaProjectAD/inst/extdata/genomicRegions/boca-hg38-consensus-ATAC.RData"))
tbl.tagHap <- read.table("~/github/ADvariantExplorer/explore/adamts4-study/haploreg-rs4575098-0.2.tsv",
                          sep="\t", as.is=TRUE, header=TRUE)
dim(tbl.tagHap)

avx <- ADvariantExplorer$new(targetGene, "chr1", 161198314, 161215395)
#----------------------------------------------------------------------------------------------------
createTMS <- function()
{
    tms <- TrenaMultiScore(tpad, targetGene, tbl.fimo, tbl.oc, quiet=FALSE)

    wdth(100)
    scoreMotifHitsForOpenChromatin(tms)
    scoreMotifHitsForConservation(tms)
    scoreMotifHitsForGeneHancer(tms)
    addDistanceToTSS(tms)
    addGenicAnnotations(tms)
    addChIP(tms)
    tbl.tms <- getMultiScoreTable(tms)

    names(etx$get.rna.matrix.codes())

    code <- "max-tcx"
    mtx.rna <- etx$get.rna.matrix(code)
    addGeneExpressionCorrelations(tms, mtx.rna, featureName=code, method="spearman")
    tbl.tms <- getMultiScoreTable(tms)
    fivenum(tbl.tms[, code])    # -0.64 -0.31 -0.12  0.10  0.67

    code <- "max-rosmap"
    mtx.rna <- etx$get.rna.matrix(code)
    addGeneExpressionCorrelations(tms, mtx.rna, featureName=code, method="spearman")
    tbl.tms <- getMultiScoreTable(tms)
    fivenum(tbl.tms[, code])    #  -0.39 -0.14 -0.04  0.05  0.31

    code <- "gtex.v8.Brain_Cortex"
    mtx.rna <- etx$get.rna.matrix(code)
    addGeneExpressionCorrelations(tms, mtx.rna, featureName=code, method="spearman")
    tbl.tms <- getMultiScoreTable(tms)
    fivenum(tbl.tms[, code])    #  -0.67 -0.40 -0.14  0.15  0.75

    dim(tbl.tms)  # 817786 18  -0.67 -0.40 -0.14  0.15  0.75

    dim(subset(tbl.tms, chip & gh > 5 & fimo_pvalue < 1e-5 & abs(gtex.v8.Brain_Cortex) > 0.3))

    checkTrue(mean(tbl.tms$gh) > 0)
    table(tbl.tms$oc)
    table(tbl.tms$chip)

    dim(tbl.eqtl)
    tbl.eqtl.01 <- subset(tbl.eqtl, gene==targetGene & grepl("brain", study, ignore.case=TRUE))
    dim(tbl.eqtl.01)   # 4053 8
    table(tbl.eqtl.01$study)

    rsids <- subset(tbl.eqtl.01, pvalue <= 0.00001 & abs(beta) >= 0.15)$rsid
    length(rsids)

    add.eQTLs(tms, tbl.eqtl.01, pval=0.01, abs.beta=0.15, eqtl.title="brain.eqtl")


    tbl.tms <- getMultiScoreTable(tms)
    dim(tbl.tms)
    tbl.tms.sub <- subset(tbl.tms, (gh > 0 | oc) & brain.eqtl > 5 & abs(gtex.v8.Brain_Cortex) > 0.4 & fimo_pvalue < 1e-3)
    tbl.tms.sub <- subset(tbl.tms, (chip & oc) | brain.eqtl > 5)
    tbl.tms.sub <- subset(tbl.tms, brain.eqtl > 5 & oc)
    dim(tbl.tms.sub) # 104 19
    tfs <- tbl.tms.sub$tf
    length(unique(tfs)) # 51
    table(tbl.tms.sub$tf)
    tfs.unique <- unique(tbl.tms.sub$tf)
    length(tfs.unique) # 51
    save(tbl.tms, file="tbl.tms.21dec2021.817786x21.RData")


    #tbl.tsne <- unique(tbl.tms.sub[, c("oc", "old.rosmap", "gh", "brain.eqtl", "tf")])
    #tfs <- tbl.tsne$tf
    #length(unique(tfs))
    #tsne.out <- Rtsne(tbl.tsne, perplexity=5)
    #plot(tsne.out$Y, xlim=c(-70, 50), ylim=c(-90, 70))
    # text(jitter(tsne.out$Y, factor=5, amount=20), labels=tfs)

} # createTMS
#----------------------------------------------------------------------------------------------------
runTrena <- function()
{
   noquote(names(etx$get.rna.matrix.codes()))

   codes <- names(etx$get.rna.matrix.codes())
   keepers <- c(grep("gtex", codes), grep("^max", codes))
   codes <- codes[keepers]
   deleters <- grep("Spinal", codes)
   codes <- codes[-deleters]
   tbls.trena <- list()

   for(code in codes){
      printf("-------------- code: %s", code)
      mtx.rna <- etx$get.rna.matrix(code)
      mtx.rna[is.na(mtx.rna)] <- 0
      dim(mtx.rna)
      var <- apply(mtx.rna, 1, var)
      deleters <- which(is.na(var))
      length(deleters)
      if(length(deleters) > 0)
          mtx.rna <- mtx.rna[-deleters,]
      solver <- EnsembleSolver(mtx.rna,
                               targetGene=targetGene,
                               candidateRegulators=tfs.unique,
                               solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"),
                               geneCutoff=1.0)
      tbl.out <- run(solver)
      new.order <- order(abs(tbl.out$spearman), decreasing=TRUE)
      head(tbl.out[new.order,], n=20)
      dim(tbl.out)
      tfbs.counts <- unlist(lapply(tbl.out$gene, function(gene) nrow(subset(tbl.tms.sub, tf==gene))))
      tbl.out$tfbs <- tfbs.counts
      tbl.out <- tbl.out[order(tbl.out$rfScore, decreasing=TRUE),]
      rownames(tbl.out) <- NULL
      for(cname in colnames(tbl.out)[2:7]) tbl.out[, cname] <- round(tbl.out[, cname], digits=2)
      tbl.final <- head(tbl.out, n=50)
      goi <- c(targetGene, tbl.final$gene)
      length(goi)
      tbl.ctAll <- avx$getCellTypes()
      ct <- unlist(lapply(goi, function(g)
          paste(sort(unique(subset(tbl.ctAll, gene==g & p_val_adj < 0.05)$ct)), collapse=",")))
      target.gene.ct <- ct[1]
      tbl.final$cellTypes <- ct[2:51]
      tbl.final$mtx <- code
      tbl.final$rank <- seq_len(nrow(tbl.final))
      tbls.trena[[code]] <- tbl.final
      }
   tbl.trena <- do.call(rbind, tbls.trena)
   deleters <- which(nchar(tbl.trena$cellTypes) == 0)
   tbl.trena <- tbl.trena[-deleters,]
   rownames(tbl.trena) <- NULL
   dim(tbl.trena)
   tbl.tfs <- as.data.frame(sort(table(tbl.trena$gene), decreasing=TRUE))
   tfs.main <- names(sort(table(tbl.trena$gene), decreasing=TRUE))

   tbl.trena.strong <- subset(tbl.trena, rank <= 10)
   tbl.tf.counts <- head(as.data.frame(sort(table(tbl.trena.strong$gene), decreasing=TRUE)), n=20)
   tbl.tf.counts$Var1 <- as.character(tbl.tf.counts$Var1)



} # runTrena
#----------------------------------------------------------------------------------------------------
break.motifs <- function()
{
   snps <- c("rs4575098", "rs11585858")
   snps <- tbl.tagHap$rsid

   snps.gr <- snps.from.rsid(rsid=snps,
                              dbSNP=SNPlocs.Hsapiens.dbSNP151.GRCh38,
                              search.genome=BSgenome.Hsapiens.UCSC.hg38)
   tfs.oi <- tbl.tf.counts[,1]
   mdb <- query(MotifDb, "Hsapiens", c("jaspar2018", "hocomoco-core-A"))
   mdb.indices <- match(tfs.oi, mcols(mdb)$geneSymbol)
   length(mdb.indices)
   mdb.indices <- mdb.indices[-which(is.na(mdb.indices))]
   length(mdb.indices)  # 18

   motifs <- mdb[mdb.indices]
   length(motifs)

   bpparam <- MulticoreParam(workers=4)
   results <- motifbreakR(snpList = snps.gr,
                           filterp = TRUE,
                           pwmList = motifs,
                           show.neutral=FALSE,
                           method = c("ic", "log", "notrans")[1],
                           bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                           BPPARAM = bpparam,
                           verbose=TRUE)
   tbl.breaks <- as.data.frame(results, row.names=NULL)
   tbl.breaks$pctDelta <- with(tbl.breaks, pctAlt - pctRef)
   new.order <- order(abs(tbl.breaks$pctDelta), decreasing=TRUE)
   tbl.breaks <- tbl.breaks[new.order,]

   coi <- c(6:8,5,  11,12,13, 16, 17, 25)
   save(tbl.tms, tbl.trena, tbl.breaks, file="tms.trena.breaks-15dec21.RData")

  # unique(subset(tbl.breaks, geneSymbol %in% head(tbl.out$gene))$SNP_id)


} # break.motifs
#----------------------------------------------------------------------------------------------------
viz <- function()
{
   igv <- start.igv(targetGene, "hg38")
   zoomOut(igv)
   zoomOut(igv)
   require(ghdb)
   ghdb <- GeneHancerDB()
   tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
   tbl.gh$score <- asinh(tbl.gh$combinedscore)
   track <- DataFrameQuantitativeTrack("GH", tbl.gh[, c("chrom", "start", "end", "score")],
                                       autoscale=TRUE, color="brown")
   displayTrack(igv, track)

   tbl.locs <- etx$rsidToLoc(c("rs4575098", "rs11585858"))
   tbl.track <- tbl.locs[, c("chrom", "hg38", "hg38", "rsid")]
   colnames(tbl.track) <- c("chrom", "start", "end", "rsid")
   tbl.track$start <- tbl.track$start - 1
   track <- DataFrameAnnotationTrack("tag.snp", tbl.track, color="red")
   displayTrack(igv, track)

   tbl.eqtls <- get(load("../tbl.eqtls.adamts4Study.RData"))
   dim(tbl.eqtls)    # 1304426       9

   tbl.ndufs2.eqtl.smallSample <- subset(tbl.eqtl, hg38 >= roi$start & hg38 <= roi$end &
                                                   gene=="NDUFS2" & pvalue < 0.01 &
                                                   grepl("brain", study, ignore.case=TRUE))
   save(tbl.ndufs2.eqtl.smallSample,
        file="~/github/TrenaMultiScore/inst/extdata/tbl.ndufs2.eqtl.smallSample.RData")
   tbl.track <- tbl.ndufs2.eqtl.smallSample[, c("chrom", "hg38", "hg38", "pvalue", "rsid", "study")]
   colnames(tbl.track)[1:4] <- c("chrom", "start", "end", "score")
   tbl.track$score <- -log10(tbl.track$score)
     # lots of duplicte reports, with different pval by study.  get the best, eliminate the others
   tbl.track <- tbl.track[order(tbl.track$score, decreasing=TRUE),]
   dups <- which(duplicated(tbl.track[, 1:3]))
   tbl.track <- tbl.track[-dups,]
   dim(tbl.track)  # 13 5
   track <- DataFrameQuantitativeTrack("NDUFS2 eQTLS", tbl.track, autoscale=TRUE, color="blue")
   displayTrack(igv, track)

} # viz
#----------------------------------------------------------------------------------------------------

