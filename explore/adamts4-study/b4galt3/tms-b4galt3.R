library(TrenaMultiScore)
library(EndophenotypeExplorer)
#----------------------------------------------------------------------------------------------------
targetGene <- "B4GALT3"
#----------------------------------------------------------------------------------------------------
runTMS <- function(goi)
{
   tbl.fimo <- get(load(sprintf("tbl.fimo.%s.RData", targetGene)))
   dim(tbl.fimo)  # 1230563       9
   tbl.oc   <- human.brain.mayo.oc()
   trenaProject <- TrenaProjectAD()
   mtx.rna <- get.rna.matrix("old-mayo-tcx")
   corf <- function(gene){
      cor(mtx.rna["NDUFS2",], mtx.rna[gene,], method="spearman", use="pairwise.complete")
      }

   tbl.tfs <- get(load("~/github/MotifDb/inst/extdata/tfs-1683-lambert.RData"))
   head(tbl.tfs)
   tfs.with.rna <- unique(intersect(rownames(mtx.rna), tbl.tfs$Gene))
   length(tfs.with.rna)  # 1436
   x <- lapply(tfs.with.rna, corf)
   names(x) <- tfs.with.rna
   fivenum(as.numeric(x))

   dim(mtx.rna)

   tms <- TMS$new(trenaProject, targetGene, tbl.fimo, tbl.oc, quiet=FALSE)
   tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
   tms$add.tf.mrna.correlations(mtx.rna, featureName="cor.all")
   tbl.tms <- tms$getTfTable()
   fivenum(tbl.tms$cor.all)
   checkEquals(nrow(tbl.tms), nrow(tbl.fimo))

   dim(subset(tbl.tms, fimo_pvalue <= 1e-6 & chip))

   dim(subset(tbl.tms, fimo_pvalue < 1e-5 & abs(cor.all) > 0.8 & chip & gh > 600))
   dim(subset(tbl.tms, fimo_pvalue < 1e-3 & abs(cor.all) > 0.4 & oc))

   tfs <- subset(tbl.tms, fimo_pvalue < 1e-3 & abs(cor.all) > 0.3 & oc)$tf
   tfs <- unique(tfs)
   tfs <- tfs.with.rna
   printf("running trena with candidate tfs: %d", length(unique(tfs)))

   tms$addRBP()
   tms$add.rbp.mrna.correlations(mtx.rna.asinh, featureName="cor.all")   # added to tbl.rbp
   tbl.rbp <- tms$getRbpTable()
   dim(tbl.rbp)

   rbps <- unique(subset(tbl.rbp, abs(cor.all) > 0.3)$gene)
   printf("candidate rbps: %d", length(rbps))
   tbl.trena.tf <- tms$build.trena.model(tfs, list(), mtx.rna.asinh)
   tbl.trena.both <- tms$build.trena.model(tfs, rbps, mtx.rna)
   tbl.trena.tf$target <- targetGene
   tbl.trena.both$target <- targetGene
   save(tbl.trena.tf, tbl.trena.both, file=sprintf("%s/trena.%s.RData", targetGene, targetGene))

} # runTMS
#------------------------------------------------------------------------------------------------------------------------
