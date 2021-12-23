f1 <- "tbl.eqtls.1420652x6.RData"
f2 <- "tbl.eqtls.25initialFailures.61389x6.RData"
tbl.1 <- get(load(f1))
tbl.2 <- get(load(f2))
tbl.eqtl <- rbind(tbl.1, tbl.2)
dim(tbl.eqtl) # [1] 1482041       6
pval <- 0.5
betaVal <- 0.05
length(unique(subset(tbl.eqtl, gene=="NDUFS2" & pvalue < pval & abs(beta) > betaVal)$rsid))

snps <- unique(subset(tbl.eqtl, gene=="NDUFS2" & pvalue < pval & abs(beta) > betaVal)$rsid)
length(snps) # 584
etx <- EndophenotypeExplorer$new("NDUFS2", "hg38", initialize.snpLocs=TRUE)
tbl.locs <- etx$rsidToLoc(snps)
dim(tbl.locs) # 497 4
save(tbl.locs, file=sprintf("tbl.locs.%d.eqtlVariants.pval.%f.beta.%f.RData", length(snps), pval, betaVal))
start <- min(tbl.locs$hg38)
end   <- max(tbl.locs$hg38)
span <- 1 + end - start
printf("span: %5.2fk", span/1000)   # 142.24k

igv <- start.igv("NDUFS2", "hg38")
library(ghdb)
ghdb <- GeneHancerDB()
goi <- "NDUFS2"
tbl.gh <- retrieveEnhancersFromDatabase(ghdb, goi, tissues="all")
tbl.gh$score <- (tbl.gh$combinedscore)
track <- DataFrameQuantitativeTrack("GH", tbl.gh[, c("chrom", "start", "end", "score")],
                                    autoscale=TRUE, color="brown")
displayTrack(igv, track)
start <- min(tbl.gh$start)
end <- max(tbl.gh$end)
span <- 1 + end - start
printf("gh span: %5.2fk", span/1000)  # 376k

tbl.track <- tbl.locs[, c("chrom", "hg38", "hg38", "rsid")]
colnames(tbl.track) <- c("chrom", "start", "end", "rsid")
tbl.track$start <- tbl.track$start - 1
track <- DataFrameAnnotationTrack("eqtl", tbl.track, color="black")
displayTrack(igv, track)
