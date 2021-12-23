library(EndophenotypeExplorer)
f1 <- "tbl.eqtls.1420652x6.RData"
f2 <- "tbl.eqtls.25initialFailures.61389x6.RData"
tbl.1 <- get(load(f1))
tbl.2 <- get(load(f2))
tbl.eqtl <- rbind(tbl.1, tbl.2)
dim(tbl.eqtl) # 1482041       6
head(tbl.eqtl)
snps <- unique(subset(tbl.eqtl, gene=="NDUFS2" & pvalue < 0.05 & abs(beta) > 0.1)$rsid)
length(snps) # 378


etx <- EndophenotypeExplorer$new("NDUFS2", "hg38", initialize.snpLocs=TRUE)
tbl.locs <- etx$rsidToLoc(snps)
