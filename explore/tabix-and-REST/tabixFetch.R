# tabix ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/GTEx/ge/GTEx_ge_brain_frontal_cortex.all.tsv.gz 8:27610984-27610987

tissue <- "GTEx_ge_brain_frontal_cortex"
chrom <- "8"
start <- 27610984
end <- 27610987

ftp.uri.base <- "ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/GTEx/ge"
cmd <- sprintf("tabix %s/%s.all.tsv.gz %s:%d-%d", ftp.uri.base, tissue, chrom, start, end)

text  <- system(cmd, intern=TRUE)
tmp.filename <- tempfile()
writeLines(text, con=tmp.filename)
tbl <- read.table(tmp.filename, sep="\t", as.is=TRUE)
dim(tbl)

colnames(tbl) <- c("ensg.1", "chrom", "hg38", "ref", "alt", "loc", "samples",
                   "q1", "q2", "pval", "beta", "type", "q3", "q4", "tpm",
                   "ensg.2", "ensg.3", "q5", "rsid")
