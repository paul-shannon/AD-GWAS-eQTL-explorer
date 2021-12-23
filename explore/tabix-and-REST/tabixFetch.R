# tabix ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/GTEx/ge/GTEx_ge_brain_frontal_cortex.all.tsv.gz 8:27610984-27610987

tissue <- "GTEx_ge_brain_frontal_cortex"
chrom <- "8"
start <- 27610984
end <- 27610987

ftp.uri.base <- "ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/GTEx/ge"
cmd <- sprintf("tabix %s/%s.all.tsv.gz %s:%d-%d", ftp.uri.base, tissue, chrom, start, end)

tbl <- system(cmd, intern=TRUE)
