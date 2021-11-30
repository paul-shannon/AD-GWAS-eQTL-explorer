#' @title ADvariantExplorer
#' @description A template for building documented, tested R6 classes
#' @name ADvariantExplorer

#' @import catalogueR
#' @import gwascat
#' @import EnsDb.Hsapiens.v79
#' @import AnnotationDbi
#' @import plyr
#
#' @export

ADvariantExplorer = R6Class("ADvariantExplorer",

    #--------------------------------------------------------------------------------
    private = list(targetGene=NULL,
                   loc.chrom=NULL,
                   loc.start=NULL,
                   loc.end=NULL,
                   tbl.gwascat=NULL,
                   tbl.gwascat.filtered=NULL,
                   tag.snps=NULL,
                   tbl.haplotypes=NULL,
                   tbl.eqtl.summary=NULL),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param targetGene character, a candidate gene for this AD GWAS locus
         #' @param loc.chrom character, the chromosome of the region to consider
         #' @param loc.start numeric, the start of the region to consider
         #' @param loc.end   numeric, the end of the region to consider
         #' @return a new instance of ADvariantExplorer
        initialize = function(targetGene, loc.chrom, loc.start, loc.end){
            private$targetGene <- targetGene
            private$loc.chrom <- loc.chrom
            private$loc.start <- loc.start
            private$loc.end   <- loc.end
            f <- system.file(package="ADvariantExplorer", "extdata", "gwascat-31oct2021.RData")
            gr.gwascat <- get(load(f))
            tbl.gwascat <- as.data.frame(gr.gwascat)
            colnames(tbl.gwascat)[1] <- "chrom"
            tbl.gwascat$chrom <- as.character(tbl.gwascat$chrom)
            private$tbl.gwascat <- tbl.gwascat
            private$tbl.eqtl.summary <- get(data(meta, package="catalogueR"))
            },

        #------------------------------------------------------------
        #' @description accessor for the object's targetGene field
        #' @return the current value of the targetGene member
        getTargetGene = function(){
            private$targetGene
        },

        #------------------------------------------------------------
        #' @description access to the EMBL-EBI GWAS catalog via bioc, filtered
        #'  by location and/or targetGene
        #' @return a data.frame
        getFilteredGwasTable = function(){
            tbl.sub <- subset(private$tbl.gwascat, chrom==sub("chr", "", private$loc.chrom) &
                                                   start >= private$loc.start & end <= private$loc.end |
                              grepl("CLU", MAPPED_GENE))
            soi <- grep("Alzheimer", tbl.sub$STUDY, ignore.case=TRUE)
            tbl <- tbl.sub[soi,]
            coi <- c("SNPS", "P.VALUE", "MAPPED_GENE", "FIRST.AUTHOR", "DATE", "INITIAL.SAMPLE.SIZE", "STUDY" )
            tbl <- tbl[, coi]
            tbl$STUDY <- substr(tbl$STUDY, 1, 30)
            tbl$INITIAL.SAMPLE.SIZE <- substr(tbl$INITIAL.SAMPLE.SIZE, 1, 20)
            new.order <- order(tbl$P.VALUE, decreasing=FALSE)
            tbl <- tbl[new.order,]
            rownames(tbl) <- NULL
            private$tbl.gwascat.filtered
            tbl
            },

        #------------------------------------------------------------
        #' @description access to the complete EMBL-EBI GWAS catalog via bioc
        #' @param trim.columns logical reduce to a legible 7 columns
        #' @return a data.frame
        getFullGwasTable = function(trim.columns=FALSE){
            tbl.full <- private$tbl.gwascat
            if(!trim.columns)
                invisible(tbl.full)
            coi <- c("SNPS", "P.VALUE", "MAPPED_GENE", "FIRST.AUTHOR", "DATE", "INITIAL.SAMPLE.SIZE", "STUDY" )
            tbl <- tbl.full[, coi]
            tbl$STUDY <- substr(tbl$STUDY, 1, 30)
            tbl$INITIAL.SAMPLE.SIZE <- substr(tbl$INITIAL.SAMPLE.SIZE, 1, 20)
            invisible(tbl)
            },

        #------------------------------------------------------------
        #' @description access to the complete EMBL-EBI eQTL Catalogue
        #' via the catalogueR R package, focusing on expression cis-QTLs
        #' and sQTL (splicing)
        #' @return a data.frame
        geteQTLSummary = function(){
            private$tbl.eqtl.summary
            },

        #------------------------------------------------------------
        #' @description returns the retrieval-ready study names given
        #' (part of) a qtl_group string
        #' @return a data.frame
        geteqtlStudyNamesForGroup = function(groupMatchingString){
            tbl.cat <- private$tbl.eqtl.summary
            indices <- grep(groupMatchingString, tbl.cat$qtl_group)
            if(length(indices) == 0)
                return(NA)
            sort(unique(tbl.cat[indices,]$unique_id))
            },

        #------------------------------------------------------------
        #' @description access to the complete EMBL-EBI eQTL Catalogue
        #' via the catalogueR R package, focusing on expression cis-QTLs
        #' and sQTL (splicing
        #' @param chrom character
        #' @param start numeric
        #' @param end numeric
        #' @param studyIDs character vector, one or more values from eQTL summary unique_id column
        #' @param method character, either REST or tabix, REST by default
        #' @param simplify logical, trims columns to just the crucial 4: rsid, pvalue, gene, samples
        #' @return a data.frame ordered by increasing pvalue.QTL
        geteQTLsByLocationAndStudyID = function(chrom, start, end, studyIDs, method="REST", simplify=FALSE){
           method <- tolower(method)
           stopifnot(method %in% c("rest", "tabix"))
           tbls <- list()
           for(id in studyIDs){
              message(sprintf("--- fetching %s (ge)", id))
              if(method=="rest"){
                  suppressWarnings({tbl <- fetch_restAPI(unique_id=id,
                                                         quant_method="ge",
                                                         chrom = sub("chr", "", chrom),
                                                         bp_lower=start,
                                                         bp_upper=end,
                                                         verbose=TRUE)})
                  } # rest
              else{ # tabix
                 suppressWarnings({tbl <- eQTL_Catalogue.fetch(unique_id=id,
                                                               quant_method="ge",
                                                               nThread = 1,
                                                               use_tabix=TRUE,
                                                               chrom = sub("chr", "", chrom),
                                                               bp_lower=start,
                                                               bp_upper=end,
                                                               verbose=TRUE)})
                 } # tabix
              tbl$id <- id
              tbls[[id]] <- tbl
              } # for id
           tbl.out <- do.call(rbind.fill, tbls)
           rownames(tbl.out) <- NULL
           new.order <- order(tbl.out$pvalue.QTL, decreasing=FALSE)
           tbl.out <- tbl.out[new.order,]
           coi <- c("rsid.QTL", "pvalue.QTL", "gene_id.QTL", "an.QTL", "beta.QTL", "id")
           if(simplify){
              tbl.out <- tbl.out[, coi]
              colnames(tbl.out) <- c("rsid", "pvalue", "gene", "total.alleles", "beta", "id")
              }
           map <- mapIds(EnsDb.Hsapiens.v79, tbl.out$gene, "SYMBOL", "GENEID")
           tbl.map <- data.frame(ensg=names(map), symbol=as.character(map), stringsAsFactors=FALSE)
           na.indices <- which(is.na(tbl.map$symbol))
           length(na.indices)
           tbl.map$symbol[na.indices] <- tbl.map$ensg[na.indices]
           tbl.out$gene <- tbl.map$symbol
           rownames(tbl.out) <- NULL
           invisible(as.data.frame(tbl.out))
           } # geteEQTLsByLocationAndCategory


       ) # public

    ) # class
#--------------------------------------------------------------------------------
