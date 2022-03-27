#' @title ADvariantExplorer
#' @description access to GWAS and eQTL data, slighly optimized for Alzheimer's disease
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
        #' by location and/or targetGene
        #' @param targetGeneOnly logical, default TRUE, otherwise all variants in this regeion,
        #' @param studyString character, default "alzheimer", possibly empty string or any other
        #'   study field matchabel substring
        #' affecting any gene whatsoever
        #'
        #' @return a data.frame
        getFilteredGwasTable = function(targetGeneOnly=TRUE, studyString="alzheimer", trim.columns=TRUE){
            printf("--- entering getFilteredGwasTable, targetGeneOnly = %s", targetGeneOnly)
            tbl.sub <- subset(private$tbl.gwascat,
                              chrom==sub("chr", "", private$loc.chrom) &
                              start >= private$loc.start &
                              end <= private$loc.end)
            if(targetGeneOnly){
               keepers <- grep(private$targetGene, tbl.sub$MAPPED_GENE)
               tbl.sub <- tbl.sub[keepers,]
               }
            if(nchar(studyString) > 0){
               soi <- grep(studyString, tbl.sub$STUDY, ignore.case=TRUE)
               tbl.sub <- tbl.sub[soi,]
               }
            if(trim.columns){
               coi <- c("SNPS", "P.VALUE", "MAPPED_GENE", "FIRST.AUTHOR", "DATE", "INITIAL.SAMPLE.SIZE", "STUDY" )
               tbl <- tbl.sub[, coi]
               tbl$STUDY <- substr(tbl$STUDY, 1, 30)
               tbl$INITIAL.SAMPLE.SIZE <- substr(tbl$INITIAL.SAMPLE.SIZE, 1, 20)
               }
            new.order <- order(tbl$P.VALUE, decreasing=FALSE)
            tbl <- tbl[new.order,]
            rownames(tbl) <- NULL
            private$tbl.gwascat.filtered <- tbl
            tbl
            },

        #------------------------------------------------------------
        #' @description access to the complete EMBL-EBI GWAS catalog via bioc
        #' @param trim.columns logical reduce to a legible 7 columns
        #' @return a data.frame
        getFullGwasTable = function(trim.columns=FALSE){
            tbl.full <- private$tbl.gwascat
            printf("--- getFullGwasTable, trim.columns: %s", trim.columns)
            if(!trim.columns){
               return(invisible(tbl.full))
               }
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
            as.data.frame(private$tbl.eqtl.summary)
            },

        #------------------------------------------------------------
        #' @description returns the retrieval-ready study names given
        #' @param groupMatchingString character, (part of) a qtl_group string
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
        geteQTLsByLocationAndStudyID = function(chrom, start, end, studyIDs, simplify=FALSE){
           #method <- tolower(method)
           #stopifnot(method %in% c("rest", "tabix"))
           tbls <- list()
           for(id in studyIDs){
              message(sprintf("--- fetching %s (ge)", id))
              tryCatch({
                  #if(method=="rest"){
                      suppressWarnings({tbl <- eQTL_Catalogue.fetch(unique_id=id,
                                                                    quant_method="ge",
                                                                    method="REST",
                                                                    chrom = sub("chr", "", chrom),
                                                                    bp_lower=start,
                                                                    bp_upper=end,
                                                                    verbose=TRUE)})
                  #} # rest
                  #else{ # tabix
                  #    suppressWarnings({tbl <- eQTL_Catalogue.fetch(unique_id=id,
                  #                                                  quant_method="ge",
                  #                                                  nThread = 1,
                  #                                                  use_tabix=TRUE,
                  #                                                  chrom = sub("chr", "", chrom),
                  #                                                  bp_lower=start,
                  #                                                  bp_upper=end,
                  #                                                  verbose=TRUE)})
                  #} # tabix
                  tbl$id <- id
                  tbls[[id]] <- tbl
              },
              error = function(e){
                  message(sprintf("eQTL_Catalogue.fetch failed on study %s", id))
                  print(e)
                  })
              } # for id
           tbl.out <- do.call(rbind.fill, tbls)
           rownames(tbl.out) <- NULL
           if(is.null(tbl.out))
               return(data.frame())
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
           }, # geteEQTLsByLocationAndCategory

        #------------------------------------------------------------------
        #' @description returns one gene's cell types, or the entire table
        #' from jiang et al 2020, nine cell types across 10 human (and mouse) brain regions:
        #' https://www.sciencedirect.com/science/article/pii/S2589004220309664
        #' scREAD: A Single-Cell RNA-Seq Database for Alzheimer's Disease
        #' @param geneSymbol character string, either a HUGO gene symbol, or NA
        #' @return a character vector or a data.frame (if geneSymbol param is NA)
        getCellTypes=function(geneSymbol=NA){

            tbl.cellTypes <- get(load(system.file(package="ADvariantExplorer", "extdata",
                                                  "scRead.tbl.celltypes.3276904x12.RData")))
            if(!is.na(geneSymbol))
               return(subset(tbl.cellTypes, gene==geneSymbol))

            invisible(tbl.cellTypes)
            }

       ) # public

    ) # class
#--------------------------------------------------------------------------------
