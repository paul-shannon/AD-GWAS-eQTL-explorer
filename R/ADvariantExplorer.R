#' @title ADvariantExplorer
#' @description A template for building documented, tested R6 classes
#' @name ADvariantExplorer

#' @import catalogueR
#' @import gwascat
#'
#' @examples
#'   rt <- R6Template$new(id="abc")
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
                   tbl.haplotypes=NULL
                   ),

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
            }


       ) # public

    ) # class
#--------------------------------------------------------------------------------
