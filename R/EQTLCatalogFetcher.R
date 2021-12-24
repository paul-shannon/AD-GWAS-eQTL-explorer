#' @title EQTLCatalogFetcher
#' @description a small class to fetch eqtl regions by study id from the EBI eQTL catalog
#' @name EQTLCatalogFetcher

#'
#' @examples
#'   studyIDs <- "GTEx_ge_brain_frontal_cortex"
#'   chrom <- "8"
#'   start <- 27610984
#'   end <- 27610987
#'
#'   fetcher <- EQTLCatalogFetcher$new(studyIDs, chrom, start, end, simplify=TURE)
#'
#' @export

EQTLCatalogFetcher = R6Class("EQTLCatalogFetcher",

    #--------------------------------------------------------------------------------
    private = list(studyIDs=NULL,
                   chrom=NULL,
                   start=NULL,
                   end=NULL,
                   simplify=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param studyIDs character vector, e.g. "GTEx_ge_brain_frontal_cortex"
         #' @param chrom character, chromosome identifier, with or without leading "chr"
         #' @param start numeric starting base pair of the region of interest
         #' @param end numeric ending base pair of the region of interest
         #' @parma simplify logical, reduce column names to small tractable set
         #' @return a new instance of EQTLCatalogFetcher
        initialize = function(studyIDs, chrom, start, end, simplify=FALSE){
            private$studyIDs <- studyIDs
            private$chrom <- chrom
            private$start <- start
            private$end <- end
            private$simplify <- simplify
            },
        #------------------------------------------------------------
        #' @description obtain the requested slice of eQTLs from each study
        #' @return nothing
        fetch = function(){
           ftp.uri.base <- "ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/GTEx/ge"
           cmd <- sprintf("tabix %s/%s.all.tsv.gz %s:%d-%d", ftp.uri.base,
                          private$studyIDs[1],
                          private$chrom,
                          private$start,
                          private$end)
           text  <- system(cmd, intern=TRUE)
           tbl <- data.frame()
           if(nchar(text) > 0){
              tmp.filename <- tempfile()
              writeLines(text, con=tmp.filename)
              tbl <- read.table(tmp.filename, sep="\t", as.is=TRUE)
              }
           return(invisible(tbl))
           } # fetch
       ) # public

    ) # class
#--------------------------------------------------------------------------------
