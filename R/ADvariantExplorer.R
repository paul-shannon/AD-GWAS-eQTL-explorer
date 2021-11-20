#' @title ADvariantExplorer
#' @description A template for building documented, tested R6 classes
#' @name ADvariantExplorer

#' @field id identifier for a class object
#'
#' @examples
#'   rt <- R6Template$new(id="abc")
#' @export

ADvariantExplorer = R6Class("ADvariantExplorer",

    #--------------------------------------------------------------------------------
    private = list(id=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param id character, an indentifier for this object
         #' @return a new instance of ADvariantExplorer
        initialize = function(id){
            private$id <- id
            },
        #------------------------------------------------------------
        #' @description accessor for the object's id field
        #' @return the current value of the id member
        getID = function(){
            private$id
            }
       ) # public

    ) # class
#--------------------------------------------------------------------------------
