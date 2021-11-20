library(RUnit)
library(ADvariantExplorer)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    id <- "abc"
    rt <- ADvariantExplorer$new(id)
    checkTrue(all(c("R6", "ADvariantExplorer") %in% class(rt)))
    checkEquals(rt$getID(), id)

} # test_ctor
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
