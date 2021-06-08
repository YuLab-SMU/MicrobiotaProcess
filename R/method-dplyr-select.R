##' @method select diffAnalysisClass 
##' @importFrom dplyr select
##' @export
select.diffAnalysisClass <- function(.data, ...) {
    dots <- quos(...)
    .data@result %<>% select(!!!dots, ...)
    return(.data)
}

##' @method select alphasample
##' @importFrom dplyr select
##' @export
select.alphasample <- function(.data, ...) {
    dots <- quos(...)
    .data <- as.data.frame(.data)
    .data %<>% select(!!!dots, ...)
    return(.data)
}
