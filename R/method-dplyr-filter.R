##' @method filter diffAnalysisClass
##' @importFrom rlang quos
##' @export
filter.diffAnalysisClass <- function(.data, ..., .preserve = FALSE) {
    dots <- quos(...)
    .data@result %<>% filter(!!!dots, .preserve = .preserve)
    return(.data)
}

##' @method filter alphasample
##' @export
filter.alphasample <- function(.data, ..., .preserve = FALSE){
    dots <- quos(...)
    .data <- as.data.frame(.data)
    .data %<>% filter(!!!dots, .preserve = .preserve)
    return(.data)
}
