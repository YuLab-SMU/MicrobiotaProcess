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
    .data <- internal_filter(x=.data, ..., .preserve = .preserve)
    return(.data)
}

##' @method filter otu_table
##' @export 
filter.otu_table <- function(.data, ..., .preserve = FALSE){
    taxa_are_rows <- .data@taxa_are_rows
    .data <- internal_filter(x=.data, ..., .preserve = .preserve)
    .data <- otu_table(object=.data, taxa_are_rows=taxa_are_rows)
    return(.data)
}

##' @method filter sample_data
##' @export
filter.sample_data <- function(.data, ..., .preserve = FALSE){
    .data <- internal_filter(x=.data, ..., .preserve = .preserve)
    .data <- sample_data(.data)
    return(.data)
}

##' @method filter taxonomyTable
##' @export
filter.taxonomyTable <- function(.data, ..., .preserve = FALSE){
    .data <- internal_filter(x=.data, ..., .preserve=.preserve)
    .data <- tax_table(as.matrix(.data))
    return(.data)
}

internal_filter <- function(x, ..., .preserve){
    dots <- quos(...)
    x <- as.data.frame(x)
    x %<>% filter(!!!dots, .preserve = .preserve)
    return(x)
}
