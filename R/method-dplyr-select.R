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
    .data <- internal_filter(x=.data, dots=dots, ...)
    return(.data)
}

##' @method select otu_table
##' @export
select.otu_table <- function(.data, ...){
    dots <- quos(...)
    taxa_are_rows <- .data@taxa_are_rows    
    .data <- internal_select(x=.data, dots=dots, ...)
    .data <- otu_table(object=.data, taxa_are_rows=taxa_are_rows)
    return(.data)
}

##' @method select sample_data
##' @export
select.sample_data <- function(.data, ...){
    dots <- quos(...)
    .data <- internal_select(x=.data, dots=dots, ...)
    .data <- sample_data(.data)
    return(.data)
}

##' @method select taxonomyTable 
##' @export
select.taxonomyTable <- function(.data, ...){
    dots <- quos(...)
    .data <- internal_select(x=.data, dots=dots, ...)
    .data <- tax_table(as.matrix(.data))
    return(.data)
}

internal_select <- function(x, dots, ...){
    x <- as.data.frame(x)
    x %<>% select(!!!dots, ...)
    return(x)
}
