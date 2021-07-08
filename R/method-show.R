#' @title method extensions to show for diffAnalysisClass or alphasample objects.
#' @rdname show-methods
#' @param object object, diffAnalysisClass or alphasample class
#' @importFrom methods show
#' @author Shuangbin Xu
#' @return print info
#' @export
#' @examples
#' \dontrun{
#' data(kostic2012crc)
#' kostic2012crc
#' head(phyloseq::sample_data(kostic2012crc),3)
#' kostic2012crc <- phyloseq::rarefy_even_depth(kostic2012crc,rngseed=1024)
#' table(phyloseq::sample_data(kostic2012crc)$DIAGNOSIS)
#' set.seed(1024)
#' diffres <- diff_analysis(kostic2012crc, classgroup="DIAGNOSIS",
#'                         mlfun="lda", filtermod="fdr",
#'                         firstcomfun = "kruskal.test",
#'                         firstalpha=0.05, strictmod=TRUE, 
#'                         secondcomfun = "wilcox.test",
#'                         subclmin=3, subclwilc=TRUE,
#'                         secondalpha=0.01, lda=3)
#' show(diffres)
#' }
setMethod("show", 
    "diffAnalysisClass",
    function(object){
      originalD <- object@originalD
    	cat(paste0("The original data: ", ncol(originalD),
      			 " features and ", nrow(originalD)," samples"),
      	  fill=TRUE)
      sampleda <- object@sampleda
      cat(paste0("The sample data: ", ncol(sampleda), " variables and ", nrow(sampleda), " samples"),
      	fill=TRUE)
      taxda <- object@taxda
      if(!is.null(taxda)){cat(paste0("The taxda contained ", nrow(taxda), " by ",ncol(taxda), " rank"),
      						fill=TRUE)}
      else{cat("The taxda is NULL",fill=TRUE)}

      
      firstvars <- extract_first_vars(obj=object) 
      firstfun <- extract_args(object, "firstcomfun")
      filtermod <- extract_args(object, "filtermod")
      alphafold <- extract_args(object, "firstalpha")
      cat(paste0("after first test (",firstfun,") number of feature (", filtermod,"<=",alphafold,"):", 
      length(firstvars)),fill=TRUE)
      secondvars <- get_second_true_var(object)
      #secondfun <- extract_args(object, "secondcomfun")
      secondfun <- check_second_fun(obj=object)
      cat(paste0("after second test (",secondfun,") number of significantly discriminative feature:", 
      		   nrow(secondvars)),
          fill=TRUE)
      mlres <- object@result 
      uncertain <- length(grep("__un_", mlres$f))
      mlmethod <- extract_args(object, "mlfun")
      cat(paste0("after ",mlmethod,", Number of discriminative features: ", 
      		   nrow(mlres), " (certain taxonomy classification:", 
      		   nrow(mlres) -uncertain , 
      		   "; uncertain taxonomy classication: ",uncertain,")"), 
      	fill=TRUE)
    }
)


extract_first_vars <- function(obj){
    filtermod <- extract_args(obj, "filtermod")
    firstalpha <- extract_args(obj, "firstalpha")
    kwres <- obj@kwres
    if (filtermod !="pvalue"){
        varsfirst <- kwres[kwres$fdr<=firstalpha& !is.na(kwres$fdr),,drop=FALSE]
    }else{
        varsfirst <- kwres[kwres$pvalue<=firstalpha&!is.na(kwres$pvalue),,drop=FALSE]
    }
    return (as.vector(varsfirst$f))
}

check_second_fun <- function(obj){
    subclass <- extract_args(obj, "subclass")
    classgroup <- extract_args(obj, "classgroup")
    strictmod <- extract_args(obj, "strictmod")
    clmin <- extract_args(obj, "clmin")
    clwilc <- extract_args(obj, "clwilc")
    fcfun <- extract_args(obj, "fcfun")
    secondcomfun <- extract_args(obj, "secondcomfun")
    subclmin <- extract_args(obj, "subclmin")
    subclwilc <- extract_args(obj, "subclwilc")
    if (!is.null(subclass) && strictmod){
        submin <- min(table(obj@sampleda[[subclass]]))
        if (submin >= subclmin && subclwilc){
            secondfun <- paste(secondcomfun, "and", fcfun)
        }else{
            secondfun <- fcfun
        }
    }else{
        groupmin <- min(table(obj@sampleda[[classgroup]]))
        if (groupmin >= clmin && clwilc){
            secondfun <- paste(secondcomfun, "and", fcfun)
        }else{
            secondfun <- fcfun
        }
    }
    return (secondfun)
}

#' @rdname show-methods
#' @exportMethod show
setMethod("show", "alphasample", function(object){
    print.alphasample(object)
})

print.alphasample <- function(object){
    msg <- "'alphasample' S4 object"
    sampleda <- object@sampleda
    if (!is.null(sampleda)){
        smsg <- sampleda
    }else{
        smsg <- "The 'alphasample' does not have sampleda slot."
    }
    cat (msg, fill=TRUE)
    cat ("The alpha diversity index :", fill=TRUE)
    print (head(object@alpha))
    if (is.character(smsg)){
        cat (smsg, fill=TRUE)
    }else{
        cat ("The sample information is present in 'alphasample':", fill=TRUE)
        print (head(sampleda))
    }
}

#' @rdname show-methods
#' @exportMethod show
setMethod("show", "MPSE", function(object){
    if (isTRUE(x = getOption(x = "restore_MPSE_show", default = FALSE))) {
        writeLines(formatted_out("The otutree (treedata object) of the MPSE object is: "))
        if (!is.null(object@otutree)){
            show(object@otutree)
        }else{
            writeLines("NULL")
        }
        writeLines(formatted_out("The taxatree (treedata object) of the MPSE object is: "))
        if (!is.null(object@taxatree)){
            show(object@taxatree)
        }else{
            writeLines(formatted_out("NULL"))
        }
        writeLines(formatted_out("The reference sequence (XStringSet object) of the MPSE object is: "))
        if (!is.null(object@refseq)){
            show(object@refseq)
        }else{
            writeLines("NULL")
        }
        writeLines(formatted_out("The abundance and sample data of the MPSE object are: "))
        f <- getMethod(f="show", signature = "SummarizedExperiment", where = asNamespace(ns = "SummarizedExperiment"))
        f(object)
    }else{
        object %>% print()
    }
})


#' @title print some objects
#' @name print
#' @param x Object to format or print.
#' @param ... Other arguments passed on to individual methods.
#' @param n Number of rows to show. If `NULL`, the default, will print all rows
#'   if less than option `tibble.print_max`. Otherwise, will print
#'   `tibble.print_min` rows.
#' @param width Width of text output to generate. This defaults to `NULL`, which
#'   means use `getOption("tibble.width")` or (if also `NULL`)
#'   `getOption("width")`; the latter displays only the columns that fit on one
#'   screen. You can also set `options(tibble.width = Inf)` to override this
#'   default and always print all columns.
#' @param n_extra Number of extra columns to print abbreviated information for,
#'   if the width is too small for the entire tibble. If `NULL`, the default,
#'   will print information about at most `tibble.max_extra_cols` extra columns.
#' @return print information
NULL

#' @rdname print
#' @method print MPSE
#' @export
print.MPSE <- function(x, ..., n = NULL, width = NULL, n_extra = NULL) {
    if (nrow(x) > 30){
        tmpx <- x[1:min(40, nrow(x)), min(1, ncol(x)), drop=FALSE]
    }else{
        tmpx <- x[,1:min(20, ncol(x)), drop=FALSE]
    }
    formatted_tb <- tmpx %>% 
                    as_tibble() %>% 
                    format(..., n = n, width = width, n_extra = n_extra)
    total_nrows <- dim(x)[1] * dim(x)[2]
    show_nrows <- ifelse(is.null(n), 10, n)
    
    new_head = sprintf(
      "A MPSE-tibble (MPSE object) abstraction: %s",
       total_nrows %>% format(format="f", big.mark=",", digits=2)
    )
    
    left_nrows <- total_nrows - show_nrows
    new_tail <- sprintf("%s more rows", 
       left_nrows %>% format(format="f", big.mark=",", digits=2)
    )
    formatted_mpse = 
      formatted_tb %>%
      {
        x = (.);
        x[1] = gsub("(A tibble: [0-9,]+)", new_head, x[1]);
        x[show_nrows+4] = gsub("([0-9,]+ more rows)", new_tail, x[show_nrows + 4]);
        x
      }
    formatted_mpse %>%
      append(sprintf(
        "\033[90m# OTU=%s | Samples=%s | Assays=%s | Taxanomy=%s\033[39m",
        nrow(x), 
        ncol(x),
        SummarizedExperiment::assays(x) %>% names %>% paste(collapse=", "),
        ifelse(is.null(x@taxatree), "NULL", 
               x@taxatree@data %>% 
                   select("nodeClass") %>% 
                   filter(! .data[["nodeClass"]] %in% c("OTU", "Root")) %>% 
                   unlist(use.names=FALSE) %>% 
                   unique %>% 
                   paste(collapse=", "))
      ), after = 1) %>%
      writeLines()
    invisible(x)
}

#' @method print tbl_mpse
#' @rdname print
#' @export
print.tbl_mpse <- function(x, ..., n = NULL, width = NULL, n_extra = NULL){
    formatted_tb <- x %>% format(..., n = n, width = width, n_extra = n_extra)
    if (valid_names(x, type="tbl_mpse")){
        new_head = "A tbl_mpse (which can be converted to MPSE via as.MPSE) abstraction:"
        formatted_tb_mpse <-
            formatted_tb %>% 
            {
               x = (.);
               x[1] = gsub("(A tibble:)", new_head, x[1]);
               x
            }
        writeLines(formatted_tb_mpse)
    }else{
        writeLines(formatted_tb)
    }
    invisible(x)
}

#' @method print grouped_df_mpse
#' @rdname print
#' @export
print.grouped_df_mpse <- function(x, ..., n = NULL, width = NULL, n_extra = NULL){
    formatted_tb <- x %>% format(..., n = n, width = width, n_extra = n_extra)
    if (valid_names(x, type="grouped_df_mpse")){
        new_head = "A grouped_df_mpse (which can be converted to MPSE via as.MPSE) abstraction:"
        formatted_grouped <- 
            formatted_tb %>%
            {
                x = (.);
                x[1] = gsub("(A tibble:)", new_head, x[1]);
                x
            }
        writeLines(formatted_grouped)
    }else{
        writeLines(formatted_tb)
    }
    invisible(x)
}

#' @method print rarecurve
#' @rdname print
#' @export
print.rarecurve <- function(x, ..., n = NULL, width = NULL, n_extra = NULL){
    formatted_tb <- x$data %>% format(..., n =  n, width = width, n_extra = n_extra)
    new_head = "A rarecurve (which can be visualized via ggrarecurve) abstraction:"

    format_rare <- 
        formatted_tb %>%
        {
            x = (.);
            x[1] = gsub("(A tibble:)", new_head, x[1]);
            x
        }
    writeLines(format_rare)

    invisible(x)
}
