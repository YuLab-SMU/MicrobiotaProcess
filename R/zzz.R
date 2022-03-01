#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    options(layout.radial.linetype="curved")
    pkgVersion <- packageDescription(pkgname, fields="Version")
    msg <- paste0(pkgname, " v", pkgVersion, "  ",
                  "For help: https://github.com/YuLab-SMU/MicrobiotaProcess/issues", "\n\n")

    citation <- paste0("If you use ", pkgname,
                       " in published research, please cite the paper:\n\n",
                       MP_citations())
    
    packageStartupMessage(paste0(strwrap(pillar::style_subtle(paste0(msg, citation, suppressmsg(pkgname)))), collapse="\n"))
}

MP_citations <- function(){
    paste(
        "S Xu, L Zhan, W Tang, Z Dai, L Zhou, T Feng, M Chen, S Liu, X Fu, T Wu, E Hu, G Yu.",
        "MicrobiotaProcess: A comprehensive R package for managing and analyzing microbiome and other ecological data within the tidy framework.",
        "04 February 2022, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-1284357/v1]\n\n"
        )
}

suppressmsg <- function(pkgname){
    paste0("This message can be suppressed by:\n",
        "suppressPackageStartupMessages(library(", pkgname ,"))"
    )
}
