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
    paste("Shuangbin Xu, Li Zhan, Wenli Tang, Qianwen Wang, Zehan Dai, Lang Zhou, Tingze Feng, Meijun Chen, Tianzhi Wu, Erqiang Hu, Guangchuang Yu.",
          "MicrobiotaProcess: A comprehensive R package for deep mining microbiome.",
          "The Innovation. 2023, 4(2):100388. doi: 10.1016/j.xinn.2023.100388\n\n",
          "Export the citation to BibTex by citation('MicrobiotaProcess')\n\n"
    )    
    
}

suppressmsg <- function(pkgname){
    paste0("This message can be suppressed by:\n",
        "suppressPackageStartupMessages(library(", pkgname ,"))"
    )
}
