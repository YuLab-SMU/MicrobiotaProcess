<!-- README.md is generated from README.Rmd. Please edit that file -->

# MicrobiotaProcess: an R package for analysis, visualization and biomarker discovery of microbiome

MicrobiotaProcess is an R package for analysis, visualization and
biomarker discovery of microbial datasets. It introduces a tidy
microbiome data structure paradigm and analysis grammar. It supports
calculating alpha index and provides functions to visualize rarefaction
curves. Moreover, it also supports visualizing the abundance of taxonomy
of samples. And It also provides functions to perform the PCA, PCoA and
hierarchical cluster analysis. In addition, MicrobiotaProcess also
provides a method for the biomarker discovery of metagenome or other
datasets.

## :writing\_hand: Authors

[Shuangbin Xu](https://github.com/xiangpin) and [Guangchuang
Yu](https://guangchuangyu.github.io)

School of Basic Medical Sciences, Southern Medical University

## :arrow\_double\_down: Installation

Get the released version from
[Bioconductor](https://bioconductor.org/packages/release/bioc/html/MicrobiotaProcess.html):

``` r
## try http:// if https:// URLs are not supported ## the url of mirror
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
## BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("MicrobiotaProcess")
```

the development version from github:

``` r
if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
remotes::install_github("YuLab-SMU/MicrobiotaProcess")
```

## :sparkling\_heart: Contributing

We welcome any contributions\! By participating in this project you
agree to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).
