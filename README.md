<!-- README.md is generated from README.Rmd. Please edit that file -->

# MicrobiotaProcess: A comprehensive R package for managing and analyzing microbiome and other ecological data within the tidy framework

[![](https://img.shields.io/badge/release%20version-1.8.1-green.svg)](https://www.bioconductor.org/packages/MicrobiotaProcess)
[![](https://img.shields.io/badge/devel%20version-1.9.2.990-green.svg)](https://github.com/YuLab-SMU/MicrobiotaProcess)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/MicrobiotaProcess.svg)](https://www.bioconductor.org/packages/devel/bioc/html/MicrobiotaProcess.html#since)
[![platform](http://www.bioconductor.org/shields/availability/devel/MicrobiotaProcess.svg)](https://www.bioconductor.org/packages/devel/bioc/html/MicrobiotaProcess.html#archives)

This package has many unique features:

-   **MicrobiotaProcess** introduces **MPSE** class, which make this
    package more interoperable with the existing computing ecosystem.
-   **MicrobiotaProcess** can bridge several common tools of microbiome
    analysis with several parsing functions.
-   **MicrobiotaProcess** introduces a tidy microbiome data structure
    paradigm and analysis grammar via formatted output avoiding memory
    consumption.
-   **MicrobiotaProcess** provides a wide variety of microbiome analysis
    procedures under the unified and common framework (tidy-like
    framework). This will make the related analysis become easier.

## Anatomy of a **MPSE**

<img src="./inst/figures/mpse.png" title="The structure of the MPSE class." alt="The structure of the MPSE class." style="display: block; margin: auto;" />

## Overview of the design of **MicrobiotaProcess** package

<img src="./inst/figures/mp-design.png" title="The Overview of the design of MicrobiotaProcess package" alt="The Overview of the design of MicrobiotaProcess package" style="display: block; margin: auto;" />

## :writing_hand: Authors

[Shuangbin Xu](https://github.com/xiangpin) and [Guangchuang
Yu](https://guangchuangyu.github.io)

School of Basic Medical Sciences, Southern Medical University

## :arrow_double_down: Installation

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

# :book: Vignette

For more details, please refer to the [online
vignette](https://bioconductor.org/packages/devel/bioc/vignettes/MicrobiotaProcess/inst/doc//MicrobiotaProcess.html).

## :sparkling_heart: Contributing

We welcome any contributions! By participating in this project you agree
to abide by the terms outlined in the [Contributor Code of
Conduct](CONDUCT.md).
