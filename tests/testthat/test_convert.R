context("dataframe to treedata")
test_that("dataframe to treedata",{
    dir <- system.file("extdata", package="MicrobiotaProcess")
    taxada <- readRDS(file.path(dir, "taxa_tab.rds"))
    rownames(taxada) <- paste0("OTU_", seq_len(nrow(taxada)))
    treeda <- convert_to_treedata(taxada, include.rownames=TRUE)
    expect_true(is(treeda, "treedata"))
})
