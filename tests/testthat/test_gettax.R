context("mp_cal_abundance and mp_extract_abundance work")
test_that("mp_cal_abundance and mp_extract_abundance work",{
    data(test_otu_data)
    otuda <- test_otu_data %>% 
             mp_extract_assays(.abundance=Abundance)
    taxada <- test_otu_data %>%
              mp_extract_taxonomy() %>%
              tibble::column_to_rownames(var="OTU")
    phyla1 <- test_otu_data %>% 
        mp_cal_abundance(.abundance=Abundance, force=TRUE) %>% 
        mp_extract_abundance(taxa.class=Phylum) %>% 
        tidyr::unnest(AbundanceBySample) %>% 
        select(label, Sample, RelAbundanceBySample) %>% 
        tidyr::pivot_wider(id_cols="label", names_from=Sample, values_from=RelAbundanceBySample) %>% 
        tibble::column_to_rownames(var="label")
    phyla2 <- get_ratio(otuda, taxada[, 2, drop=FALSE])
    testthat::expect_equal(phyla1, phyla2 * 100)
})  
