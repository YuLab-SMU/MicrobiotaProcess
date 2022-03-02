# MicrobiotaProcess 1.7.8.990

+ add `mp_select_as_tip` and fix the bug of `mp_diff_analysis` with specific `tip.level` (not `OTU`) argument. (2022-03-02, Mon)

# MicrobiotaProcess 1.7.8

+ supporting multiple group names and supporting numeric type for `.group` of `mp_plot_alpha`. (2022-02-15, Tue)
+ supporting multiple group names for `.group` of `mp_plot_abundance` when `plot.group=TRUE`. (2022-02-14, Mon)
+ fix the `width` of `mp_plot_abundance` with `geom="flowbar". (2022-02-01, Tue)`
 - related [issue](https://github.com/YuLab-SMU/MicrobiotaProcess/issues/41)

# MicrobiotaProcess 1.7.7

+ fix the bug about the constant variables within groups in `lda` of `MASS` (2022-01-27, Thu) 

# MicrobiotaProcess 1.7.6

+ fix the [issue](https://github.com/YuLab-SMU/MicrobiotaProcess/issues/40), that `kingdom` level of taxonomy 
  information contains `k__` or `K__`, which is unknown annotation in `kingdom`. (2021-01-14, Fri)
+ add `mp_extract_taxatree` and `mp_extract_otutree` (alias of `mp_extract_tree`). (2022-01-14, Fri)

# MicrobiotaProcess 1.7.5

+ add the message for the not integers in `mp_cal_alpha`. (2021-12-31, Fri) 
+ remove the features which variance of their abundance is zero before identify different taxa. (2021-12-28, Tue)
+ add `bar` option in `mp_plot_abundance`, default is `flowbar`, the other options are `bar` and `heatmap`
+ use corrected relative eigenvalues when the eigenvalues has negative values. (2021-12-27, Mon)
+ add new distmethod from `hopach`. (2021-12-27, Mon)
+ update `tax_table` without required `phyloseq`. (2021-12-20, Mon)
+ update `mp_diff_analysis` to support the factor type group (`.group` specified). (2021-12-20, Mon)

# MicrobiotaProcess 1.7.4

+ update `taxatree<-` and `otutree<-` which will extract the intersection between the 
  tip labels of input treedata and the rownames of `MPSE`. (2021-12-14, Tue)
+ add `taxonomy<-` for `MPSE` to assign the taxonomy information, which will be
  converted to `taxatree` automatically. (2021-12-14, Tue) 
+ add `tax_table` for `MPSE` and return `taxonomyTable` defined in `phyloseq`. (2021-12-14, Tue)
+ update `mp_import_metaphlan` to better parse the output of `MetaPhlAn2`. (2021-11-30, Tue)

# MicrobiotaProcess 1.7.3

+ update 'mp_plot_abundance' (2021-11-24, Wed)
  - support `heatmap` by setting `geom`.
  - `.group` supports multiple characters
    and `.sec.group` will be removed in 
    the next version.
+ update `mp_plot_diff_res` (2021-11-23, Tue)
  - support `otutree` and `taxatree` class by setting `tree.type`.
  - support multiple layout types of tree by setting `layout`.
  - support adjusting the gap between panel and width of panel by setting
    `offset.abun`, `pwidth.abun`, `offset.effsize`, `pwidth.effsize`
  - support whether display the relative abundance of `group` 
    instead of `sample` by setting `group.abun=TRUE` or 
    sample number > 50
+ add `mp_plot_diff_res` to visualize the result of mp_diff_analysis. (2021-11-22, Mon)

# MicrobiotaProcess 1.7.2

+ speed up the `mp_cal_abundance`, `mp_cal_venn` and `mp_cal_upset` with `dtplyr`. (2021-11-18, Thu)
+ update the guide of x axis of `ggside` in `mp_plot_ord`. (2021-11-15, Mon)
+ update `mp_plot_abundance` to visualize the abundance of taxonomy from high 
  (bottom) to low (top). (2021-11-15, Mon)
+ support multiple annotation rows or cols of `heatmap` of mp_plot_dist with 
  `.group=c(group1, group2)`, and add `set_scale_theme` to adjust the `scale` 
  or `theme` of subplot of `heatmap`. (2021-11-10, Wed)
+ fix the issue when the taxonomy info is removed with `select`. (2021-11-09, Tue)
+ update `print` for `MPSE` class. (2021-11-09, Tue)
+ update `otutree<-` for support `phylo` class. (2021-11-09, Tue)
+ speed up the integration of `mp_cal_dist` result with `action="add"`. (2021-11-09, Tue)
+ update `as.MPSE` for `biom` class to support parsing the metadata of sample. (2021-11-09, Tue)

# MicrobiotaProcess 1.7.1

+ fix the issue when using `filter` only return a `assays` contained one feature (nrow=1). (2021-11-05, Fri)
+ fix the error of `rownames<-` when `rownames` of `MPSE` is NULL. (2021-11-04, Thu)
+ update 'message' or 'stop error message' when the 'Abundance' cannot be rarefied in some 
  functions, such as `mp_cal_alpha`, `mp_cal_venn`, `mp_cal_upset`, `mp_cal_abundance`
  and `mp_cal_NRT_NTI`. (2021-10-29, Fri)
+ introduce `.sec.group` argument to specify the second group name in `mp_plot_abundance`,
  if it is provided, the nested facet will be displayed. (2021-11-02, Tue)

# MicrobiotaProcess 1.6.0

+ Bioconductor 3.14 release. (2021-10-27, Wed)

# MicrobiotaProcess 1.5.9

+ add `include.rownames` to control whether consider the `OTU` as taxonomy feature table in 
  `diff_analysis` and `get_alltaxadf` or tip labels in `as.treedata`. (2021-10-19, Tue)
+ fix rename bug, rename the taxonomy names can work now. (2021-10-12, Tue)
+ introduce `trimSample` in `mp_rrarefy` to check whether to remove the samples that
  do not have enough abundance. (2021-10-11, Mon)
+ update `MPSE` to allow `assays` supporting `data.frame` or `DFrame` class. (2021-10-08, Fri)
+ update `mp_plot_ord` to suppress the message of the third depend package. (2021-10-08, Fri)

# MicrobiotaProcess 1.5.8

+ fix the bug of `AsIs` list class in `unnest` for the `tidyr` (>= 1.1.4). (2021-10-01, Fri)
+ add `mp_aggregate` function. (2021-09-26, Sun)

# MicrobiotaProcess 1.5.7

+ fix bug of `mp_plot_upset`. (2021-09-10, Fri)
+ update the `mp_plot_ord`. (2021-09-13, Mon)

# MicrobiotaProcess 1.5.6

+ convert the type of first element of assays to `matrix` to compatible with `DESeqDataSet` 
  of `DESeq2`, `test_differential_abundance` of `tidybulk`. (2021-09-09, Thu)
+ update `show` and `print` for format output of `MPSE` class. (2021-09-08, Wed)
+ update `mp_cal_abundance` use new `tidytree`. (2021-09-07, Tue)
+ introduce `include.lowest` parameter in `mp_filter_taxa`. (2021-09-07, Tue)

# MicrobiotaProcess 1.5.5

+ update `mp_plot_ord` to display the `bioplot` for result of `cca`, `rda` and `envfit`. (2021-09-06, Mon)
+ update the vignettes of `MicrobiotaProcess`. (2021-09-04, Sat)
+ return updated `MPSE` object after the `mp_diff_analysis` is done with `action="add"`. (2021-08-31, Fri)
  - then the `taxtree` and `otutree` with the result of different analysis can be extracted with `mp_extract_tree`.
+ fix issue `print` for one line of `MPSE` and update `mp_plot_ord` to display the side `boxplot`. (2021-08-31, Tue)
+ add `mp_plot_ord` for `MPSE` or `tbl_mpse` object after one of `mp_cal_pca`, `mp_cal_pcoa`, `mp_cal_rda`, 
  `mp_cal_nmds`, `mp_cal_rda`, `mp_cal_cca`, `mp_cal_dca` or `mp_envfit` has been run with `action='add'`. (2021-08-30, Mon)
+ add `mp_plot_dist` for `MPSE` or `tbl_mpse` object after `mp_cal_dist` is performed with `action="add"`. (2021-08-28, Sat)
+ add `mp_plot_abundance`, `mp_plot_alpha`, `mp_plot_rarecurve`, `mp_plot_venn`, `mp_plot_upset` for `MPSE` after 
  the corresponding `mp_cal_abundance`, `mp_cal_alpha`, `mp_cal_rarecurve`, `mp_cal_venn`, `mp_cal_upset`
  are performed with `action="add"`. (2021-08-27, Fri)
+ fix the issue when the `rowname` or `colnames` of `SummarizedExperiment` is NULL for `as.MPSE`. (2021-08-26, Thu)

# MicrobiotaProcess 1.5.4

+ fix the `rownames` of `assays` and `colnames` of `colData` to identical for `SummarizedExperiment(1.23.3)`. (2021-08-26, Thu)
+ add `mp_extract_refseq` for `MPSE` object. (2021-08-25, Wed)
+ update `as.MPSE` for `SummarizedExperiment` object. (2021-08-24, Tue)
+ add `mp_filter_taxa` to drop the taxa that low abundance and low occurrences. (2021-08-24, Tue)
+ add `colData<-` and `left_join` for `MPSE`. (2021-08-23, Mon)
+ fix `mutate` for `MPSE` object. 
+ don't import the `parse_taxonomy_greengenes` and `parse_taxonomy_qiime` from `phyloseq`. (2021-08-17, Tue)
+ add `as.MPSE` for `TreeSummarizedExperiment` class. (2021-08-17, Tue)
+ add `mp_import_metaphlan` to parsing the output of `MetaPhlAn`. (2021-08-12, Thu)
  - add `treefile` argument to import the tree of `MetaPhlAn3` (`mpa_v30_CHOCOPhlAn_201901_species_tree.nwk`) (2021-08-13, Fri)
+ update the `print` of `MPSE` object via `pillar` package. (2021-08-06, Fri)
+ update `mp_extract_dist` by introducing `.group` argument to return a `tbl_df` for 
  visualization. (2021-08-04, Wed)
+ add `taxatree`, `taxatree<-`, `otutree`, `otutree<-`, `refseq`, `refseq<-` for `MPSE`. (2021-08-04, Wed)
+ add `mp_extract_rarecurve` to extract the result of `mp_cal_rarecurve` from 
  `MPSE` or `tbl_mpse` object. (2021-08-04, Wed)
+ add `mp_stat_taxa` to count the number and total number taxa for each sample at 
  different taxonomy levels (Kingdom, Phylum, Class, Order, Family, Genus, Species, OTU). (2021-08-03, Tue)

# MicrobiotaProcess 1.5.3

+ rename `mp_extract_abundance` to `mp_extract_assays` from `MPSE` or `tbl_mpse`. (2021-07-31, Sat)
+ update the method to save the result of `mp_cal_clust` by introducing `action` argument. (2021-07-29, Thu).
+ update `as.phyloseq` for `MPSE` or `tbl_mpse` object. (2021-07-28, Wed)
+ add `mp_diff_analysis` for `MPSE` or `tbl_mpse` object. (2021-07-27, Tue)
+ add `dr_extract` for the visualization of the result of ordination. (2021-07-26, Mon)
+ comment out the function for `phyloseq` and add rd of the function for `MPSE` or `tbl_mpse`. (2021-07-24, Sat)
+ update the function to parsing the result of `rda`, `cca`, `envfit`. (2021-07-23, Fri)
+ add `tidydr` to convert the result of `reduce dimension` to `tbl_df`
  - such `pca`, `pcoa`, `nmds`, `rda`, `cca`. (2021-07-22, Thu)
+ optimize the `print` for `MPSE`. (2021-07-22, Thu)

# MicrobiotaProcess 1.5.2

+ add `mp_mantel` and `mp_mrpp` for `MPSE` or `tbl_mpse` object. (2021-07-19, Mon)
+ add `mp_envfit` and update `mp_cal_dist` to support the distance calculation with continuous environment 
  factors and rename `mp_cal_adonis` to `mp_adonis`, `mp_cal_anosim` to `mp_anosim`. (2021-07-17, Sat)
+ add `mp_cal_rda`, `mp_cal_cca`, `mp_cal_adonis` and `mp_cal_anosim` for `MPSE` or `tbl_mpse` object. (2021-07-16, Fri)
+ add `mp_cal_dca`, `mp_cal_nmds` and `mp_extract_internal_attr`. (2021-07-15, Thu)
+ add `mp_cal_pca`, `mp_cal_pcoa` and `mp_extract_abundance`. (2021-07-14, Wed)
+ add `mp_cal_clust` to perform the hierarchical cluster analysis of samples and `mp_extract_dist` to 
  extract the `dist` object from `MPSE` object or `tbl_mpse` object. (2021-07-13, Thu)
+ add `mp_cal_dist` to calculate the distance between samples with `MPSE` or `tbl_mpse` object. (2021-07-12, Mon)
+ add `mp_extract_sample`, `mp_extract_taxonomy`, `mp_extract_feature` to extract the `sample`, `taxonomy`
  and `feature` (`OTU`) information and return `tbl_df` or `data.frame`. (2021-07-09, Fri)
+ add `mp_cal_venn` to build the input for `venn plot` (2021-07-09, Fri)
+ `mp_cal_rarecurve` add `action` argument to control whether the 
  result will be added to `MPSE` and `tbl_mpse` or return directly. (2021-07-08, Thu) 
+ add `mp_cal_upset` to get the input of `ggupset`. (2021-07-08, Thu)
+ add `mp_extract_tree` to extract the `otutree` or `taxatree` from `MPSE` or `tbl_mpse` object. (2021-07-07, Wed)
+ add `pull` and `slice` to support the `MPSE` object. (2021-07-06, Tue)
+ add `mp_cal_rarecurve` to calculate the `rarecurve` of each sample with `MPSE` or `tbl_mpse`. (2021-07-06, Tue)
+ add `mp_cal_abundance` to calculate the relative abundance of each taxonomy class with `MPSE` or `tbl_mpse`. (2021-07-05, Mon)
+ add `mp_decostand` provided several standardization methods for `MPSE`, `tbl_mpse` and `grouped_df_mpse`. (2021-07-04, Sun)
+ add `mp_import_qiime` to parse the output of `qiime` old version. (2021-07-03, Sat)
+ add `taxatree` slot to `MPSE`. (2021-06-30, Wed)
+ add `mp_cal_alpha` function for `MPSE` or `mpse` object. (2021-07-01, Thu)
+ add `rownames<-` to support renaming the names of feature. (2021-07-01, Thu)
+ add `mp_import_qiime2` and `mp_import_dada2` to 
   parse the output of `dada2` or `qiime2` and return `MPSE` object. (2021-07-02, Fri)

+ update `print` information for `MPSE`, `tbl_mpse` and `grouped_df_mpse`. (2021-06-29, Tue)
+ add `[` to the accessors of `MPSE`. (2021-06-29, Tue)

+ use `MPSE` object. (2021-06-28, Mon)
 - add `as.MPSE` to convert `phyloseq` or `tbl_mpse` to `MPSE` class.
 - Formatted output.

+ tidy framework for `MPSE` object.
  - `as_tibble` to convert `MPSE` and `phyloseq` to `tbl_mpse`. (2021-06-28, Mon)
  - `filter` to subset a data frame from `tbl_mpse`. (2021-06-28, Mon)
  - `group_by` to do some data operations on groups for `tbl_mpse`. (2021-06-28, Mon)
  - `arrange` to order the rows of a data frame for `tbl_mpse`. (2021-06-28, Mon)
  - `mutate` to adds new variables and preserves existing ones for `tbl_mpse`. (2021-06-28, Mon)
  - `select` to select variables in `tbl_mpse`. (2021-06-28, Mon)
  - `distinct` to select only unique/distinct rows in `tbl_mpse`. (2021-06-28, Mon)
  - `rename` to rename the variable names in `tbl_mpse`. (2021-06-28, Mon)
  - `nest` to create a list-column of `tbl_mpse`, it will convert `tbl_mpse` to `tbl_mpse_nest`. (2021-06-28, Mon)
  - `unnest` to convert the `tbl_mpse_nest` to `tbl_mpse`. (2021-06-28, Mon)
  - `as.treedata` to convert `tbl_mpse` to `treedata`, then we can explore 
    the data with `treedata`.
    - add `tiplevel` argument to control whether use `OTU` as tip label,
      default is `OTU`. (2021-06-28, Mon)
  - `left_join` to mutate joins based the left `tbl_mpse` structure. (2021-06-28, Mon)
+ changed `clustplotClass` to `treedata`. (2021-06-28, Tue)
+ add `mp_rrarefy` method to rarefy species richness. (2021-06-29, Tue)
  - it supports `MPSE`, `tbl_mpse`, `grouped_df_mpse` object via wrapping `vegan::rrarefy`.
+ ~~update `as.MPSE` and `as.treedata` for `grouped_df_mpse` object. (2021-06-29, Tue)~~
  ~~- This feature is useful to explore the microbiome data in taxa tree.~~ 
  This feature has been replaced by the `taxatree` slot
  
# MicrobiotaProcess 1.5.1

+ add `ellipse_linewd` and `ellipse_lty` in `ggordpoint` to control 
  the width and line type of ellipse line. (2021-05-24, Mon)
+ fixed the regular expression match for the internal function 
  to print the results of `diff_analysis`. (2021-06-06, Sun)
+ add `filter` function to filter the result of `diff_analysis`. (2021-06-07, Mon)
+ more accessor function for result of `diff_analysis`. (2021-06-07, Mon)
  - `head`
  - `tail`
  - `[`
  - `[[]]`
  - `$`
  - `dim`
+ add `get_NRI_NTI` to calculate the `NRI` and `NTI`. (2021-06-08, Tue)

# MicrobiotaProcess 1.4.0

+ new version released. (2021-05-20, Thu)

# MicrobiotaProcess 1.3.11

+ fill `ggclust` bug to map `color` and `shape`. (2021-05-12, Wed)

# MicrobiotaProcess 1.3.9

- more layouts for `ggdiffclade`. (2021-04-16, Fri)
+ remove unclassified, ambiguous taxonomy names. (2021-03-30, Tue)
+ <https://github.com/YuLab-SMU/MicrobiotaProcess/issues/23>

# MicrobiotaProcess 1.3.8 

+ rename files of code. (2021-03-23, Tue)
+ add aliases for `ggbartaxa` and `ggdiffbartaxa`. (2021-03-23, Tue)

# MicrobiotaProcess 1.3.7

+ optimized import for installation. (2021-03-15, Mon) 
+ add `sampleLevels` in `ggbartax` to adjust the order of axis label. (2021-03-12, Fri)
+ update `import_qiime2` to avoid error when only feature table is provided. (2021-02-26, Fri)

# MicrobiotaProcess 1.3.6

+ convert `svg` dev to `pdf` dev. (2021-02-04, Thu)

# MicrobiotaProcess 1.3.5

+ fix an error for example of `ggrarecurve`. (2021-01-07, Thu)
  `factorNames="Group"` to `factorNames="group"`

# MicrobiotaProcess 1.3.4

+ fix a bug for numeric sample name. (2020-11-26, Thu)
+ `geom_tiplab` also support circular layout, so remove `geom_tiplab2`. (2020-11-26, Thu)

# MicrobiotaProcess 1.3.3

+ add `as.treedata` for `taxonomyTable` class. (2020-11-23, Mon)

# MicrobiotaProcess 1.3.2

+ `ggrarecurve` can be set color with variable of group for each samples. (2020-11-11, Tue)
  - using `shadow=FALSE` and providing `factorNames`
  - <https://github.com/YuLab-SMU/MicrobiotaProcess/issues/21>
+ add `get_rarecurve` to avoid repeated calculation when displaying rare curve. (2020-11-17, Tue)
  + `rareres <- get_rarecurve(obj, chunks=400)`
    `p <- ggrarecurve(rareres)`

# MicrobiotaProcess 1.3.1

+ `ggordpoint` add `showsample` to show the labels of sample. (2020-10-29, Thu)
+ the point of `ggordpoint` use the points of [`ggstar`](https://github.com/xiangpin/ggstar). (2020-10-30, Fri)
+ to obtain the dynamic arguments of `diff_analysis`, the `call` was changed to `someparams`.
  `someparams` contained the arguments used in other functions. (2020-11-09, Mon) 
  - <https://github.com/YuLab-SMU/MicrobiotaProcess/issues/20>

# MicrobiotaProcess 1.2.0

+ Bioconductor 3.12 release (2020-10-28, Wed)

# MicrobiotaProcess 1.1.13

+ removed `retrieve_seq` and `mapply_retrieve_seq` function, since these need internet. 
  Which might cause time out when check. (2020-10-16, Fri)

# MicrobiotaProcess 1.1.12

+ modified a bug in diff_analysis.phyloseq: change `tax_table(ps)` to `ps@tax_table` to 
  avoid generate error when tax_table is NULL. (2020-10-15, Thu)

# MicrobiotaProcess 1.1.11

+ update the examples of `drop_taxa`. (2020-10-14, Wed)

# MicrobiotaProcess 1.1.10

+ update `ggdiffclade` to support data.frame input when `reduce` is `TRUE`. (2020-08-28, Fri)

# MicrobiotaProcess 1.1.9

+ update `ggordpoint` to fit the usage when user want to set mapping by manually. (2020-08-25, Tue)

# MicrobiotaProcess 1.1.8

+ `get_taxadf`, `get_alltaxadf` and `diff_analysis` has supported function datasets or other type datasets. (2020-08-17, Mon)

# MicrobiotaProcess 1.1.7

+ bugfix: `cladetext` argument has been omitted in `ggdiffclade`, now it has been fixed. (2020-08-14, Fri)
+ deprecated argument: the `size` argument controlled the width of line of tree has been deprecated.
  The `linewd` replace it (2020-08-14, Fri).

# MicrobiotaProcess 1.1.6

+ `removeUnkown` argument has been replaced with `removeUnknown` in `ggdiffbox`,
  `ggeffectsize`, `ggdifftaxbar` and `ggdiffclade`. (2020-08-12, Wed)
+ `class` argument has been replaced with `classgroup` in `diff_analysis`. (2020-08-12, Wed)
+ add `inward_circular` layout in `ggdiffclade`. (2020-08-12, Wed)

# MicrobiotaProcess 1.1.5

+ `ggdifftaxbar` supports `png`, `tiff` format. (2020-08-10, Mon)
+ add stop information to state the class argument in `diff_analysis`. (2020-08-10, Mon)

# MicrobiotaProcess 1.1.4

+ add `tax_table` information to result of `get_taxadf`. (2020-08-07, Fri)

# MicrobiotaProcess 1.1.3

+ change according to dplyr (v=1.0.0) (2020-08-05, Wed)
  + remove `rename_` and `group_by_`
+ modified the angle to 90 in `ggdiffclade` when layout is `slanted` or `rectangular` (2020-08-05, Wed)

# MicrobiotaProcess 1.1.2

+ fix a bug. When the first rank taxa level (Kingdom) is NA.


# Changes in version 0.99.1 (2020-03-03)

+ First release of Bioconductor

# Changes in version 0.99.0 (2019-08-14)

+ Submitted to Bioconductor
