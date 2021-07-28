# MicrobiotaProcess 1.5.2.990

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
