# MicrobiotaProcess 1.3.6.991

+ add `sampleLevels` in `ggbartax` to adjust the order of axis label. (2021-03-12, Fri)

# MicrobiotaProcess 1.3.6.990

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
