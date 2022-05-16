#' \code{clamqc} package
#'
#' @docType package
#' @name clamqc
#'
#' @description Calico Lipidomics and Metabolomics QC-associated R functions features methods for working with .mzrollDB mass spec datasets.
#'
#' @importFrom dplyr %>%
#' @importFrom stats rt
#' @importFrom utils data
#' @import ggplot2
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
utils::globalVariables(c(
  "annotatedRt",
  "basepeak",
  "category",
  "compound_max_ic",
  "data",
  "desc",
  "eic",
  "filename",
  "fraction_max",
  "groupId",
  "ic",
  "ic_adjusted",
  "join_col",
  "mean_peakAreaTop",
  "mean_peakMz",
  "mean_peakQuality",
  "mean_peakRt",
  "methodId",
  "name",
  "observedRt",
  "retentionTime",
  "rt",
  "rt_type",
  "sample_standards",
  "sample_type",
  "samplename",
  "samplename_factor",
  "standard",
  "standard_data",
  "standardType",
  "y_offset"
))