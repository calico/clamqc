#' Standard QC for Peakgroups
#'
#' @param mzroll_db_con a connection to a mzroll database
#' @param clamr_config a named list of mass spec parameters with special formatting of instrument tolerances
#' @param mass_spec_standards_con A connection to an SQL database which stores
#'   experimental mass spectrometry data.
#' @inheritParams qc_extract_standards
#'
#' @return the most plausible peakgroup label for each standard
#'
#' @export
standard_qc_peakgroups <- function(mzroll_db_con, clamr_config, mass_spec_standards_con, query_methodId) {
  possible_qcStandards <- dplyr::tbl(mass_spec_standards_con, dbplyr::sql("SELECT * FROM qcStandards")) %>%
    dplyr::collect() %>%
    dplyr::filter(methodId == query_methodId)

  if (nrow(possible_qcStandards) == 0) {
    warning("Zero standards present for ", query_methodId)
    return(NULL)
  }

  clamr::require_tolerances(clamr_config, 1L)
  variable_tolerances <- tibble::tibble(
    variable = "mz",
    tolerance = clamr_config$MS1tol$tol,
    relative_or_absolute = clamr_config$MS1tol$absolute_or_relative
  )

  # summarize peakgroups
  peakgroups_aug <- clamr::augment_peakgroups(mzroll_db_con)

  # match all peakgroups to standards based on mass-tolerance
  peak_standard_join <- peakgroups_aug %>%
    dplyr::select(groupId, mz = mean_peakMz, mean_peakRt, mean_peakAreaTop, mean_peakQuality) %>%
    clamr::join_matching_standards(possible_qcStandards, variable_tolerances, threshold = 1)

  # annotate a standard as the largest peaks within mass tolerance
  peak_standard_join %>%
    dplyr::arrange(desc(mean_peakAreaTop)) %>%
    dplyr::group_by(name, standard, standardType) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}
