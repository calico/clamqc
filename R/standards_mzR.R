#' QC Extract Standards
#'
#' @param samples_df a tibble containing samplename and filename
#' @param query_methodId the analytical method being used
#' @param clamr_config a named list of mass spec parameters with special formatting of instrument tolerances
#' @param mzroll_db_path path to a .mzrollDB SQLite file.
#' @param mass_spec_standards_con A connection to an SQL database which stores
#'   experimental mass spectrometry data.
#'
#' @return Extracted data for each internal and external standard
#'
#' @export
qc_extract_standards <- function(samples_df, clamr_config, mzroll_db_path = NULL, mass_spec_standards_con, query_methodId) {
  debugr::dwatch(msg = "Starting [(standards_mzR.R) clamqc::qc_extract_standards]")

  clamr::require_tolerances(clamr_config, 1L)
  stopifnot(class(query_methodId) == "character", length(query_methodId) == 1)

  possible_qcStandards <- dplyr::tbl(mass_spec_standards_con, dbplyr::sql("SELECT * FROM qcStandards")) %>%
    dplyr::collect() %>%
    dplyr::filter(methodId == query_methodId)

  # external standards

  external_standards <- samples_df %>%
    # match external standards
    dplyr::mutate(standard = stringr::str_extract(samplename, "std[0-9]{3}[A-Z]+")) %>%
    dplyr::filter(!is.na(standard))

  undefined_external_standards <- unique(external_standards$standard[!(external_standards$standard %in% possible_qcStandards$standard)])
  if (length(undefined_external_standards) != 0) {
    warning("internal standard composition is unknown for ", paste(undefined_external_standards, collapse = ", "))
  }

  matched_external_standards <- external_standards %>%
    dplyr::inner_join(possible_qcStandards, by = "standard") %>%
    dplyr::mutate(standard = paste0("External-", standard))

  # internal standards

  possible_internal_standards <- possible_qcStandards %>%
    dplyr::filter(standardType == "Internal") %>%
    dplyr::mutate(standard = "Internal")

  if (nrow(possible_internal_standards) != 0) {
    matched_internal_standards <- samples_df %>%
      # remove external standards
      dplyr::filter(!stringr::str_detect(samplename, "std[0-9]{3}[A-Z]+")) %>%
      dplyr::mutate(methodId = query_methodId) %>%
      # match to method-matched internal standards
      dplyr::left_join(possible_internal_standards, by = "methodId")
  } else {
    matched_internal_standards <- tibble::tibble()
  }

  before <- Sys.time()

  matched_standards <- dplyr::bind_rows(
    matched_external_standards,
    matched_internal_standards
  ) %>%
    tidyr::nest(sample_standards = c(-samplename, -filename)) %>%
    dplyr::mutate(data = furrr::future_map2(filename, sample_standards,
      clamr::tidy_mzR_EICs,
      mzroll_db_path = mzroll_db_path,
      clamr_config = clamr_config
    )) %>%
    dplyr::select(-sample_standards) %>%
    tidyr::unnest(data) %>%
    tidyr::nest(standard_data = c(-standard, -methodId, -standardType))

  after <- Sys.time()
  debugr::dwatch(msg = paste("determined matched_standards in", difftime(after, before, units = "secs"), "[(standards_mzR.R) clamqc::qc_extract_standards]"))

  return(matched_standards)
}

#' QC Standard Summaries
#'
#' @param extracted_standards Output of \link{qc_extract_standards}
#'
#' @return a table of sample ~ standard ion counts
#'
#' @export
qc_standard_summaries <- function(extracted_standards) {
  distinct_samples <- extracted_standards %>%
    tidyr::unnest(standard_data) %>%
    dplyr::distinct(samplename) %>%
    dplyr::arrange(samplename) %>%
    dplyr::mutate(samplename_factor = factor(samplename, levels = samplename))

  standard_peaks <- extracted_standards %>%
    tidyr::unnest(standard_data) %>%
    tidyr::unnest(eic) %>%
    dplyr::group_by(samplename, name) %>%
    dplyr::arrange(desc(ic)) %>%
    dplyr::slice(1) %>%
    dplyr::left_join(distinct_samples, by = "samplename") %>%
    dplyr::ungroup() %>%
    dplyr::arrange(samplename_factor) %>%
    dplyr::ungroup() %>%
    dplyr::select(samplename, name, ic) %>%
    tidyr::spread(name, ic, fill = 0)

  return(standard_peaks)
}
