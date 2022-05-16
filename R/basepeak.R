#' Basepeak
#'
#' Find the largest ion count for each MS1 scan
#'
#' @param ms_file_path Path to an mzXML, mzML, netCDF or mzData file
#'
#' @return a tibble with one row per MS1 scan and scan, rt, ion count columns
#'
#' @export
mzR_basepeak <- function(ms_file_path) {
  if (!file.exists(ms_file_path)) {
    stop("Sample not found at: ", ms_file_path)
  }

  an_ms_file <- clamr::tidy_mzR_from_msfile(ms_file_path, msLevels = 1L)

  basepeak <- an_ms_file$scan_data %>%
    dplyr::group_by(scan) %>%
    dplyr::arrange(desc(ic)) %>%
    dplyr::slice(1) %>%
    dplyr::left_join(an_ms_file$header %>%
      dplyr::select(scan, rt = retentionTime) %>%
      dplyr::mutate(rt = rt / 60),
    by = "scan"
    ) %>%
    dplyr::ungroup()

  basepeak
}

#' Plot Basepeak of Samples Compared to Reference
#'
#' @inheritParams qc_extract_standards
#' @param reference_sample path to a reference sample
#' @param n_plot_samples number of samples to plot
#'
#' @export
basepeak_comparison <- function(samples_df, reference_sample, n_plot_samples = 12) {
  stopifnot(all(class(reference_sample) == "character"), length(reference_sample) == 1)
  stopifnot(class(n_plot_samples) %in% c("numeric", "integer"), length(n_plot_samples) == 1, n_plot_samples > 0)

  if (n_plot_samples > nrow(samples_df)) {
    n_plot_samples <- nrow(samples_df)
  }

  if (!file.exists(reference_sample)) {
    stop("Reference sample not found at: ", reference_sample)
  }

  missing_samples <- samples_df$samplename[!file.exists(samples_df$filename)]
  if (length(missing_samples) != 0) {
    stop(length(missing_samples), " samples could not be located: ", paste(missing_samples, collapse = " & "))
  }


  # extract basepeak
  sample_basepeak_df <- samples_df %>%
    dplyr::sample_n(n_plot_samples) %>%
    dplyr::mutate(basepeak = purrr::map(filename, mzR_basepeak)) %>%
    dplyr::select(-filename) %>%
    tidyr::unnest(basepeak)

  reference_basepeak <- mzR_basepeak(reference_sample)

  basepeak_comparison_plot <- plot_basepeak_comparison(sample_basepeak_df, reference_basepeak)

  basepeak_comparison_plot
}

plot_basepeak_comparison <- function(sample_basepeak_df, reference_basepeak) {

  # TO DO: replace this block with dplyr::crossing
  reference_expansion <- sample_basepeak_df %>%
    dplyr::distinct(samplename) %>%
    dplyr::mutate(join_col = 1) %>%
    dplyr::left_join(
      reference_basepeak %>%
        dplyr::mutate(join_col = 1),
      by = "join_col"
    ) %>%
    dplyr::select(-join_col)

  sample_basepeaks <- dplyr::bind_rows(
    reference_expansion %>%
      dplyr::mutate(
        category = "reference",
        ic = -1 * ic
      ),
    sample_basepeak_df %>%
      dplyr::mutate(category = "sample")
  )

  ggplot(sample_basepeaks, aes(x = rt, y = ic, color = category)) +
    geom_path() +
    facet_wrap(~samplename, scale = "free_y") +
    theme_bw()
}
