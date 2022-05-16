#' QC Standard Plots
#'
#' @param standard_data produced by \code{\link{qc_extract_standards}}
#'
#' @return a tibble containing standard QC grobs
#'
#' @export
qc_standard_plots <- function(standard_data) {
  distinct_samples <- standard_data %>%
    dplyr::distinct(samplename) %>%
    dplyr::arrange(samplename) %>%
    dplyr::mutate(
      samplename_factor = factor(samplename, levels = samplename),
      sample_type = dplyr::case_when(
        stringr::str_detect(tolower(samplename), "blank") ~ "blank",
        stringr::str_detect(tolower(samplename), "std|standard") ~ "standard",
        TRUE ~ "sample"
      )
    )

  standard_peaks <- standard_data %>%
    tidyr::unnest(eic) %>%
    dplyr::group_by(samplename, name) %>%
    dplyr::arrange(desc(ic)) %>%
    dplyr::slice(1) %>%
    dplyr::left_join(distinct_samples, by = "samplename") %>%
    dplyr::ungroup()

  # summary of expected RT
  standard_reference_values <- standard_peaks %>%
    dplyr::distinct(name, observedRt, annotatedRt) %>%
    tidyr::gather("rt_type", "rt", -name) %>%
    dplyr::filter(!is.na(rt))

  # summaries for each compound x sample

  if ("samplename_factor" %in% colnames(standard_peaks)) {
    max_ic_plot <- ggplot(standard_peaks, aes(x = samplename_factor, y = ic, color = sample_type, group = name)) +
      geom_point() +
      facet_wrap(~name, scale = "free_y") +
      scale_x_discrete("Sample Name") +
      scale_y_continuous("Max ion count") +
      scale_color_manual("Sample Type", values = c("blank" = "azure3", "standard" = "blue2", "sample" = "forestgreen")) +
      expand_limits(y = 0) +
      theme_bw() +
      theme(text = element_text(size = 12), axis.text.x = element_text(angle = 70, hjust = 1), legend.position = "bottom")

    rt_consistency_plot <- standard_peaks %>%
      dplyr::group_by(name) %>%
      dplyr::mutate(fraction_max = ic / max(ic)) %>%
      dplyr::ungroup() %>%
      ggplot(aes(x = samplename_factor, y = rt, fill = sample_type, group = name)) +
      geom_point(aes(alpha = fraction_max), shape = 21, stroke = 0.1, size = 2) +
      geom_hline(data = standard_reference_values, aes(yintercept = rt, color = rt_type)) +
      facet_wrap(~name, scale = "free_y") +
      expand_limits(y = 0) +
      scale_x_discrete("Sample Name") +
      scale_y_continuous("Retention Time") +
      scale_color_brewer("RT type", palette = "Set2") +
      scale_alpha_continuous("Fraction of Max Signal for Standard") +
      scale_fill_manual("Sample Type", values = c("blank" = "azure3", "standard" = "blue2", "sample" = "forestgreen")) +
      theme_bw() +
      theme(text = element_text(size = 12), axis.text.x = element_text(angle = 70, hjust = 1), legend.position = "bottom")
  } else {
    max_ic_plot <- NULL
    rt_consistency_plot <- NULL
  }

  # add y offsets for visualization

  y_adjustments <- standard_peaks %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(
      compound_max_ic = max(ic),
      y_offset = compound_max_ic * 0.4 * (as.integer(samplename_factor) - 1) / (nrow(distinct_samples) - 1)
    ) %>%
    dplyr::select(samplename, name, y_offset)

  unnested_eic <- standard_data %>%
    tidyr::unnest(eic) %>%
    # add y offset
    dplyr::left_join(y_adjustments, by = c("samplename", "name")) %>%
    dplyr::mutate(ic_adjusted = ic + y_offset) %>%
    # filter to within 5 minutes of expected RT
    dplyr::filter(
      is.na(observedRt) | abs(rt - observedRt) < 5,
      is.na(annotatedRt) | abs(rt - annotatedRt) < 5
    ) %>%
    dplyr::left_join(distinct_samples, by = "samplename")

  eic_plot <- ggplot(data = unnested_eic, aes(x = rt, y = ic_adjusted, color = sample_type, group = samplename)) +
    geom_path(alpha = 0.5) +
    facet_wrap(~name, scale = "free") +
    scale_x_continuous("Retention Time") +
    scale_y_continuous("Ion Count (with samples offset)") +
    scale_color_manual("Sample Type", values = c("blank" = "azure3", "standard" = "blue2", "sample" = "forestgreen")) +
    theme_minimal() +
    theme(text = element_text(size = 12), axis.text.x = element_text(angle = 70, hjust = 1), legend.position = "bottom")

  tibble::tribble(
    ~plot_type, ~grob,
    "IC", max_ic_plot,
    "RT", rt_consistency_plot,
    "EIC", eic_plot
  )
}
