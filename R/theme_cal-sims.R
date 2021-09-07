#' [ggplot2] theme for cal-sims
#'
theme_cal <- function() {
  ggplot2::theme_bw() +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "white"),
    strip.text = ggplot2::element_text(colour = "black")
  )
}
