#' Louis' preferred ggplot2 theme
#'
#'`theme_le` may be used as a `theme()` function. The theme is a slightly modified version of `theme_classic`.
#' Changes included an increased default font size and default serif font.
#'
#' @param ... see `theme()`



theme_le <- function(...){
  theme_classic() +
    theme(...,
    text = element_text(family = "serif", size = 16)
  )
}
