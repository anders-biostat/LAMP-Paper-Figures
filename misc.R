# Use these color scales for CT-values.
# You can still change the breaks manually by passing the breaks argument
scale_color_ct <- function(...) {
  scale_color_gradientn(
    ...,
    guide = guide_colourbar(reverse = TRUE),
    breaks = rev(c(10, 15, 20, 25, 30, 35, 40, 45)),
    labels = rev(c(10, 15, 20, 25, 30, 35, 40, "neg")),
    colors = c(viridis::magma(100, direction=1L)[1:75],         # range of CP values
               rep("#FFFFFF", 5),                               # white spacer
               viridis::magma(100, direction=1L)[rep(90, 10)])  # negative
  )
}

scale_fill_ct <- function(...) {
  scale_fill_gradientn(
    ...,
    guide = guide_colourbar(reverse = TRUE),
    breaks = rev(c(10, 15, 20, 25, 30, 35, 40, 45)),
    labels = rev(c(10, 15, 20, 25, 30, 35, 40, "neg")),
    colors = c(viridis::magma(100, direction=1L)[1:75],         # range of CP values
               rep("#FFFFFF", 5),                               # white spacer
               viridis::magma(100, direction=1L)[rep(90, 10)])  # negative
  )
}