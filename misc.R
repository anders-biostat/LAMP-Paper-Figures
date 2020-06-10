# Sets the default ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(family = 'Arial'),
                             panel.grid.minor = element_blank(),  # remove axis grid
                             panel.grid.major = element_blank(),
                             strip.background = element_blank(),  # remove silly gray bg of facet title
                             strip.text = element_text(hjust = 0), # left-justify facets
                             plot.tag = element_text(size = 20, face = "bold")
))


# Use these color scales for CT-values.
# You can still change the breaks manually by passing the breaks argument
scale_color_ct <- function(...) {
  scale_color_gradientn(
    ...,
    guide = guide_colourbar(reverse = TRUE),
    limits  = c(8, 47),
    breaks = rev(c(10, 15, 20, 25, 30, 35, 40, 45)),
    labels = rev(c(10, 15, 20, 25, 30, 35, 40, "neg")),
    colors = c(viridis::plasma(100, direction=-1L)[1:85],      # range of CP values
               rep("#FFFFFF", 5),                               # white spacer
               rep("grey55", 15))                               # negative 
  )
}

scale_fill_ct <- function(...) {
  scale_fill_gradientn(
    ...,
    guide = guide_colourbar(reverse = TRUE),
    limits  = c(8, 47),
    breaks = rev(c(10, 15, 20, 25, 30, 35, 40, 45)),
    labels = rev(c(10, 15, 20, 25, 30, 35, 40, "neg")),
    colors = c(viridis::plasma(100, direction=-1L)[1:85],      # range of CP values
               rep("#FFFFFF", 5),                               # white spacer
               rep("grey55", 15))                               # negative 
  )
}

# considered log-scale
logap_trans <- function( a=0 ) {
  scales::trans_new("logap",
                    function(x) log10( x + 10^a), 
                    function(x) 10^x - 10^a)
                    #format = scales::label_math(10^.x, format = log10)
                    # ^^ This here does not work, removes zeros
}

# reverse complement
rc <- function(x) stringi::stri_reverse(stringi::stri_trans_char(x, "GATC", "CTAG"))

