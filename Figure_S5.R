library(tidyverse)
library(patchwork)
library(ggrepel)

source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "ccliccccdd" ) -> tecan
read_tsv( "data/ngs_counts.tsv" ) -> ngs
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00003", "CP00005", "CP00006", "CP00008", 
                    "CP00009", "CP00010", "CP00011", "CP00012", "CP00013", "CP00016" )

lamp_thresholds <- c(.3)
qpcr_thresholds <- c(30, 42)
ngs_thresholds <- c(200, 1e4)

lamp_colors <- c("positive" = "#00D302", "negative" = "#C7007C", "too few" = "black")

tbl <- ngs %>%
  full_join( tecan ) %>%
  filter( plate %in% plates_to_use ) %>%
  filter( !(plate == "CP00003" & plateRemark != "2")) %>% 
  filter( minutes == 30, gene=="N" ) %>% 
  filter( ! ( plate == "CP00005" & well >= "G" ) ) %>%
  left_join( tblCT )%>%
  filter( !is.na(CT) ) %>%
  filter( is.na(wellRemark) ) %>%
  replace_na(list(matchedTRUE = 0, matchedFALSE = 0)) 

set.seed(42)
panels_data <-  tbl %>%
  mutate( NGS = case_when(
    matchedTRUE + matchedFALSE <= ngs_thresholds[[1]] ~ "too few",
    matchedTRUE > ngs_thresholds[[2]]                 ~ "positive",
    TRUE                                              ~ "negative")) %>%
  mutate(CT = ifelse( CT>40, runif( n(), 43, 47 ), CT ) ) %>%
  mutate(well = fct_reorder(well, -CT))

mybreaks <- c( 0, 1, 10, 100, 1000, 1e4, 1e5, 1e6 )
mylabels <- c( 0, 1, 10, expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6) )
panel_a <- panels_data %>%
  mutate(NGSshape = if_else(NGS != "undetected", "detected", NGS)) %>%
  ggplot(aes(x = matchedTRUE + matchedFALSE, y = matchedTRUE)) +
  geom_vline(xintercept = ngs_thresholds[[1]], color = "lightgray" ) +
  geom_hline(yintercept = ngs_thresholds[[2]], color = "lightgray" ) +
  geom_point(aes( fill = CT ), colour = "black", alpha = .6, size = 2.2, shape = 21 ) +
  geom_text_repel(aes(label = well), data = filter(panels_data, plate == "CP00012", str_detect(well, "^A")), size = 2.5, nudge_x = -.1, nudge_y = .6) +
  scale_x_continuous( trans=logap_trans(), breaks = mybreaks, labels=mylabels ) +
  scale_y_continuous( trans=logap_trans(), breaks = mybreaks, labels=mylabels ) +
  scale_fill_ct( name="RT-qPCR\nCT value") +
  labs(x = "total UMI count",
       y = "virus-matching UMI count") +
  annotate("text", x = 3e6, y = 5e2, label = "negative", angle = 90, colour = lamp_colors[["negative"]]) +
  annotate("text", x = 3e6, y = 1e5, label = "positive", angle = 90, colour = lamp_colors[["positive"]]) +
  annotate("text", color = "gray70", x = 1e6, y = 4e3, label = "LAMP-sequencing result", angle = 90) +
  annotate("text", color = "gray40", x = 20, y = 0, label = "too few", angle = 0) +
  annotate("text", color = "gray40", x = 3e3, y = 0, label = "sufficient", angle = 0) +
  annotate("text", color = "gray70", x = 2e2, y = 1, label = "UMI count", angle = 0) 
panel_a

panel_b <- rsvg::rsvg("SVGs/Figure_S5b.svg")

panel_c <- 
panels_data %>%
  ggplot(aes(x = absBlue-absYellow, y = matchedTRUE)) +
  geom_point(aes(fill=NGS), alpha = .6, size = 2.2, shape = 21 ) +
  scale_y_continuous( trans=logap_trans(), breaks = mybreaks, labels=mylabels ) +
  labs(y = "virus-matching UMI count",
       x = "ΔOD", col="LAMP-sequencing") +
  scale_fill_manual(name  = "LAMP-sequencing", values = lamp_colors) +
  geom_hline(yintercept = ngs_thresholds[[2]], color = "lightgray" ) +
  geom_vline(xintercept = lamp_thresholds, color = "lightgray" ) +
  annotate("text", color = "gray40", x = .15, y = 0, label = "negative", angle = 0) +
  annotate("text", color = "gray40", x = .45, y = 0, label = "positive", angle = 0) +
  annotate("text", color = "gray70", x = .3, y = 1, label = "RT-LAMP result", angle = 0) +
  annotate("text", x = .62, y = 1e3, label = "negative", angle = 90, colour = lamp_colors[["negative"]]) +
  annotate("text", x = .62, y = 1e5, label = "positive", angle = 90, colour = lamp_colors[["positive"]]) +
  annotate("text", x = .62, y = 1e1, label = "too few UMIs", angle = 90, colour = lamp_colors[["too few"]]) +
  annotate("text", color = "gray70", x = .55, y = 1e3, label = "LAMP-sequencing result", angle = 90) +
  theme(legend.position = "none")
#  theme(legend.position = c(0.3, 0.85), legend.background = element_blank(), legend.box.background = element_blank(), legend.key=element_blank())
panel_c

fig_layout <- '
AB
CX'
panel_a + 
  wrap_elements(plot =  grid::rasterGrob(panel_b))  +
  panel_c + plot_spacer() +
  plot_layout(design = fig_layout) +
  plot_annotation(tag_levels = "A")

# Export figures
ggsave("SVGs/Figure_S5_tmp.svg", width=20, height=18, units="cm")
ggsave("Figure_S5_tmp.png", width=20, height=18, units="cm", dpi=300)

