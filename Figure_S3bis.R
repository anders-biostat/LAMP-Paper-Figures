library(tidyverse)
library(patchwork)

source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "ccliccccdd" ) -> tecan
read_tsv( "data/ngs_counts.tsv" ) -> ngs
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00003", "CP00005", "CP00006", "CP00008", 
                    "CP00009", "CP00010", "CP00011", "CP00012", "CP00013", "CP00016" )

lamp_thresholds <- c(.3)
qpcr_thresholds <- c(30, 42)
ngs_threshold <- c(100, 1e4)

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
    matchedTRUE > ngs_threshold[[2]] ~ "positive",
    between(matchedTRUE, ngs_threshold[[1]], ngs_threshold[[2]]) ~ "negative",
    matchedTRUE <= ngs_threshold[[1]]           ~ "too_low")) %>%
  mutate(CT = ifelse( CT>40, runif( n(), 43, 47 ), CT ) ) %>%
  mutate(plate = str_extract(plate, "\\d{2}$")) %>%
  mutate(well = fct_reorder(well, -CT))

mybreaks <- c( 0, 1, 10, 100, 1000, 1e4, 1e5, 1e6 )
mylabels <- c( 0, 1, 10, expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6) )
panel_b <- panels_data %>%
  mutate(NGSshape = if_else(NGS != "undetected", "detected", NGS)) %>%
  ggplot(aes(x = matchedTRUE + matchedFALSE, y = matchedTRUE)) +
  geom_vline(xintercept = ngs_threshold[[1]], color = "lightgray" ) +
  geom_hline(yintercept = ngs_threshold[[2]], color = "lightgray" ) +
  geom_point(aes( col = CT ), alpha = .6, size = 1.2 ) +
  scale_x_continuous( trans=logap_trans(), breaks = mybreaks, labels=mylabels ) +
  scale_y_continuous( trans=logap_trans(), breaks = mybreaks, labels=mylabels ) +
  scale_color_ct() +
  labs(x = "total read count",
       y = "virus-matching read count") +
  annotate("text", color = "gray40", x = 3e6, y = 1e5, label = "positive", angle = 90) +
  annotate("text", color = "gray40", x = 3e6, y = 1e3, label = "negative", angle = 90) +
  annotate("text", color = "gray70", x = 1e6, y = 1e4, label = "LAMP-seq result:", angle = 90) +
  annotate("text", color = "gray40", x = 20, y = 0, label = "too low", angle = 0) +
  annotate("text", color = "gray40", x = 3e2, y = 0, label = "ok", angle = 0) +
  annotate("text", color = "gray70", x = 1e2, y = 1, label = "read count", angle = 0) +
  coord_fixed()
panel_b


panel_c <- 
panels_data %>%
  ggplot(aes(x = absBlue-absYellow, y = matchedTRUE)) +
  geom_point(aes(col=NGS), alpha = .6, size = 1.2 ) +
  scale_y_continuous( trans=logap_trans(), breaks = mybreaks, labels=mylabels ) +
  labs(y = "LAMP-sequencing: virus-matching read count",
       x = "RT-LAMP: ΔOD", col="LAMP-sequencing") +
  geom_hline(yintercept = ngs_thresholds[[2]], color = "lightgray" ) +
  geom_vline(xintercept = lamp_thresholds, color = "lightgray" ) +
  annotate("text", color = "gray40", x = .2, y = 0, label = "negative", angle = 0) +
  annotate("text", color = "gray40", x = .4, y = 0, label = "positive", angle = 0) +
  annotate("text", color = "gray70", x = .3, y = 1, label = "LAMP result", angle = 0) +
  annotate("text", color = "gray40", x = .6, y = 1e3, label = "negative", angle = 90) +
  annotate("text", color = "gray40", x = .6, y = 1e5, label = "positive", angle = 90) +
  annotate("text", color = "gray40", x = .6, y = 1e1, label = "too few reads", angle = 90) +
  annotate("text", color = "gray70", x = .55, y = 1e3, label = "LAMP-sequencing result", angle = 90) 
panel_c


panel_a <- wrap_elements(panel = grid::textGrob('Here goes the figrue with the primers')) +
  theme(plot.margin = margin(t = 15, r = 5.5, b = 5.5, l = 15, "pt"))


fig_layout <- '
AAAAA
BBBCC
BBBCC
BBBCC
'
wrap_elements(full = panel_a + theme(plot.margin = margin(t = 15, r = 5.5, b = 5.5, l = 15, "pt"))) +
  panel_b + panel_c + 
  plot_layout(design = fig_layout) +
  plot_annotation(tag_levels = "a")

# Export figures
ggsave("SVGs/Figure_S3bis.svg", width=20, height=20, units="cm")
ggsave("Figure_S3bis.png", width=20, height=23, units="cm", dpi=300)

