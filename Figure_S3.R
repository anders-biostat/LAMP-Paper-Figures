library(dplyr)
library(readr)
library(scales)
library(stringr)
library(ggplot2)
library(patchwork)
library(ggsci)

source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "ccliccccdd" ) -> tecan
read_tsv( "data/ngs_counts.tsv" ) -> ngs
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00003", "CP00005", "CP00006", "CP00008", 
                    "CP00009", "CP00010", "CP00011", "CP00012", "CP00013", "CP00016" )

lamp_thresholds <- c(.3)
qpcr_thresholds <- c(30, 42)
ngs_threshold <- 3000

## Figure S2
tbl <- ngs %>%
  filter( !is.na(matchedTRUE)) %>%
  left_join( tecan ) %>%
  filter( plate %in% plates_to_use ) %>%
  filter( !(plate == "CP00003" & plateRemark != "2")) %>%
  filter( minutes == 30, gene=="N" ) %>% 
  filter( ! ( plate == "CP00005" & well >= "G" ) ) %>%
  left_join( tblCT )%>%
  filter( !is.na(CT) ) %>%
  filter( is.na(wellRemark) )

set.seed(42)
panels_data <-  tbl %>%
  mutate( CT = ifelse( CT>40, runif( n(), 43, 47 ), CT ) ) %>%
  mutate(plate = str_extract(plate, "\\d{2}$")) %>%
  mutate(well = fct_reorder(well, desc(CT)))

## Figure S2a
panel_a <-  panels_data %>%
  ggplot() +
  geom_hline(yintercept = lamp_thresholds, color = "lightgray" ) +
  geom_vline(xintercept = ngs_threshold, color = "lightgray" ) +
  geom_point( aes( x = matchedTRUE, y = absBlue - absYellow, fill = CT, group = well ),
              colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
  scale_x_continuous(trans = "log10", labels = label_math(10^.x, format = log10) ) +
  labs(subtitle = str_interp( "${unique(tbl$minutes)} min at 65°C\n${length(tbl$minutes)} samples on ${length(unique(tbl$plate))} plates"), 
       x = expression(paste("Multiplexed sequencing (", log[10](CPM), ")")),
       y = expression(paste("RT-LAMP (", Delta, OD, ")"))) +
  scale_fill_ct( name="RT-qPCR\nCT value") +
  #annotate("text", x = 0.1, y = -.26, label = str_glue("negative"), angle = 90, col="grey50") +
  annotate("text", x = 0.1, y = 0, label = str_glue("negative"), angle = 90, col="grey50") +
  #annotate("text", x = 0.1, y = .125, label = str_glue("inconclusive"), angle = 90, col="grey50") +
  annotate("text", x = 0.1, y = .425, label = str_glue("positive"), angle = 90, col="grey50")
panel_a

# Figure S2b
panel_b <-  panels_data %>%
  ggplot() +
  geom_hline(yintercept = lamp_thresholds, color = "lightgray" ) +
  geom_vline(xintercept = ngs_threshold, color = "lightgray" ) +
  geom_point( aes( x = matchedTRUE, y = absBlue - absYellow, fill = plate, group = well ),
              colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
  scale_x_continuous(trans = "log10", labels = label_math(10^.x, format = log10) ) +
  labs(subtitle = str_interp( "${unique(tbl$minutes)} min at 65°C\n${length(tbl$minutes)} samples on ${length(unique(tbl$plate))} plates"), 
       x = expression(paste("Multiplexed sequencing (", log[10](CPM), ")")),
       y = expression(paste("RT-LAMP (", Delta, OD, ")"))) +
  scale_fill_d3(palette="category20", guide = guide_legend(label = FALSE)) +
  #annotate("text", x = 0.1, y = -.26, label = str_glue("negative"), angle = 90, col="grey50") +
  annotate("text", x = 0.1, y = 0, label = str_glue("negative"), angle = 90, col="grey50") +
  #annotate("text", x = 0.1, y = .125, label = str_glue("inconclusive"), angle = 90, col="grey50") +
  annotate("text", x = 0.1, y = .425, label = str_glue("positive"), angle = 90, col="grey50")
panel_b

panel_a + panel_b +
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(face = "bold"))

# tbl %>%
#   ggplot() +
#   geom_hline(yintercept = ngs_threshold, color = "lightgray" ) +
#   geom_point( aes( x = matchedFALSE, y = matchedTRUE , fill = matchedTRUE > 3000 ),
#               colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
#  # scale_fill_ct( name="RT-qPCR\nCT value") +
#   scale_x_continuous(trans = "log10", labels = label_math(10^.x, format = log10) ) +
#   scale_y_continuous(trans = "log10", labels = label_math(10^.x, format = log10) )

# Export figures
ggsave("SVGs/Figure_S3.svg", width=22, height=10, units="cm")
ggsave("Figure_S3.png", width=22, height=10, units="cm", dpi=300)
