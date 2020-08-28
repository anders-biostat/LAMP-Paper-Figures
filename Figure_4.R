library(dplyr)
library(tidyr)
library(forcats)
library(readr)
library(ggplot2)
library(stringr)
library(patchwork)

source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "ccldccccdd" ) -> tecan
read_tsv( "data/ngs_counts.tsv" ) -> ngs
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00003", "CP00005", "CP00006", "CP00008", 
                    "CP00009", "CP00010", "CP00011", "CP00012", "CP00013", "CP00016" )

lamp_thresholds <- c(.3)
qpcr_thresholds <- c(30, 42)
ngs_thresholds <- c(200, 1e4)

panel_a <- rsvg::rsvg("SVGs/Figure_4a.svg")

tbl <- ngs %>%
  full_join( tecan ) %>%
  filter( plate %in% plates_to_use ) %>%
    mutate(minutes = ifelse(plate == "CP00012" & minutes == 45, 40, minutes)) %>%
  filter( !(plate == "CP00003" & plateRemark != "2")) %>% 
  filter( minutes %in% c(30, 40), gene=="N" ) %>% 
  filter( ! ( plate == "CP00005" & well >= "G" ) ) %>%
  left_join( tblCT )%>%
  filter( !is.na(CT) ) %>%
  filter( is.na(wellRemark) ) %>%
  replace_na(list(matchedTRUE = 0, matchedFALSE = 0)) 


set.seed(2020)
panels_data <-  tbl %>%
  mutate( NGS = case_when(
    matchedTRUE + matchedFALSE <= ngs_thresholds[[1]] ~ "too few UMIs",
    matchedTRUE > ngs_thresholds[[2]]                 ~ "positive",
    TRUE                                              ~ "negative")) %>%
  mutate(CT = ifelse( CT>40, runif( n(), 43, 47 ), CT ) ) %>%
  mutate(well = fct_reorder(well, -CT))

lamp_colors <- c("positive" = "#00D302", "negative" = "#C7007C", "too few UMIs" = "black")

panel_b <- panels_data %>%
  arrange(NGS) %>%
ggplot() +
  geom_hline(yintercept = lamp_thresholds, color = "lightgray" ) +
  geom_vline(xintercept = qpcr_thresholds, color = "lightgrey" ) +
  geom_point( aes( x = CT, y = absBlue - absYellow, fill = NGS ), colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
  scale_x_continuous( breaks = c( 20, 30, 40, 45 ), labels = c( 20, 30, 40, "neg" ), trans = "reverse" ) +
  labs(subtitle = str_interp( "${nrow(filter(tbl, minutes == 30))} samples on ${length(unique(tbl$plate))} plates" ),
      x = "RT-qPCR (CT value)",
       y = "RT-LAMP (ΔOD)") +
  scale_fill_manual(name  = "LAMP-sequencing:", values = lamp_colors) +
  facet_grid(cols = vars(minutes), labeller = as_labeller(function(x) str_c(x, " min at 65°C"))) +
  annotate("text", x = 50, y = 0, label = str_glue("negative"), angle = 90, col="grey50") +
  annotate("text", x = 50, y = .425, label = str_glue("positive"), angle = 90, col="grey50") +
  coord_cartesian(xlim = c(11.75, 49.5)) +
  theme(plot.subtitle = element_text(hjust = .5), legend.position = "bottom")

fig_layout <- '
A
B
'
wrap_elements(plot =  grid::rasterGrob(panel_a)) +
  panel_b +
  plot_layout(design = fig_layout) +
  plot_annotation(tag_levels = "A")
  
# Export figures
ggsave("SVGs/Figure_4tmp.svg", width=20, height=22, units="cm")
ggsave("Figure_4tmp.png", width=20, height=22, units="cm", dpi=300)

# count table S3
ct_breaks <- c(0, 25, 30, 35, 40, Inf)

confusion_matrix <- tbl %>%
  filter(minutes == 30) %>%
  mutate( CTbin = cut( CT, ct_breaks ) ) %>% 
  mutate_at( "CTbin", recode, 	"neg" = "(40,Inf]" ) %>%
  mutate( LAMPres = 
            cut( absBlue-absYellow, c( -Inf, lamp_thresholds, Inf ) ) %>%
            as.integer %>%
            { c( "neg", "pos" )[.] } ) %>%
  mutate( NGSres = case_when(
    matchedTRUE + matchedFALSE <= ngs_thresholds[[1]] ~ "too_few",
    matchedTRUE > ngs_thresholds[[2]]                 ~ "pos",
    TRUE                                              ~ "neg")) %>%
  mutate_at( "LAMPres", fct_relevel, "pos", "neg" ) %>%
  mutate_at( "NGSres", fct_relevel, "pos", "neg", "too_few" ) %>%
  count( LAMPres, NGSres, CTbin ) %>%
  pivot_wider( names_from = LAMPres, values_from = n, values_fill = c(n=0) ) %>%
  arrange(NGSres, CTbin)
confusion_matrix

confusion_matrix %>% write_tsv( "NGS_confMatrix.tsv" )

