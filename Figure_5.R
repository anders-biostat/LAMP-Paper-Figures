library(dplyr)
library(readr)
library(ggplot2)

read_tsv( "data/tecan_values.tsv", guess_max = 1e4 ) -> tecan
read_tsv( "data/ngs_counts.tsv" ) -> ngs
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00001", "CP00003", "CP00005", "CP00006", "CP00008", 
                    "CP00009", "CP00010", "CP00011", "CP00012", "CP00013", "CP00016" )

lamp_thresholds <- c(-.05, .3)
qpcr_thresholds <- c(30, 42)
ngs_threshold <- 3000



## Figure 5B
tbl <- ngs %>%
  filter( !is.na(matchedTRUE)) %>%
  left_join( tecan ) %>% 
  filter( plate %in% plates_to_use ) %>%
  filter( minutes==30, gene=="N" ) %>% 
  filter( ! ( plate == "CP00005" & well >= "G" ) ) %>%
  left_join( tblCT )%>%
  filter( !is.na(CT) ) %>%
  mutate( CT = ifelse( CT>40, runif( n(), 43, 46 ), CT ) ) %>%
  mutate( NGS = case_when(
    matchedTRUE > ngs_threshold ~ "positive",
    matchedTRUE <= ngs_threshold ~ "negative",
    TRUE                ~ "undetected")) %>%
  mutate_at( "NGS", fct_relevel, "positive", "negative" )
ggplot(tbl) +
  geom_hline(yintercept = lamp_thresholds, color = "lightgray" ) +
  geom_vline(xintercept = qpcr_thresholds[[1]], color = "lightgray" ) +
  geom_vline(xintercept = qpcr_thresholds[[2]], color = "grey20" ) +
  geom_point( aes( x = CT, y = absBlue - absYellow, fill = NGS ), colour = "black", alpha = .6, shape = 21, size = 2 ) + 
  #scale_x_continuous( breaks = c( 20, 30, 40, 45 ), labels = c( 20, 30, 40, "neg" ) ) +
  scale_x_continuous( breaks = c( 20, 30, 40, 45 ), labels = c( 20, 30, 40, "neg" ), trans = "reverse" ) +
  labs(subtitle = str_interp( "${nrow(tbl)} samples on ${tbl%>%select(plate)%>%unique%>%nrow} plates" ),
       x = expression(paste("RT-qPCR (CT value)")),
       y = expression(paste("RT-LAMP assay (", Delta, OD, ")"))) +
  scale_fill_manual(values = c("positive" = "black", "negative" = "white")) +
  #facet_grid(cols = vars(minutes), labeller = as_labeller(function(x) str_c(x, " min at 65°C"))) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank()) +
  annotate("text", x = 11, y = -.30, label = str_glue("negative"), angle = 90) +
  annotate("text", x = 11, y = .06, label = str_glue("inconclusive"), angle = 90) +
  annotate("text", x = 11, y = .4, label = str_glue("positive"), angle = 90)

dev.copy( svg, "figs/Figure_5b.svg", width=8, height=6 )
dev.off()
