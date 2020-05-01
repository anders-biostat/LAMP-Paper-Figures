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
  mutate( CT = ifelse( CT>40, runif( n(), 43, 47 ), CT ) ) %>%
  mutate( NGS = case_when(
    matchedTRUE > ngs_threshold ~ "positive",
    matchedTRUE <= ngs_threshold ~ "negative",
    TRUE                ~ "undetected")) %>%
  mutate_at( "NGS", fct_relevel, "positive", "negative" )

ggplot(tbl) +
  geom_hline(yintercept = lamp_thresholds, color = "lightgray" ) +
  geom_vline(xintercept = qpcr_thresholds, color = "darkgrey" ) +
  geom_point( aes( x = CT, y = absBlue - absYellow, fill = NGS ), colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
  scale_x_continuous( breaks = c( 20, 30, 40, 45 ), labels = c( 20, 30, 40, "neg" ) ) +
  #scale_x_continuous( breaks = c( 20, 30, 40, 45 ), labels = c( 20, 30, 40, "neg" ), trans = "reverse" ) +
  labs(title = str_interp( "30 min at 65°C\n${nrow(tbl)} samples on ${tbl%>%select(plate)%>%unique%>%nrow} plates" ),
       x = "RT-qPCR (CT value)",
       y = "RT-LAMP assay (ΔOD)") +
  scale_fill_manual(values = c("positive" = "black", "negative" = "white")) +
  #facet_grid(cols = vars(minutes), labeller = as_labeller(function(x) str_c(x, " min at 65°C"))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank()) +
  annotate("text", x = 11.25, y = -.26, label = str_glue("negative"), angle = 90) +
  annotate("text", x = 11.25, y = .125, label = str_glue("inconclusive"), angle = 90) +
  annotate("text", x = 11.25, y = .425, label = str_glue("positive"), angle = 90) +
  coord_cartesian(xlim = c(11.75, 47.5))

# Export figures
ggsave("figs/Figure_5b.svg", width=14, height=10, units="cm")
ggsave("figs/Figure_5b.png", width=14, height=10, units="cm", dpi=300)
