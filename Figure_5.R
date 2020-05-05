library(dplyr)
library(tidyr)
library(forcats)
library(readr)
library(ggplot2)

source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "ccldcccdd" ) -> tecan
read_tsv( "data/ngs_counts.tsv" ) -> ngs
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00003", "CP00005", "CP00006", "CP00008", 
                    "CP00009", "CP00010", "CP00011", "CP00012", "CP00013", "CP00016" )

lamp_thresholds <- c(.3)
qpcr_thresholds <- c(30, 42)
ngs_threshold <- 3000

## Figure 5B
tbl <- ngs %>%
  filter( !is.na(matchedTRUE)) %>%
  left_join( tecan ) %>% 
  mutate(minutes = ifelse(plate == "CP00012" & minutes == 45, 40, minutes)) %>%
  filter( plate %in% plates_to_use ) %>%
  filter( !(plate == "CP00003" & plateRemark != "2")) %>%
  filter( minutes %in% c(30, 40), gene=="N" ) %>% 
  filter( ! ( plate == "CP00005" & well >= "G" ) ) %>%
  left_join( tblCT )%>%
  filter( !is.na(CT) )

set.seed(42)
tbl %>%
  mutate( CT = ifelse( CT>40, runif( n(), 43, 47 ), CT ) ) %>%
  mutate( NGS = case_when(
    matchedTRUE > ngs_threshold ~ "positive",
    matchedTRUE <= ngs_threshold ~ "negative",
    TRUE                ~ "undetected")) %>%
  mutate_at( "NGS", fct_relevel, "positive", "negative" ) %>%
ggplot() +
  geom_hline(yintercept = lamp_thresholds, color = "lightgray" ) +
  geom_vline(xintercept = qpcr_thresholds, color = "darkgrey" ) +
  geom_point( aes( x = CT, y = absBlue - absYellow, fill = NGS ), colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
  #scale_x_continuous( breaks = c( 20, 30, 40, 45 ), labels = c( 20, 30, 40, "neg" ) ) +
  scale_x_continuous( breaks = c( 20, 30, 40, 45 ), labels = c( 20, 30, 40, "neg" ), trans = "reverse" ) +
  labs(subtitle = str_interp( "${nrow(filter(tbl, minutes == 30))} samples on ${length(unique(tbl$plate))} plates" ),
      x = "RT-qPCR (CT value)",
       y = "RT-LAMP (ΔOD)") +
  scale_fill_manual(name  = "LAMP-sequencing", values = c("positive" = "black", "negative" = "white")) +
  facet_grid(cols = vars(minutes), labeller = as_labeller(function(x) str_c(x, " min at 65°C"))) +
  #facet_grid(cols = vars(facets)) +
  annotate("text", x = 50, y = 0, label = str_glue("negative"), angle = 90, col="grey50") +
  #annotate("text", x = 50, y = -.26, label = str_glue("negative"), angle = 90, col="grey50") +
  #annotate("text", x = 50, y = .125, label = str_glue("inconclusive"), angle = 90, col="grey50") +
  annotate("text", x = 50, y = .425, label = str_glue("positive"), angle = 90, col="grey50") +
  coord_cartesian(xlim = c(11.75, 49.5)) +
  theme(plot.subtitle = element_text(hjust = .5), legend.position = "bottom")

# Export figures
ggsave("SVGs/Figure_5b.svg", width=20, height=10, units="cm")
ggsave("Figure_5b.png", width=20, height=10, units="cm", dpi=300)

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
    matchedTRUE > ngs_threshold ~ "pos",
    matchedTRUE <= ngs_threshold ~ "neg",
    TRUE                ~ "undet")) %>%
  mutate_at( "LAMPres", fct_relevel, "pos", "neg" ) %>%
  mutate_at( "NGSres", fct_relevel, "pos", "neg" ) %>%
  count( LAMPres, NGSres, CTbin ) %>%
  pivot_wider( names_from = LAMPres, values_from = n, values_fill = c(n=0) ) %>%
  arrange(NGSres, CTbin)
confusion_matrix

confusion_matrix %>% write_tsv( "NGS_confMatrix.tsv" )
