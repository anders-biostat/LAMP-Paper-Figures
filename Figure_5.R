library(dplyr)
library(readr)
library(ggplot2)

read_tsv( "data/tecan_values.tsv", guess_max = 1e4 ) -> tecan
read_tsv( "data/ngs_counts.tsv" ) -> ngs
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

lamp_thresholds <- c(-.15, .25)
qpcr_thresholds <- c(30, 42)
ngs_threshold <- 3000

## Figure 5B
tbl <- ngs %>%
  left_join( tecan ) %>% 
  left_join( tblCT ) %>%
  filter( gene == "N" ) %>%
  filter( !is.na(CT) ) %>%
  filter( is.na(plateRemark) | plateRemark == "A" ) %>% 
  filter( is.na(wellRemark) ) %>%
  filter( minutes %in% c(30, 40) ) %>%
  filter( plate > "CP00004" , plate <= "CP00016") %>%
  filter( !(plate %in% c("CP00014", "CP00015")) ) %>%   ## this could be omited if these plates would be annotated as "bead" in 'plateRemark' column 
  mutate( CT = ifelse( CT>40, runif( n(), 43, 46 ), CT ) ) %>%
  mutate( NGS = case_when(
    matchedTRUE > ngs_threshold ~ "positive",
    matchedTRUE <= ngs_threshold ~ "negative",
    TRUE                ~ "undetected")) %>%
  mutate_at( "NGS", fct_relevel, "positive", "negative" ) %>%
  mutate( minutes = str_c(minutes, " min"))
ggplot(tbl) +
  geom_hline(yintercept = lamp_thresholds, color = "lightgray" ) +
  geom_vline(xintercept = qpcr_thresholds, color = "lightgray" ) +
  geom_point( aes( x = CT, y = absBlue - absYellow, fill = NGS ), colour = "black", alpha = .6, shape = 21, size = 2 ) + 
  scale_x_continuous( breaks = c( 20, 30, 40, 45 ), labels = c( 20, 30, 40, "neg" ) ) +
  labs(title = str_glue("(n = {nrow(filter(tbl, minutes == '30 min'))})"),
       x = "RT-qPCR\nCT value gene E",
       y = "LAMP assay\nΔOD(434nm - 560nm)") +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(cols = vars(minutes)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank()) +
  annotate("text", x = 13, y = -.30, label = str_glue("negative"), angle = 90) +
  annotate("text", x = 13, y = .06, label = str_glue("inconclusive"), angle = 90) +
  annotate("text", x = 13, y = .4, label = str_glue("positive"), angle = 90)

dev.copy( svg, "figs/Figure_5b.svg", width=5, height=4 )
dev.off()
