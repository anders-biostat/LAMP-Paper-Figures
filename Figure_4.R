library(tidyverse)
library(patchwork)
library(ggtext)
library(glue)
library(ggsci)

read_tsv( "data/tecan_values.tsv" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00001", "CP00003", "CP00005", "CP00006", "CP00008", 
   "CP00009", "CP00010", "CP00011", "CP00012", "CP00013", "CP00016" )

lamp_thresh <- c(-.05, .3)

## Figure 4a
tecan %>% 
  filter( plate %in% plates_to_use ) %>%
  filter( minutes==30, gene=="N" ) %>% 
  filter( ! ( plate == "CP00005" & well >= "G" ) ) %>%
  left_join( tblCT ) %>% 
  filter( !is.na(CT) ) -> tbl

fig4a <- tbl %>%
  mutate(CT = ifelse(CT > 40, runif(n(), 43, 47), CT)) %>%
  ggplot() +
  geom_vline(xintercept = c(30, 41.5), color = "darkgray" ) +
  geom_hline(yintercept = lamp_thresh, color = "lightgray" ) +
  geom_point(aes(x = CT, y = absBlue - absYellow, color = plate), size = .7) +
  scale_x_continuous(breaks = c(20, 30, 40, 45), labels = c(20, 30, 40, "negative")) +
  labs(title = glue("30 min at 65°C\n{nrow(tbl)} samples on {tbl%>%distinct(plate)%>%nrow} plates"),
       x = glue("RT-qPCR (CT value)"),
       y = "RT-LAMP assay (ΔOD)") +
  scale_color_d3(palette="category20") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #annotate("text", x = 15.5, y = -.42, label = glue("n = {nrow(lamp_ct_fig4)}")) +
  annotate("text", x = 13, y = -.25, label = glue("negative"), angle = 90) +
  annotate("text", x = 13, y = .105, label = glue("inconclusive"), angle = 90) +
  annotate("text", x = 13, y = .4, label = glue("positive"), angle = 90) +
  coord_cartesian(xlim = c(14, 49))

fig4a


## Confusion matrix
tbl %>%
mutate( CTbin = cut( CT, c( 0, 25, 30, 35, 40, Inf ) ) ) %>% 
mutate_at( "CTbin", recode, 	"neg" = "(40,Inf]" ) %>%
mutate( LAMPres = 
   cut( absBlue-absYellow, c( -Inf, lamp_thresh, Inf ) ) %>%
   as.integer %>%
   { c( "neg", "incl", "pos" )[.] } ) %>%
group_by( CTbin, LAMPres ) %>%
count() %>%
pivot_wider( names_from = LAMPres, values_from = n, values_fill = c(n=0) )









