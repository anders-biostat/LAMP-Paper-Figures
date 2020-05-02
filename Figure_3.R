library( tidyverse )
library( patchwork )
source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "cclicccdd" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT


## Figure 3b 
tecan %>% 
left_join( tblCT ) %>%
mutate( CT = ifelse( CT>40, 46.5, CT ) ) %>%
filter( !is.na(CT) ) %>%
filter( plate %in% "CP00001", gene=="N", plateRemark=="2" ) %>% 
ggplot +
  geom_line( aes( x=minutes, y=absBlue-absYellow, group=well, col=CT ), size=.5 ) +
  scale_color_ct() +
  xlab( "incubation time [minutes]" ) + 
  ylab( expression( "RT-LAMP (ΔOD)" ) ) -> plot3b

plot3b

## Figure 3c
tbl_3c <- tecan %>% 
  filter( plate == "CP00001" ) %>%
  filter( ( gene=="N" & minutes == 30 ) | ( gene=="1a" & minutes==40 ) ) %>%
  left_join( tblCT ) %>%
  filter( !is.na(CT) ) %>%
  mutate( CT = ifelse( CT>40, runif( n(), 43, 47 ), CT ) ) %>%
  mutate( facet = str_c( gene, " gene, ", minutes, " min" ) )

tbl_3c %>%
  filter(gene == "N") %>%
  ggplot() +
  geom_vline( xintercept = 42, col="darkgray" ) +
  geom_point( aes( x=CT, y=absBlue-absYellow ), size=.7, shape=19, alpha=.6, fill="lightgray" ) +
  scale_x_reverse( breaks = c( 10, 20, 30, 40, 45 ), labels = c( 10, 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value)" ) + 
  ylab( "RT-LAMP (ΔOD)" ) -> plot3c

## Supplementary Figure S3
tbl_3c %>%
  filter(gene == "1a") %>%
  ggplot() +
  geom_vline( xintercept = 42, col="darkgray" ) +
  geom_point( aes( x=CT, y=absBlue-absYellow ), size=.7, shape=19, alpha=.6, fill="lightgray" ) +
  scale_x_reverse( breaks = c( 10, 20, 30, 40, 45 ), labels = c( 10, 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value)" ) + 
  ylab( "RT-LAMP (ΔOD)" ) -> plotS3

ggplot() / (plot3b | plot3c) +
  plot_annotation(tag_levels = "a")

ggsave("SVGs/Figure_3.svg", width=20, height=14, units="cm")
ggsave("Figure_3.png", width=20, height=14, units="cm", dpi=300)


plotS3

ggsave("SVGs/Figure_S3.svg", width=9, height=7, units="cm")
ggsave("Figure_S3.png", width=9, height=7, units="cm", dpi=300)
