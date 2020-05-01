library( tidyverse )
library( patchwork )
source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "cclicccdd" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT


## Figure 3b 
tecan %>% 
left_join( tblCT ) %>%
mutate( CT = ifelse( CT>40, runif( n(), 43, 46 ), CT ) ) %>%
filter( !is.na(CT) ) %>%
filter( plate %in% "CP00001", gene=="N", plateRemark=="2" ) %>% 
ggplot +
  geom_line( aes( x=minutes, y=absBlue-absYellow, group=well, col=CT ) ) +
  scale_color_ct() +
  xlab( "incubation time [minutes]" ) + 
  ylab( expression( "LAMP signal (ΔOD)" ) ) -> plot3b

print(plot3b)

## Figure 3cd

tecan %>% 
filter( plate == "CP00001" ) %>%
filter( ( gene=="N" & minutes == 30 ) | ( gene=="1a" & minutes==40 ) ) %>%
left_join( tblCT ) %>%
filter( !is.na(CT) ) %>%
mutate( CT = ifelse( CT>40, runif( n(), 43, 46 ), CT ) ) %>%
mutate( facet = str_c( "gene ", gene, ", ", minutes, " min" ) ) %>%
mutate_at( vars(facet), fct_rev ) %>%
ggplot +
  geom_vline( xintercept = 41.5, col="darkgray" ) +
  geom_point( aes( x=CT, y=absBlue-absYellow ), size=.7, shape=19, alpha=.6, fill="lightgray" ) +
  scale_x_reverse( breaks = c( 20, 30, 40, 44.5 ), labels = c( 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value)" ) + 
  ylab( "RT-LAMP (ΔOD)" ) +
  facet_wrap( . ~ facet ) -> plot3c

( plot_spacer() | plot3b ) / plot3c
dev.copy( svg, "Figure_3.svg", width=9, height=6 )
dev.off()
