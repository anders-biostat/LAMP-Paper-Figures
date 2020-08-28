library( tidyverse )
library( patchwork )
source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "ccliccccdd" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

set.seed(42)
tecan %>% 
  filter( plate == "CP00001" ) %>%
  filter( ( gene=="1a" & minutes==40 ) ) %>%
  left_join( tblCT ) %>%
  filter( !is.na(CT) ) %>%
  mutate( CT = ifelse( CT>40, runif( n(), 43, 47 ), CT ) ) %>%
  ggplot() +
  geom_vline( xintercept = 42, col="darkgray" ) +
  geom_point( aes( x = CT, y = absBlue - absYellow ), fill = "black", colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
  scale_x_reverse( breaks = c( 10, 20, 30, 40, 45 ), labels = c( 10, 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value)" ) + 
  ylab( "RT-LAMP (ΔOD)" ) -> fig
fig

ggsave("SVGs/Figure_S3.svg", width=9, height=7, units="cm")
ggsave("Figure_S3.png", width=9, height=7, units="cm", dpi=300)