library( tidyverse )

read_tsv( "data/tecan_values.tsv" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT


## Figure 3b 
tecan %>% 
left_join( tblCT ) %>%
mutate( CT = ifelse( CT>40, runif( n(), 43, 46 ), CT ) ) %>%
filter( !is.na(CT) ) %>%
filter( plate %in% "CP00001", gene=="N", plateRemark=="2" ) %>% 
ggplot +
  geom_line( aes( x=minutes, y=absBlue-absYellow, group=well, col=CT ) ) +
  scale_color_gradientn( 
    colours=viridis::magma(100)[ c( 1:80, rep(90, 10 ) ) ],
    breaks = c( 20, 30, 40, 46 ), labels = c( 20, 30, 40, "neg" ) ) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
  xlab( "incubation time [minutes]" ) + 
  ylab( expression( OD["434 nm"] -  OD["560 nm"] ) ) +
  ggtitle( " " )

dev.copy( svg, "figs/Figure_3b.svg", width=3.8, height=2.5 )
dev.off()

## Figure 3c

tecan %>% 
filter( plate == "CP00001", minutes > 30 ) %>%
left_join( tblCT ) %>%
filter( !is.na(CT) ) %>%
mutate( CT = ifelse( CT>40, runif( n(), 43, 46 ), CT ) ) %>%
ggplot +
  geom_vline( xintercept = 41.5, col="darkgray" ) +
  geom_point( aes( x=CT, y=absBlue-absYellow, col=gene ) ) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + 
  scale_x_continuous( breaks = c( 20, 30, 40, 44.5 ), labels = c( 20, 30, 40, "neg" ) ) +
  xlab( "CT value from qPCR" ) + 
  ylab( expression( OD["434 nm"] -  OD["560 nm"] ) ) 

dev.copy( svg, "figs/Figure_3c.svg", width=3.8, height=2.5 )
dev.off()
  

## Figure 3d

tecan %>% 
filter( plate %in% c( "CP00001" ), gene=="1a", minutes == 40 ) %>%
left_join( tblCT ) %>%
mutate( CT = ifelse( CT>40, runif( n(), 43, 46 ), CT ) ) %>%
ggplot +
  geom_vline( xintercept = 41.5, col="darkgray" ) +
  geom_point( aes( x=CT, y=absBlue-absYellow ) ) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + 
  scale_x_continuous( breaks = c( 20, 30, 40, 44.5 ), labels = c( 20, 30, 40, "neg" ) ) +
  xlab( "CT value from qPCR" ) + 
  ylab( expression( OD["434 nm"] -  OD["560 nm"] ) ) +
  ggtitle( "40 min, gene 1a" )

dev.copy( svg, "figs/Figure_3d.svg", width=3.8, height=2.5 )
dev.off()

  