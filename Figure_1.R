library( tidyverse )

read_tsv( "data/tecan_values.tsv" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT


## Figure 1a 

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
  ylab( expression( OD["434 nm"] -  OD["560 nm"] ) )


## Figure 1b

tecan %>% 
filter( plate %in% c( "CP00001" ), gene=="N", plateRemark=="2", minutes == 30 ) %>%
left_join( tblCT ) %>%
mutate( CT = ifelse( CT>40, runif( n(), 43, 46 ), CT ) ) %>%
ggplot +
  geom_vline( xintercept = 41.5, col="darkgray" ) +
  geom_point( aes( x=CT, y=absBlue-absYellow ) ) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + 
  scale_x_continuous( breaks = c( 20, 30, 40, 44.5 ), labels = c( 20, 30, 40, "neg" ) ) +
  ggtitle( "30 min, gene N" )
  

## Figure 1c

tecan %>% 
filter( plate %in% c( "CP00001" ), gene=="1a", plateRemark=="2", minutes == 30 ) %>%
left_join( tblCT ) %>%
mutate( CT = ifelse( CT>40, runif( n(), 43, 46 ), CT ) ) %>%
ggplot +
  geom_vline( xintercept = 41.5, col="darkgray" ) +
  geom_point( aes( x=CT, y=absBlue-absYellow ) ) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + 
  scale_x_continuous( breaks = c( 20, 30, 40, 44.5 ), labels = c( 20, 30, 40, "neg" ) ) +
  ggtitle( "30 min, gene 1a" )

  