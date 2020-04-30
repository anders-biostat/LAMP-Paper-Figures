library( tidyverse )

read_tsv( "data/tecan_values.tsv" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

tecan %>% 
left_join( tblCT ) %>%
filter( !is.na(CT) ) %>%
filter( plate %in% "CP00001", gene=="N", plateRemark=="2" ) %>%
ggplot +
  geom_line( aes( x=minutes, y=absBlue-absYellow, group=well, col=CT ) )

tecan %>% 
filter( plate %in% c( "CP00001", "CP00003" ), minutes == 30 ) %>%
left_join( tblCT ) %>%
mutate( CT = ifelse( CT>40, runif( n(), 43, 46 ), CT ) ) %>%
ggplot +
  geom_point( aes( x=CT, y=absBlue-absYellow, col=plate ) )

  