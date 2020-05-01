library( tidyverse )
library( patchwork )

read_tsv( "data/tecan_values.tsv", col_types = "ccliccdd" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00035", "CP00036", "CP00037" ) 

tecan %>% 
filter( plate %in% plates_to_use ) %>%
left_join( tblCT ) %>%
mutate( CT = ifelse( CT>40, runif( n(), 43, 46 ), CT ) ) %>%
filter( !is.na(CT) ) %>%
ggplot + 
  geom_line( aes( x=minutes, y=absBlue-absYellow, group=str_c(plate,well), col=CT ) ) +
  facet_grid( . ~ heat95 )

tecan %>% 
filter( plate %in% plates_to_use, well=="A1", minutes==30 )
