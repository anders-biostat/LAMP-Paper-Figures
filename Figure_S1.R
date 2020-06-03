library( tidyverse )
source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "ccliccccdd" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT
read_tsv( "data/qPCR_CP00001_E_and_N.tsv" ) -> tblCT_extra

tblCT %>% 
filter( plate == "CP00001" ) %>%
left_join( tblCT_extra ) %>%
mutate_at( c( "vNCOVE", "vNCOVN" ), 
   ~ ifelse( .=="negativ", Inf, suppressWarnings(as.numeric(.)) ) ) %>%  
assertr::verify( ( !is.finite(CT) & !is.finite(vNCOVE) ) | abs( CT - vNCOVE ) < .2 ) ->
qPCRpl1  

read_tsv( "data/tecan_values.tsv" ) %>%
  filter( plate == "CP00001", gene == "N", minutes==30 ) %>%
left_join( qPCRpl1 ) -> tbl

tbl %>%
mutate_at( c( "vNCOVE", "vNCOVN" ), 
   ~ ifelse( .=="negativ", runif( n(), 43, 47 ), suppressWarnings(as.numeric(.)) ) ) %>%  
ggplot() +
  geom_vline( xintercept = 42, col="lightgrey" ) +
  geom_point( aes( x = vNCOVE, y = absBlue - absYellow ), 
    fill = "black", colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
  scale_x_reverse( breaks = c( 10, 20, 30, 40, 45 ), labels = c( 10, 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value)" ) + 
  ylab( "RT-LAMP (ΔOD)" ) + theme_bw() + theme(panel.grid = element_blank())

tbl %>%
mutate_at( c( "vNCOVE", "vNCOVN" ), 
   ~ ifelse( .=="negativ" | is.na(.), runif( n(), 43, 47 ), suppressWarnings(as.numeric(.)) ) ) %>%  
ggplot() +
  geom_vline( xintercept = 42, col="lightgray" ) +
  geom_point( aes( x = vNCOVN, y = absBlue - absYellow ), 
    fill = "black", colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
  scale_x_reverse( breaks = c( 10, 20, 30, 40, 45 ), labels = c( 10, 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value for gene 'N')" ) + 
  ylab( "RT-LAMP (ΔOD)" ) + theme_bw() + theme(panel.grid = element_blank())

