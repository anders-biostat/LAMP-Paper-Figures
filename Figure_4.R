library( tidyverse )

read_tsv( "data/tecan_values.tsv" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00001", "CP00002", "CP00003", "CP00005", "CP00006", "CP00008", "CP00009",
   "CP00010", "CP00011", "CP00012", "CP00013", "CP00016" )
  
## Figure 4a
tecan %>% 
filter( plate %in% plates_to_use ) %>%
filter( minutes==30, gene=="N" ) %>% 
left_join( tblCT ) %>% 
mutate( CT = ifelse( CT>40, runif( n(), 43, 46 ), CT ) ) %>%
filter( !is.na(CT) ) %>%
ggplot +
  geom_vline( xintercept = 41.5, col="darkgray" ) +
  geom_vline( xintercept = 30, col="lightgray" ) +
  geom_hline( yintercept = .25, col="lightgray" ) +
  geom_hline( yintercept = -.15, col="lightgray" ) +
  geom_point( aes( x=CT, y=absBlue-absYellow, col=plate ) ) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + 
  scale_x_continuous( breaks = c( 20, 30, 40, 44.5 ), labels = c( 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value)" ) + 
  ylab( "LAMP (ΔOD)" )  

dev.copy( svg, "figs/Figure_4a.svg", width=3.8, height=2.5 )
dev.off()
  
