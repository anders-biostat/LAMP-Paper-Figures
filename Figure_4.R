library( tidyverse )

read_tsv( "data/tecan_values.tsv" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00001", "CP00003", "CP00005", "CP00006", "CP00008", 
   "CP00009", "CP00010", "CP00011", "CP00012", "CP00013", "CP00016" )

thresh_lower <- -.15  
thresh_upper <- .25

## Figure 4a
tecan %>% 
filter( plate %in% plates_to_use ) %>%
filter( minutes==30, gene=="N" ) %>% 
filter( ! ( plate == "CP00005" & well >= "G" ) ) %>%
left_join( tblCT ) %>% 
filter( !is.na(CT) ) -> tbl

tbl %>%
mutate( CT = ifelse( CT>40, runif( n(), 43, 46 ), CT ) ) %>%
ggplot +
  geom_vline( xintercept = 41.5, col="darkgray" ) +
  geom_vline( xintercept = 30, col="lightgray" ) +
  geom_hline( yintercept = .25, col="lightgray" ) +
  geom_hline( yintercept = -.15, col="lightgray" ) +
  geom_point( aes( x=CT, y=absBlue-absYellow, col=plate ) ) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + 
  scale_x_continuous( breaks = c( 20, 30, 40, 44.5 ), labels = c( 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value)" ) + 
  ylab( "LAMP (ΔOD)" ) +
  ggtitle( "", subtitle = str_interp( "${nrow(tbl)} samples on ${tbl%>%select(plate)%>%unique%>%nrow} plates" )  )

tbl %>%
mutate( CTbin = cut( CT, c( 0, 25, 30, 35, 40, Inf ) ) ) %>% 
mutate_at( "CTbin", recode, 	"neg" = "(40,Inf]" ) %>%
mutate( LAMPres = 
   cut( absBlue-absYellow, c( -Inf, thresh_lower, thresh_upper, Inf ) ) %>%
   as.integer %>%
   { c( "neg", "incl", "pos" )[.] } ) %>%
group_by( CTbin, LAMPres ) %>%
count() %>%
pivot_wider(

dev.copy( svg, "figs/Figure_4a.svg", width=3.8, height=2.5 )
dev.off()
  
