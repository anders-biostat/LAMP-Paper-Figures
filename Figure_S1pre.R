library( tidyverse )
source( "~/w/repos/LAMP-Paper-Figures/misc.R" )

readxl::read_excel( "~/w/repos/covid_test/tidy/plate_contents/CP00001+2.xlsx" ) %>%
  filter( plate == "CP00001" ) -> plate_content
read_csv( "~/w/repos/covid_test/tidy/qpcr_results/alles_CT_bis20200402.csv", 
  col_types = cols(.default = "c") ) -> qPCRres

plate_content %>%
left_join( qPCRres, by=c( "labID"="SCHEINID" ) ) %>%
filter( CODE1 %in% c( "vNCOVE", "vNCOVN" ) ) %>%
select( plate, row, col, labID, CODE1, ERGEBNIST1 ) %>%
pivot_wider( names_from = CODE1, values_from = ERGEBNIST1 ) %>%
unite( well, row, col, sep="" ) -> qPCRpl1

read_tsv( "~/w/repos/LAMP-Paper-Figures/data/tecan_values.tsv" ) %>%
  filter( plate == "CP00001", gene == "N", minutes==30 ) %>%
left_join( qPCRpl1 ) -> tbl

tbl %>%
mutate( resE = ifelse( vNCOVE=="negativ", 44 , 
  suppressWarnings(as.numeric(vNCOVE)) ) ) %>%
mutate( rnd = runif( n(), 43, 47 ) ) %>%
mutate_at( c( "vNCOVE", "vNCOVN" ), 
   ~ ifelse( .=="negativ", rnd, suppressWarnings(as.numeric(.)) ) ) %>%  
mutate( resE = ifelse( is.na(vNCOVN), NA, resE ) ) %>%
ggplot() +
  geom_vline( xintercept = 42, col="darkgray" ) +
  geom_point( aes( x = vNCOVE, fill=resE, y = absBlue - absYellow ), 
    colour = "black", alpha = .8, shape = 21, size = 1.6 ) + 
  scale_x_reverse( breaks = c( 10, 20, 30, 40, 45 ), labels = c( 10, 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value for E gene)" ) + 
  ylab( "RT-LAMP (ΔOD)" ) +
  scale_fill_ct( limits=c(12,47),  name="CT_E" )

tbl %>%
mutate( resE = ifelse( vNCOVE=="negativ", 44 , 
  suppressWarnings(as.numeric(vNCOVE)) ) ) %>%
mutate( rnd = runif( n(), 43, 47 ) ) %>%
mutate_at( c( "vNCOVE", "vNCOVN" ), 
   ~ ifelse( .=="negativ", rnd, suppressWarnings(as.numeric(.)) ) ) %>%  
mutate( resE = ifelse( is.na(vNCOVN), NA, resE ) ) %>%
ggplot() +
  geom_vline( xintercept = 42, col="darkgray" ) +
  geom_point( aes( x = vNCOVE, fill=resE, y = absBlue - absYellow ), 
    colour = "black", alpha = .8, shape = 21, size = 1.6 ) + 
  scale_x_reverse( breaks = c( 10, 20, 30, 40, 45 ), labels = c( 10, 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value for N gene)" ) + 
  ylab( "RT-LAMP (ΔOD)" ) +
  scale_fill_ct( limits=c(12,47),  name="CT_E" )

tbl %>%
mutate_at( c( "vNCOVE", "vNCOVN" ), as.numeric ) %>%
ggplot() + 
  geom_point( aes(x=vNCOVE,y=vNCOVN,fill=vNCOVE), 
    colour = "black", alpha = .8, shape = 21, size = 1.6 ) +
  coord_fixed() +
  geom_abline() +
  scale_fill_ct( limits=c(12,47),  name="CT_E" ) +
  xlab("RT-qPCR (CT value for E gene)" ) +
  ylab("RT-qPCR (CT value for N gene)" ) 




#####

tbl %>%
mutate_at( c( "vNCOVE", "vNCOVN" ), 
   ~ ifelse( .=="negativ", runif( n(), 43, 47 ), suppressWarnings(as.numeric(.)) ) ) %>%  
ggplot() +
  geom_vline( xintercept = 42, col="darkgray" ) +
  geom_point( aes( x = vNCOVE, y = absBlue - absYellow ), 
    fill = "black", colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
  scale_x_reverse( breaks = c( 10, 20, 30, 40, 45 ), labels = c( 10, 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value)" ) + 
  ylab( "RT-LAMP (ΔOD)" ) + theme_bw() + theme(panel.grid = element_blank())

tbl %>%
mutate_at( c( "vNCOVE", "vNCOVN" ), 
   ~ ifelse( .=="negativ" | is.na(.), runif( n(), 43, 47 ), suppressWarnings(as.numeric(.)) ) ) %>%  
ggplot() +
  geom_vline( xintercept = 42, col="darkgray" ) +
  geom_point( aes( x = vNCOVN, y = absBlue - absYellow ), 
    fill = "black", colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
  scale_x_reverse( breaks = c( 10, 20, 30, 40, 45 ), labels = c( 10, 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value for gene 'N')" ) + 
  ylab( "RT-LAMP (ΔOD)" ) + theme_bw() + theme(panel.grid = element_blank())
