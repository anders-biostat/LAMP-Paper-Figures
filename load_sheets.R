library( tidyverse )
library( tecan )

read_csv( "data/tecan_sheets.csv" ) -> sheet_list

1:nrow(sheet_list) %>%
map_dfr( function(i) read_tecan( 
  file.path( "data/tecan_scans", sheet_list$filename[i] ), 
  sheet_list$sheetname[i] ) ) %>%
left_join( sheet_list ) -> a
a %>%
mutate_at( "remark", replace_na, "" ) %>% 
mutate( heat95 = heat==95 ) %>%
select( plate, gene, heat95, minutes, plateRemark=remark, well, 
  absBlue=blue_absorbance, absYellow = yellow_absorbance ) -> tecan

tecan %>% write_tsv( "data/tecan_values.tsv" ) 
