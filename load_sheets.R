library( tidyverse )
library( tecan )

read_tsv( "data/tecan_sheets.tsv" ) -> sheet_list

1:nrow(sheet_list) %>%
map_dfr( function(i) read_tecan( 
  file.path( "data/tecan_scans", sheet_list$filename[i] ), 
  sheet_list$sheetname[i] ) ) %>%
left_join( sheet_list ) %>%
select( plate, gene, heat, minutes, plateRemark=remark, well, 
  absBlue=blue_absorbance, absYellow = yellow_absorbance ) -> tecan

tecan %>% write_tsv( "data/tecan_values.tsv" ) 
