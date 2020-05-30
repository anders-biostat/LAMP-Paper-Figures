library( tidyverse )
library( tecan )

read_csv( "data/tecan_sheets.csv" ) -> sheet_list

1:nrow(sheet_list) %>%
map_dfr( function(i) read_tecan( 
  file.path( "data/tecan_scans", sheet_list$filename[i] ), 
  sheet_list$sheetname[i] ) ) -> a
a %>% 
left_join( sheet_list ) %>%
mutate_at( "remark", replace_na, "(none)" ) %>% 
mutate_at( "exclude", replace_na, "no" ) %>% 
mutate( heat95 = heat==95 ) %>%
select( plate, gene, heat95, minutes, extraction, plateRemark=remark, exclude, well, 
  absBlue=blue_absorbance, absYellow = yellow_absorbance ) -> tecan

tecan %>%
group_by( plate, minutes, heat95, extraction, plateRemark, exclude, well ) %>%
count() %>%
group_by( plate, minutes, heat95, extraction, plateRemark, exclude ) %>%
summarise( n=max(n) ) %>%
assertr::verify( n==1 ) 

tecan %>% write_tsv( "data/tecan_values.tsv" ) 

