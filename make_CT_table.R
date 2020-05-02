# This script requires files stored on my laptop at ~/w/repos/covid_test/tidy/
# The output of the script is the file data/plates_with_CTs.tsv

setwd( "~/w/repos/covid_test/tidy/" )

list.files( "qpcr_results", "xlsx$" ) %>%
map_chr( ~ file.path( "qpcr_results/", . ) ) %>%
map_dfr( readxl::read_excel, col_types="text" ) %>%
filter( CODE1 != "vNCOV" ) %>%
filter( !str_detect( CODE1, "IC" ) ) %>%
select( orderID = AUFTRAGNR, labID = SCHEINID, datum = LABEINDAT, test = CODE1, CP=ERGEBNIST1 ) -> qpcr_pre

# Get results for E only
qpcr_pre %>%
filter( !is.na(CP) ) %>%
filter( !is.na(suppressWarnings(as.numeric(CP))) | CP == "negativ" ) %>%
filter( str_starts( test, "vNCOVE" ) ) %>%
select( -datum, -test ) %>%
group_by( labID ) %>%
filter( CP == first(CP) ) %>%
summarise( CP = first(CP) ) %>%
# Add special file
bind_rows( 
  read_csv( "qpcr_results/CP00004.csv" ) %>%
  mutate( CP = str_replace( CP, "neg", "negativ" ) )
) %>%
bind_rows( 
  read_csv( "qpcr_results/CP00003+5+6.csv" ),
  read_csv( "qpcr_results/CP00014.csv" ) %>% select( labID, CP )
) -> resultsE


# Read all plate content files
bind_rows(
  list.files( "plate_contents", "^CP.*xlsx$" ) %>%
    map_chr( ~ file.path( "plate_contents", . ) ) %>%
    map_dfr( readxl::read_excel, col_types="text" ) %>% 
    mutate_at( "remark", replace_na, ""  ),
  list.files( "plate_contents", "CP.....\\.csv" ) %>%
    set_names( str_remove( ., ".csv" ) ) %>%
    map_chr( ~ file.path( "plate_contents", . ) ) %>%
    map_dfr( read_csv, col_types=cols(.default="c") , .id="plate" ) %>%
    select( -CP ) ) %>% 
  mutate( col = sprintf( "%02d", as.integer(col) ) ) %>% 
mutate( labID = ifelse( str_detect( labID, "^\\d{5}$" ), 
  str_c( "V20", labID ), labID ) )  -> plateContentsTbl

# Add plate contents files from barcodes
read_csv( "plate_contents/CP00007.o.csv" ) %>%
add_column( plate = "CP00007" ) %>% 
mutate_at( "orderID", str_remove, "97$" ) %>%
mutate( col = sprintf( "%02d", as.integer(col) ) ) %>% 
left_join( qpcr_pre %>% select( orderID, labID ) %>% unique() ) %>%
#filter( is.na( labID) )
select( -orderID ) %>%
bind_rows( plateContentsTbl ) -> plateContentsTbl

# Add plate 20
bind_rows( 
  plateContentsTbl,
  read_csv( "plate_contents/CP00020.s.csv" ) %>% 
  add_column( plate = "CP00020" ) %>% 
  rename( col=column ) %>%
  mutate( col = sprintf( "%02d", as.integer(col) ) ) %>% 
  mutate( Beads = str_c( "Beads ", Beads ) ) %>%
  unite( remark, `remarks H20 for dil.`, Material, Beads, sep=" / " )
) -> plateContentsTbl

# Add controls and bead wells
bind_rows( 
  plateContentsTbl,
  read_csv( "plate_contents/CP00003+5+6.s.csv" ) %>% 
    mutate( col = sprintf( "%02d", as.integer(col) ) ) ) -> plateContentsTbl

# Replace any order IDs with lab IDs, then remove still missing labIDs
plateContentsTbl %>% 
rename( "labOrOrderID" = "labID" ) %>%
mutate_at( "labOrOrderID", str_remove, "(?<=799.{5})97$" ) %>%
left_join( qpcr_pre %>% select( orderID, labID ), by = c( "labOrOrderID" = "orderID" ) ) %>% 
mutate( labID = coalesce( labID, labOrOrderID ) ) %>%
select( -labOrOrderID ) %>%
filter( !is.na( labID ) ) %>%
unique() -> plateContentsTbl

plateContentsTbl %>% 
group_by( plate, row, col ) %>%
summarise( n = n() ) %>%
assertr::verify( n == 1 ) 

## END taken fromn load_data.R

## BEGIN taken from ZeroLysis.R

bind_rows(

  # Plate without number, of 20-04-13, now named CP10020 and CP10021:
  read_csv( "plate_contents/CP10020+CP10021.s.csv") %>%
  rename( labID = sample ) %>%
  mutate( col = sprintf( "%02d", as.integer(col) ) ),
  
  # Plate 21
  plateContentsTbl %>% filter( plate == "CP00021" ) %>% left_join( resultsE ),
  
  # Plate 22
  readxl::read_excel( "tests/200416_Layout_CP00022.xlsx", skip=4 ) %>%
  rename( row=Reihe, col=Spalte, barcode=Barcode, remark=Remark ) %>%
  add_column( plate = "CP00022" ) %>%
  mutate( col = sprintf( "%02d", as.integer(col) ) ),
    
  # Plate 24-26
  bind_rows(
    read_csv( "plate_contents/CP00024.s.csv", col_types="cccc" ) %>% 
      add_column( plate="CP00024" ),
    read_csv( "plate_contents/CP00025.s.csv", col_types="cccc" ) %>% 
      add_column( plate="CP00025" ) %>%
      mutate( barcode = str_c( "V20", barcode ) ),
    readxl::read_excel( "plate_contents/200419_Layout_CP00026.xlsx", skip=4, col_types = "text" ) %>%
      rename( CP = "...4" ) %>%
      rename( barcode = Barcode ) %>%
      rename( row = Reihe, col = Spalte ) %>%
      add_column( plate = "CP00026" ) ) %>%
  mutate( col = sprintf( "%02d", as.integer(col) ) ) %>%
  select( plate, everything() ),
  
) %>%
mutate_at( "CP", str_replace, "negative?", "Inf" ) %>%
mutate_at( "CP", ~ ifelse( .=="unknown", NA, . ) ) %>%
mutate_at( "CP", as.numeric ) %>%
select( plate, row, col, CP, everything() ) -> platesA  

# Plates 35 and 36
read_csv( "plate_contents/CP00035.s.csv" ) %>%
mutate( col = sprintf( "%02d", as.integer(col) ) ) %>%
mutate( index = case_when(
  row == "D" & col == "11" ~ "#AA",
  row == "D" & col == "12" ~ "#AB",
  TRUE ~ index ) ) %>%
filter( row <= "D" ) -> a

bind_rows(
  a,
  a %>% rowwise() %>% mutate( row = intToUtf8( utf8ToInt(row) + 4 ) ) ) %>%
left_join(
  read_csv( "other_data/CP00035_CPs.csv", skip=2 ) %>%
    mutate( index = str_extract( Name, "(?<=CP00035_)#\\d\\d?" ) ) %>%
    mutate_at( "Cp", replace_na, Inf ) %>%
    select( CP = Cp, index ) %>%
    filter( !is.na(index) ) %>%
    # Add Loan's Zettel
    bind_rows( tibble( 
      index = c( "#AA", "#AB" ),
      CP = c( 35.8, 35.9 ) ) ) ) %>%
add_column( plate = "CP00035" ) %>%
mutate( col = sprintf( "%02d", as.integer(col) ) ) -> CP35

read_csv( "plate_contents/CP00036.s.csv") %>%
add_column( plate = "CP00036" ) %>%
mutate( col = sprintf( "%02d", as.integer(col) ) ) %>%
rename( "barcode" = "Barcode" ) %>%
mutate_at( "CP", ~ ifelse( .=="neg", Inf, suppressWarnings(as.numeric(.)) ) ) %>%
filter( !is.na(barcode) ) -> CP36

read_csv( "plate_contents/CP00037.s.csv") %>%
add_column( plate = "CP00037" ) %>%
mutate( col = sprintf( "%02d", as.integer(col) ) ) %>%
mutate_at( "CP", ~ ifelse( .=="Neg", Inf, suppressWarnings(as.numeric(.)) ) ) %>%
filter( !is.na(barcode) ) -> CP37

bind_rows(
  plateContentsTbl,
  platesA,
  CP35,
  CP36,
  CP37 ) -> plates

plates %>% pull( plate ) %>% unique %>% sort

plates %>%
left_join( resultsE, by="labID" ) %>% 
mutate( CT = coalesce( CP.y, as.character(CP.x) ) ) %>%
mutate( barcode = coalesce( labID, barcode ) ) %>%
arrange( plate, row, col ) %>%
mutate( well = str_c( row, str_remove( col, "^0" ) ) ) %>%
select( plate, well, barcode, CT, wellRemark=remark ) %>%
mutate( CT = ifelse( str_starts( CT, "neg" ), Inf, suppressWarnings(as.numeric(CT)) ) ) -> plateCTs

plateCTs %>%
mutate( barcode_prefix = str_sub( barcode, 1, 3 ) ) %>%
mutate( sensitive_barcode = barcode_prefix %in% c( "799", "V20", "204" ) ) %>%
mutate( barcode = ifelse( sensitive_barcode, 
  str_c( barcode_prefix, "_R", sample( 100000, n() ) ),
  barcode ) ) %>%
select( -barcode_prefix, -sensitive_barcode )
  
setwd( "~/w/repos/LAMP-Paper-Figures/" )
plateCTs %>% write_tsv( "data/plates_with_CTs.tsv" )
