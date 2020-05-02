library( tidyverse )
library( patchwork )
library( ggtext )
library( glue )
library( ggsci )
library( binom )
source( "misc.R" )

# Load data
read_tsv( "data/tecan_values.tsv", col_types = "cclicccdd" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_w_repl <- paste0("CP000", c(35, 36, 37))

deltaOD_cutoff_zeroL <- 0.3

tecan %>%
  filter( plate %in% plates_w_repl ) %>%
  left_join( tblCT, by = c("plate", "well") ) %>% 
  filter( !is.na(CT) & minutes == 30 & heat95 == TRUE ) -> tbl

duplicated_samples <- tbl %>%
  group_by(plate, barcode) %>%
  tally() %>% ungroup() %>%
  filter(n == 2) %>%
  pull(barcode)

duplicated_samples

## No duplicates found :(


zero_dupl <- lamp_df %>%
  filter(minutes == 30 & celsius == 95 & plate %in% plates_w_repl & barcode %in% duplicated_samples) %>%
  mutate(delta_od = blue_absorbance - yellow_absorbance) %>%
  dplyr::select(-blue_absorbance, -yellow_absorbance, -well, -row, -col, -remark, -index, -run) %>%
  arrange(plate, barcode) %>%
  group_by(plate, barcode) %>%
  mutate(replicate_n = paste0("replicate_", row_number())) %>% ungroup() %>%
  pivot_wider(names_from = replicate_n, values_from = delta_od) %>%
  mutate(CP_jitter = ifelse( CP > 40, runif( n(), 44, 48 ), CP ),
         lamp_max = if_else(replicate_1 > replicate_2, replicate_1, replicate_2))

zero_dupl %>%
  mutate(CP = if_else(CP == Inf, 45, CP)) %>%
  ggplot(aes(x = replicate_1, y = replicate_2, color = CP)) +
  geom_point() +
  scale_color_gradientn(
    guide = guide_colourbar(reverse = TRUE),
    colors = cp_colors,
    breaks = rev(c(20, 25, 30, 35, 40, 43.5)), labels = rev(c(20, 25, 30, 35, 40, "neg"))) +
  labs(x = "LAMP ΔOD of 1st replicate",
       y = "LAMP ΔOD of 2nd replicate",
       color = "CT-value",
       title = "30min 95°C zero lysis") +
  theme_bw() +
  coord_fixed()