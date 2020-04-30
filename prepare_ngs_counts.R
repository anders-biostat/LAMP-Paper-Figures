# ##The count table is rather large (~1.6 GB). The folowing was done on my workstation to unify by UMI. This can be revisited later.
# counts_raw <- read_tsv("Counted/head100k-Undetermined_S0_R1_001_pass_counted.tsv") %>% #read_tsv("Counted/Undetermined_S0_R1_001_pass_counted.tsv") %>%
#   mutate(WellRow = str_sub(Well, 1, 1),
#          WellCol = as.integer(str_sub(Well, 2)))
# counts <- counts_raw %>%
#   rowwise() %>%
#   mutate(matched = sum(Seq1, Seq2, Seq3) > 0) %>%
#   ungroup() %>%
#   group_by(Plate, Well, WellRow, WellCol, matched) %>%
#   summarize(count = n()) %>%
#   group_by(Plate) %>%
#   mutate(cpm = count / sum(count) * 1e6) %>%
#   ungroup()
# save(counts, file = "counts.Rdata")

library(dplyr)
library(readr)

load( "data/counts.Rda" )

counts %>% 
  filter(!(Plate %in% sprintf("Plt%02d", c(4, 17:20)))) %>% # remove all plates which haven't been sequenced
  mutate( plate = str_replace( Plate, "Plt", "CP000" ) ) %>%
  rename(row = WellRow, col = WellCol) %>%
  select(plate, row, col, matched, cpm) %>%
  pivot_wider( names_from = "matched", values_from = "cpm", names_prefix = "matched" ) %>%
  write_tsv( "data/ngs_counts.tsv"  )