library(readr)
library(dplyr)
library(stringr)

count_file <- tail(commandArgs(), 1)

if(!file.exists(count_file)) stop(str_interp("File '${count_file}' does not exist. Exiting..."))

counts_raw <- read_tsv(count_file) %>%
    mutate(WellRow = str_sub(Well, 1, 1),
           WellCol = as.integer(str_sub(Well, 2)))
counts <- counts_raw %>%
    rowwise() %>%
    mutate(matched = sum(Seq1, Seq2, Seq3) > 0) %>%
    ungroup() %>%
    group_by(Plate, Well, WellRow, WellCol, matched) %>%
    summarize(count = n()) %>%
    ungroup() %>%
    write_tsv("counts.tsv")
save(counts, file = "counts.Rda")
