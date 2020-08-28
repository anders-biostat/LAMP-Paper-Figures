library(dplyr)
library(tidyr)
library(readr)
library(scales)
library(forcats)
library(purrr)
library(stringr)
library(ggplot2)
library(patchwork)
library(ggforce)

source("misc.R")

panel_a <- rsvg::rsvg("SVGs/Figure_S4a.svg")
panel_b <- rsvg::rsvg("SVGs/Figure_S4b.svg")

lamp_product <- data.frame(pos = 28515, start = 28515, stop = 28751)
N_gene <- data.frame(pos = 28274, start = 28274, stop = 29533)
depth <-
  read_tsv("data/ngs_SARS-CoV-2_10Mreads_mapped_depth.tsv",
                  col_names = c("ref", "pos", "reads")) %>%
  mutate(pos_bin = cut(pos, seq(1, nrow(.), 5))) %>%   ## bin over 5 consecutive bp
  filter(!is.na(pos_bin)) %>%
  group_by(pos_bin) %>%
  summarise(reads = mean(reads)) %>%
  ungroup() %>%
  mutate(pos = as.numeric(str_match(pos_bin, "\\((.+),")[,2])) %>%
  mutate(cumratio = cumsum(reads) / sum(reads)) %>%
  mutate(zoom = between(pos, 2.75e4, 3e4))
cumratio_limit <- max(depth$reads)

panel_c <- ggplot(depth) +
  geom_hline(yintercept = cumratio_limit * .8/ 1e6, color = "lightgrey", linetype = 2 ) +
  geom_area(aes(pos/1e3, cumratio * .8 * max(reads) / 1e6 ), fill = "#b3cde3", colour = "#b3cde3", alpha = .3) +
  geom_rect(aes(xmin = start/1e3, xmax = stop/1e3, ymin = -.65, ymax = -.05), data = N_gene, fill = "#888888") +
  geom_rect(aes(xmin = start/1e3, xmax = stop/1e3, ymin = -.65, ymax = -.05), data = lamp_product, fill = "#fff2ae") +
  geom_text(aes(x = (start/1e3 + stop/1e3)/2, y = -.3), label = "N gene", data = N_gene, colour = "white") +
  #annotate("text", x = 29, y = -.3, label = "N gene", colour = "lightgrey") +
  geom_area(aes(pos/1e3, reads/1e6), fill = "darkgray", colour = "black", alpha = .3) +
  facet_zoom(x = between(pos, 2.75e4, 3e4), zoom.size = 3) +
  labs(title = bquote(paste("mapped reads (", .((matched_total_frac)*100), "%)")),
       x = "SARS-CoV-2 genome position (kbp)",
       y = expression(paste("mapped reads / ", 10^6, ""))) +
  scale_y_continuous(limits = c(-.7, cumratio_limit/1e6+.2), breaks = c(0, 3, 6), expand = c(0,0),
                     sec.axis = sec_axis(~ . / (.8 * max(.) - .2) * 100, breaks = c(0, 50, 100), name = "cumulative mapped\nreads (%)")) +
  theme_light() + theme( text = element_text(family = 'Arial'), plot.title = element_text(hjust = .5),
                         panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                        panel.grid = element_line(colour = "grey92"), 
                        panel.border = element_rect(colour = "grey20", inherit.blank = TRUE))
panel_c

primer <- read_tsv("data/LAMP-primer-gene_N.tsv")
primer <- bind_rows(
  primer,
  primer %>% mutate(name = str_c(name, "_rc"), sequence = rc(sequence))
)

file_kmers <- "data/ngs_SARS-CoV-2_10Mreads_unmapped_9kmers.tsv"
kmers <- read_tsv(file_kmers) %>%
  rename(kmer = Kmer, count = Count) 

match_primer <- function(kmer){
  res <- which(Biostrings::vcountPattern(kmer, primer$sequence, max.mismatch = 2, with.indels = TRUE) > 0)
  ifelse(length(res) > 0, res, integer(0))
}

# this is computationally slightly extensive, that's why only do it if needed
file_kmers_primers <- "ngs_SARS-CoV-2_10Mreads_unmapped_9kmers_primers.tsv"
if(!file.exists(file_kmers_primers) & (file.info(file_kmers)[["mtime"]] > file.info(file_kmers_primers)[["mtime"]])){
tbl_kmers <- kmers %>%
  mutate(rank = row_number(), cumsum = cumsum(count), cumratio = cumsum/max(cumsum)) %>%
  mutate(match = map_int(kmer, match_primer)) %>%
  write_tsv(file_kmers_primers)
} else {
  tbl_kmers <- read_tsv(file_kmers_primers)
}

matched_total_frac <- 0.806

#match_colors <- c("primer" = "#b3e2cd", "primer (r. c.)" = "#cbd5e8", "no match" = "#cccccc")
match_colors <- c("primer" = "#fdbf6f", "no match" = "#cccccc")
tbl_primer <- tbl_kmers %>% group_by(match) %>%
  summarise(counts = sum(count)) %>%
  ungroup %>%
  mutate(perc =  counts / sum(counts) * 100, perc_total = perc * (1-matched_total_frac)) %>%
  mutate(primer = map_chr(match, function(.){ifelse(is.na(.), "no match", primer$name[[.]])})) %>%
  mutate(group = case_when(
    str_detect(primer, "_rc") ~ "primer (r. c.)",
    primer == "no match" ~ "no match",
    TRUE ~ "primer"),
    group = fct_relevel(group, c("no match", "primer (r. c.)", "primer"))) %>%
  mutate(name = map_chr(str_split(primer, "_"), `[[`, 1))

unmatched_primer_frac <- round(sum(filter(tbl_primer, group != "no match")$perc_total)/(1-matched_total_frac), 1)
panel_d <- tbl_primer %>%
  group_by(name) %>%
  summarise(perc = sum(perc), perc_total = sum(perc_total)) %>%
  ungroup() %>%
  mutate(group = if_else(name == "no match", "no match", "primer")) %>%
  arrange(desc(group), desc(perc_total)) %>%
  mutate(label = if_else(perc_total > .5, str_remove(name, "GeneN-A-"), "")) %>%
  ggplot(aes("a", perc, fill = group)) +
  geom_col(colour= "black") +
  geom_hline(yintercept = unmatched_primer_frac ) +
  geom_text(aes(label = label), position = position_stack(vjust = .5), size = 3) +
  scale_fill_manual(name  = "k-mer match:", values = match_colors, guide = "none") +# guide_legend(reverse = TRUE)) +
  scale_y_continuous(name = "fraction of unmapped reads (%)",
                     breaks = c(0,25,50,75,100,unmatched_primer_frac),
                     sec.axis = sec_axis(~.*(1-matched_total_frac), name = "fraction of all reads (%)",
                                         breaks = c(0,5,10,15,20, round(unmatched_primer_frac*(1- matched_total_frac), 1) ))) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(title = bquote(paste("unmapped reads (", .((1-matched_total_frac)*100), "%)"))) +
  theme(legend.position = "bottom") +
coord_flip() +
  theme(plot.title = element_text(hjust = .5), panel.border = element_blank(), axis.line.x = element_line(colour = "grey20"))
panel_d

wrap_elements(plot =  grid::rasterGrob(panel_a)) +
  wrap_elements(plot =  grid::rasterGrob(panel_b)) +
  wrap_elements(plot = panel_c) +
  panel_d +
  plot_layout(nrow = 4, heights = c(4,7,9,1)) +
  plot_annotation(tag_levels = "A")


# Export figures
ggsave("SVGs/Figure_S4_tmp.svg", width=20, height=27, units="cm")
ggsave("Figure_S4_tmp.png", width=20, height=27, units="cm", dpi=300)
