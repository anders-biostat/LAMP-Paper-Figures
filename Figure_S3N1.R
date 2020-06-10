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

panel_a <- rsvg::rsvg("SVGs/Figure_S3N1a.svg")

lamp_product <- data.frame(pos = 28515, start = 28515, stop = 28751)
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

panel_b <- ggplot(depth) +
  geom_area(aes(pos/1e3, cumratio * .8 * max(reads) / 1e6 ), fill = "#b3cde3", colour = "#b3cde3", alpha = .3) +
  geom_hline(yintercept = cumratio_limit * .8/ 1e6, color = "lightgrey", linetype = 2 ) +
  geom_rect(aes(xmin = start/1e3, xmax = stop/1e3, ymin = -.45, ymax = -.05), data = lamp_product, fill = "#fff2ae") +
  geom_area(aes(pos/1e3, reads/1e6), fill = "darkgray", colour = "black", alpha = .3) +
  facet_zoom(x = between(pos, 2.75e4, 3e4), zoom.size = 3) +
  labs(title = bquote(paste("mapped reads (", .((matched_total_frac)*100), "% of", ~10^6, " reads)")),
       x = "genomic position (kbp)",
       y = expression(paste("mapped reads / ", 10^6, ""))) +
  scale_y_continuous(limits = c(-.5, cumratio_limit/1e6+.2), breaks = c(0, 3, 6), expand = c(0,0),
                     sec.axis = sec_axis(~ . / (.8 * max(.) - .2) * 100, breaks = c(0, 50, 100), name = "cumulative mapped\nreads (%)")) +
  theme_light() + theme( text = element_text(family = 'Arial'),
                         panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                        panel.grid = element_line(colour = "grey92"), 
                        panel.border = element_rect(colour = "grey20", inherit.blank = TRUE))
panel_b

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

# this is computationally slightly extensive, that's why only do if needed
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

match_colors <- c("primer" = "#b3e2cd", "primer (r. c.)" = "#cbd5e8", "other" = "#cccccc")
tbl_primer <- tbl_kmers %>% group_by(match) %>%
  summarise(counts = sum(count)) %>%
  ungroup %>%
  mutate(perc =  counts / sum(counts) * 100, perc_total = perc * (1-matched_total_frac)) %>%
  mutate(primer = map_chr(match, function(.){ifelse(is.na(.), "other", primer$name[[.]])})) %>%
  mutate(group = case_when(
    str_detect(primer, "_rc") ~ "primer (r. c.)",
    primer == "other" ~ "other",
    TRUE ~ "primer"),
    group = fct_relevel(group, c("other", "primer (r. c.)", "primer"))) %>%
  mutate(name = map_chr(str_split(primer, "_"), `[[`, 1))

panel_c <- tbl_primer %>%
  arrange(desc(group), desc(perc_total)) %>%
  #       xx = x - .5 * perc_total) %>%
  mutate(label = if_else(perc_total > 1, name, "")) %>%
  ggplot(aes("a", perc, fill = group)) +
  geom_col(colour= "black") +
  geom_text(aes(label = label), position = position_stack(vjust = .5), angle = 90) +
  scale_fill_manual(name  = "k-mer match", values = match_colors, guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(name = "fraction unmapped reads (%)", sec.axis = sec_axis(~.*(1-matched_total_frac), name = "fraction all reads (%)")) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(title = bquote(paste("unmapped reads (", .((1-matched_total_frac)*100), "% of", ~10^6, " reads)"))) +
  theme(legend.position = "bottom") +
coord_flip()
panel_c

fig_layout <- '
A
A
B
B
B
B
C
C
'
wrap_elements(plot =  grid::rasterGrob(panel_a)) + 
  wrap_elements(plot = panel_b) + 
  panel_c +
  plot_layout(design = fig_layout) +
  plot_annotation(tag_levels = "a")


# Export figures
ggsave("SVGs/Figure_S3N1_tmp.svg", width=20, height=22, units="cm")
ggsave("Figure_S3N1_tmp.png", width=20, height=22, units="cm", dpi=300)
