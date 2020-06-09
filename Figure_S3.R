library(dplyr)
library(tidyr)
library(readr)
library(scales)
library(forcats)
library(stringr)
library(ggplot2)
library(patchwork)
library(ggsci)
library(ggforce)

source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "ccliccccdd" ) -> tecan
read_tsv( "data/ngs_counts.tsv" ) -> ngs
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00003", "CP00005", "CP00006", "CP00008", 
                    "CP00009", "CP00010", "CP00011", "CP00012", "CP00013", "CP00016" )

lamp_thresholds <- c(.3)
qpcr_thresholds <- c(30, 42)
ngs_threshold <- c(100, 1e4)

tbl <- ngs %>%
  full_join( tecan ) %>%
  filter( plate %in% plates_to_use ) %>%
  filter( !(plate == "CP00003" & plateRemark != "2")) %>%
  filter( minutes == 30, gene=="N" ) %>% 
  filter( ! ( plate == "CP00005" & well >= "G" ) ) %>%
  left_join( tblCT )%>%
  filter( !is.na(CT) ) %>%
  filter( is.na(wellRemark) )

lamp_product <- data.frame(pos = 28515, start = 28515, stop = 28751)
depth <-
  read_tsv("~/Desktop/LAMP/NGS/TMP/Alignments/SARS-CoV2-ASM985889v3_Undetermined_S0_R1_001_pass_trimmed_10M_sorted_depth.tsv",
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

panel_a <- ggplot(depth) +
  geom_area(aes(pos/1e3, cumratio * .8 * max(reads) / 1e6 ), fill = "#b3cde3", colour = "#b3cde3", alpha = .3) +
  geom_hline(yintercept = cumratio_limit * .8/ 1e6, color = "lightgrey", linetype = 2 ) +
  geom_rect(aes(xmin = start/1e3, xmax = stop/1e3, ymin = -.45, ymax = -.05), data = lamp_product, fill = "#fed9a6") +
  geom_area(aes(pos/1e3, reads/1e6), fill = "darkgray", colour = "black", alpha = .3) +
  facet_zoom(x = between(pos, 2.75e4, 3e4), zoom.size = 3) +
  labs(title = str_interp("virus-mapped reads (${(matched_total_frac)*100}% of 10 M reads)"),
       x = "genomic position (kbp)",
       y = expression(paste("mapped reads / ", 10^6, ""))) +
  scale_y_continuous(limits = c(-.5, cumratio_limit/1e6+.2), breaks = c(0, 3, 6), expand = c(0,0),
                     sec.axis = sec_axis(~ . / (.8 * max(.) - .2) * 100, breaks = c(0, 50, 100), name = "cumulative mapped\nreads (%)")) +
  theme_light() + theme( text = element_text(family = 'Arial'),
                         panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                        panel.grid = element_line(colour = "grey92"), 
                        panel.border = element_rect(colour = "grey20", inherit.blank = TRUE))
panel_a

panel_b <- wrap_elements(panel = grid::textGrob('Here there will be plot showing that\n18.1% of reads contain primer sequences.'))

panel_c <- wrap_elements(panel = grid::textGrob('Here there will be the N gene\nsequence with primer annotation.'))

set.seed(42)
panels_data <-  tbl %>%
  replace_na(list(matchedTRUE = 0, matchedFALSE = 0)) %>%
  mutate( NGS = case_when(
    matchedTRUE > ngs_threshold[[2]] ~ "positive",
    between(matchedTRUE, ngs_threshold[[1]], ngs_threshold[[2]]) ~ "negative",
    matchedTRUE <= ngs_threshold[[1]]           ~ "undetected")) %>%
  mutate(CT = ifelse( CT>40, runif( n(), 43, 47 ), CT ) ) %>%
  mutate(plate = str_extract(plate, "\\d{2}$")) %>%
  mutate(well = fct_reorder(well, -CT))

mybreaks <- c( 0.1, 10, 100, 1000, 1e4, 1e5, 1e6 )
lamp_colors <- c("positive" = "#4daf4a", "negative" = "#ff7f00", "undetected" = "black")
panel_d <- panels_data %>%
  mutate(NGSshape = if_else(NGS != "undetected", "detected", NGS)) %>%
  ggplot(aes(x = matchedTRUE, y = matchedTRUE + matchedFALSE)) +
  geom_vline(xintercept = ngs_threshold[[1]], color = "lightgray", linetype = 2 ) +
  geom_vline(xintercept = ngs_threshold[[2]], color = "lightgray" ) +
  geom_point(aes( fill = NGS, shape = NGSshape ), colour = "black", alpha = .6, size = 1.2 ) +
  scale_fill_manual(name  = "LAMP-sequencing", values = lamp_colors) +
  scale_shape_manual(values = c("detected" = 21, "undetected" = 4)) +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 21, 4)), title.hjust = 0.5),
         shape = "none") +
  scale_x_continuous( trans=logap_trans(), breaks = mybreaks ) +
  scale_y_continuous( trans=logap_trans(), breaks = mybreaks ) +
  labs(x = "virus-matching read count",
       y = "total read count") +
  theme(legend.position = c(0.8, 0.2), legend.background = element_rect(fill = "white"), legend.box.background = element_blank(), legend.key=element_blank())
panel_d

# x <- tbl %>%
#   mutate(matchedTRUE = if_else(is.na(matchedTRUE), .1, matchedTRUE), matchedFALSE = if_else(is.na(matchedFALSE), .1, matchedFALSE)) %>%
#   mutate( NGS = case_when(
#     matchedTRUE > ngs_threshold ~ "positive",
#     matchedTRUE <= ngs_threshold ~ "negative",
#     TRUE                ~ "undetected")) %>%
#   filter(minutes == 30)
# ggplot(x) +
#   geom_vline(xintercept = ngs_threshold/1e4, color = "lightgray" ) +
#   geom_point( aes( x = matchedTRUE/1e4, y = matchedFALSE/1e4, fill = NGS),
#               colour = "black", alpha = .6, shape = 21, size = 1.2 ) +
#   ggrepel::geom_text_repel(aes(x = matchedTRUE/1e4, y = matchedFALSE/1e4, label = well, color = (matchedTRUE==.1|matchedFALSE==.1)),
#                            data = filter(x, plate == "CP00012")) +
#   guides(color = FALSE) +
#   scale_x_log10() + scale_y_log10() +
#   scale_fill_manual(name  = "LAMP-sequencing", values = lamp_colors) +
#   labs(x = expression(atop("virus-matching seq.", paste("(CPM / ", 10^4, ")"))),
#        y = expression(atop("unmatched seq.", paste("(CPM / ", 10^4, ")")))) #+
# 

panel_e <-  panels_data %>%
  ggplot(aes(x = matchedTRUE, y = absBlue - absYellow )) +
  geom_hline(yintercept = lamp_thresholds, color = "lightgray" ) +
  geom_vline(xintercept = ngs_threshold[[1]], color = "lightgray", linetype = 2 ) +
  geom_vline(xintercept = ngs_threshold[[2]], color = "lightgray" ) +
  geom_point(aes( fill = CT ),
             data = filter(panels_data, NGS != "undetected"),
              colour = "black", shape = 21, alpha = .6, size = 1.2 ) +
  geom_point(data = filter(panels_data, NGS == "undetected"), shape = 4,
             colour = "black", fill = lamp_colors["undetected"], alpha = .6, size = 1.2) +
  scale_x_continuous(trans=log1p_trans(), breaks = mybreaks ) +
  labs(subtitle = str_interp( "${unique(tbl$minutes)} min at 65°C\n${length(tbl$minutes)} samples on ${length(unique(tbl$plate))} plates"), 
       x = "LAMP-sequencing\n(read count in sample)",
       y = expression(paste("RT-LAMP (", Delta, OD, ")"))) +
  scale_fill_ct( name="RT-qPCR\nCT value") +
  annotate("text", x = 0.1, y = 0, label = str_glue("negative"), angle = 90, col="grey50") +
  annotate("text", x = 0.1, y = .425, label = str_glue("positive"), angle = 90, col="grey50")
panel_e

panel_f <-  panels_data %>%
  ggplot(aes(x = matchedTRUE, y = absBlue - absYellow )) +
  geom_hline(yintercept = lamp_thresholds, color = "lightgray" ) +
  geom_vline(xintercept = ngs_threshold[[1]], color = "lightgray", linetype = 2 ) +
  geom_vline(xintercept = ngs_threshold[[2]], color = "lightgray" ) +
  geom_point(aes( fill = plate ),
             data = filter(panels_data, NGS != "undetected"),
             colour = "black", shape = 21, alpha = .6, size = 1.2 ) +
  geom_point(data = filter(panels_data, NGS == "undetected"), shape = 4,
             colour = "black", fill = lamp_colors["undetected"], alpha = .6, size = 1.2) +
  scale_x_continuous(trans=log1p_trans(), breaks = mybreaks ) +
  labs(subtitle = str_interp( "${unique(tbl$minutes)} min at 65°C\n${length(tbl$minutes)} samples on ${length(unique(tbl$plate))} plates"), 
       x = "LAMP-sequencing\n(read count in sample)",
       y = expression(paste("RT-LAMP (", Delta, OD, ")"))) +
  scale_fill_d3(palette="category20", guide = guide_legend(label = FALSE)) +
  annotate("text", x = 0.1, y = 0, label = str_glue("negative"), angle = 90, col="grey50") +
  annotate("text", x = 0.1, y = .425, label = str_glue("positive"), angle = 90, col="grey50")
panel_f


fig_layout <- '
AABB
AABB
CCDD
CCDD
EEFF
EEFF
'
wrap_elements(full = panel_a + theme(plot.margin = margin(t = 15, r = 5.5, b = 5.5, l = 15, "pt"))) +
  panel_b + panel_c + panel_d + panel_e + panel_f +
  plot_layout(design = fig_layout) +
  plot_annotation(tag_levels = "a")

# Export figures
ggsave("SVGs/Figure_S3.svg", width=20, height=20, units="cm")
ggsave("Figure_S3.png", width=20, height=23, units="cm", dpi=300)

