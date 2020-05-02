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

plates_to_use <- c( "CP10020", "CP10021", "CP00024", "CP00025",# "CP00030", 
  "CP00035", "CP00036", "CP00037" ) 

deltaOD_cutoff_zeroL <- 0.3

tecan %>%
  filter( plate %in% plates_to_use ) %>%
  filter( exclude != "yes" ) %>%
  filter( !( plate == "CP10020" & as.integer(str_sub(well,2,-1)) > 4 ) ) %>% 
  left_join( tblCT, by = c("plate", "well") ) %>% 
  filter( !is.na(CT) & minutes == 30 ) -> tbl

# Check that we have each plate only once
tbl %>%
  group_by( plate, heat95, well ) %>%
  count() %>% group_by(n) %>% filter(n>1) 

# Figure 8a
n_65 <- tbl %>% filter( ! heat95) %>% nrow()
n_65_dstnct <- tbl %>% filter( ! heat95) %>% distinct(barcode) %>% nrow()
n_plates_65 <- tbl %>% filter( ! heat95) %>% distinct(plate) %>% nrow()
n_95 <- tbl %>% filter(heat95) %>% nrow()
n_95_dstnct <- tbl %>% filter(heat95) %>% distinct(barcode) %>% nrow()
n_plates_95 <- tbl %>% filter(heat95) %>% distinct(plate) %>% nrow()

fig8a <- tbl %>%
  mutate( celsius = if_else(heat95,
                            glue("5 min 95°C prior to testing, 30 min at 65°C\n",
                              "{n_95} aliquots from {n_95_dstnct} samples on {n_plates_95} plates"),
                            glue("direct testing, 30 min at 65°C\n",
                              "{n_65} aliquots from {n_65_dstnct} samples on {n_plates_65} plates"))) %>%
  mutate( CT = ifelse( CT > 40, runif( n(), 43, 47 ), CT ) ) %>% 
  mutate_at( "plate", fct_relevel,  "CP10020", "CP10021" ) %>% # to get plates in chronological order
  sample_frac() %>%  # randomize order of points so that it's not one specific plate plotted on top of all others
  ggplot() +
  geom_hline( yintercept = deltaOD_cutoff_zeroL, col="lightgray" ) +
  geom_vline( xintercept = c( 30, 42 ), col="lightgray" ) +
  geom_point( aes( x = CT, y = absBlue - absYellow, fill = plate), colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
  facet_wrap( ~ fct_rev(celsius)) +
  scale_x_reverse(breaks = c(20, 30, 40, 45), labels = c(20, 30, 40, "negative")) +
  scale_fill_d3(palette="category10", labels=rep("", 99)) +
  labs(y = expression( "RT-LAMP ( ΔOD"["30 min"]~")" ),
       x = "RT-qPCR (CT value)") +
  annotate("text", color = "gray50", x = 50, y = 0, label = "negative", angle = 90) +
  annotate("text", color = "gray50", x = 50, y = .47, label = "positive", angle = 90)


# Fig 8b
# Sensitivity + Specificity of Zero Lysis plates
# Treating replicates as independent samples

# classify samples into false positives, TP, FN, TN...
ct_breaks <- c(0, 25, 30, 35, 40, Inf)
lamp_cls <- tbl %>%
  mutate(delta_abs = absBlue - absYellow,
         lamp_result = if_else(delta_abs > deltaOD_cutoff_zeroL, "positive", "negative"),
         qpcr_result = if_else(CT == Inf, "negative", "positive")) %>%
  mutate(result = case_when(lamp_result == "negative" & qpcr_result == "negative" ~ "TN",
                            lamp_result == "positive" & qpcr_result == "positive" ~ "TP",
                            lamp_result == "positive" & qpcr_result == "negative" ~ "FP",
                            lamp_result == "negative" & qpcr_result == "positive" ~ "FN")) %>%
  mutate(ct_bin = cut(CT, ct_breaks)) 


ss_binned <- lamp_cls %>%
  group_by(result, heat95, ct_bin) %>%
  tally() %>% ungroup() %>%
  pivot_wider(names_from = result, values_from = n, values_fill = list(n=0)) %>%
  mutate(sensitivity = TP / (TP + FN),
         sensitivity_ci_upper = binom.confint(TP, (TP + FN), method="wilson")$upper,
         sensitivity_ci_lower = binom.confint(TP, (TP + FN), method="wilson")$lower,
         specificity = TN / (TN + FP),
         specificity_ci_upper = binom.confint(TN, (TN + FP), method="wilson")$upper,
         specificity_ci_lower = binom.confint(TN, (TN + FP), method="wilson")$lower) %>%
  mutate(ct_bin = if_else(str_detect(ct_bin, "Inf"), "negative",
                          paste(str_extract(ct_bin, "(?<=,)\\d+"), str_extract(ct_bin, "\\d+"), sep="-")))

p_pos_zl <- ss_binned %>%
  filter(!is.na(sensitivity)) %>%
  mutate(celsius = fct_rev(if_else(!heat95, "direct testing", "95°C treatment"))) %>%
  ggplot(aes(x = fct_rev(ct_bin), y = sensitivity, ymin = sensitivity_ci_lower, ymax = sensitivity_ci_upper, color = celsius, group = celsius)) +
  geom_crossbar(fill="white", position = position_dodge(width=.6), width=.5) +
  labs(x = "RT-qPCR CT value", color = "") +
  theme(legend.position = c(0.17, 0.85), legend.background = element_blank(), legend.box.background = element_blank(), legend.key=element_blank()) +
  #guides(color = guide_legend(label.position = "left", label.hjust = 1)) +
  plot_layout(tag_level = "new")

p_neg_zl <- ss_binned %>%
  filter(is.na(sensitivity)) %>%
  mutate(celsius = fct_rev(if_else(!heat95, "direct testing", "95°C treatment"))) %>%
  ggplot(aes(x = ct_bin, y = specificity, ymin = specificity_ci_lower, ymax = specificity_ci_upper, color = celsius, group = celsius)) +
  geom_crossbar(fill="white", position = position_dodge(width=.6), width=.5) +
  labs(x="") +
  theme(legend.position = "none")

fig8b <- (p_neg_zl + p_pos_zl + plot_layout(widths = c(1, 4))) &
  theme(panel.grid.major.y = element_line(colour = "lightgrey")) &
  coord_cartesian(ylim = c(0, 1)) &
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1), breaks=0:5/5) &
  scale_color_brewer(palette = "Set1", direction = -1L)


fig8a / fig8b +
  plot_layout(heights = c(3, 2)) +
  plot_annotation(tag_levels = "a")

ggsave("Figure_8.png", width=20, height=14.5, units="cm", dpi=400)
ggsave("SVGs/Figure_8.svg", width=20, height=14.5, units="cm")


lamp_cls %>% 
  group_by( heat95, ct_bin, lamp_result ) %>%
  tally() %>%
  pivot_wider( names_from = lamp_result, values_from = n, values_fill = c(n=0) ) %>%
  left_join( ss_binned ) %>% write_csv( "LAMP_zeroLysis_confMatrix.tsv" )
