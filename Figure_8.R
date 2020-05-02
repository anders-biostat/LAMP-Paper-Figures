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

plates_to_use <- c( "CP10020", "CP10021", "CP00024", "CP00025", "CP00030", 
  "CP00035", "CP00036", "CP00037" ) 

deltaOD_cutoff_zeroL <- 0.25

tecan %>%
  filter( plate %in% plates_to_use ) %>%
  filter( !( plate == "CT00036" & minutes==30 & plateRemark != "30 mins 2nd scan" ) ) %>%
  filter( !( plate == "CT10020" & as.integer(str_sub(well,1,-1)) > 4 ) ) %>%
  left_join( tblCT, by = c("plate", "well") ) %>%
  filter( !is.na(CT) & minutes == 30 ) -> tbl


# Figure 8a
n_65 <- tbl %>% filter( ! heat95) %>% nrow()
n_plates_65 <- tbl %>% filter( ! heat95) %>% distinct(plate) %>% nrow()
n_95 <- tbl %>% filter(heat95) %>% nrow()
n_plates_95 <- tbl %>% filter(heat95) %>% distinct(plate) %>% nrow()

fig8a <- tbl %>%
  mutate( celsius = if_else(heat95,
                            glue("5 min 95°C prior to testing, 30 min at 65°C\n{n_95} samples on {n_plates_95} plates"),
                            glue("direct testing, 30 min at 65°C\n{n_65} samples on {n_plates_65} plates"))) %>%
  mutate( CT = ifelse( CT > 40, runif( n(), 43, 47 ), CT ) ) %>% 
  mutate_at( "plate", fct_relevel,  "CP10020", "CP10021" ) %>% # to get plates in chronological order
  sample_frac() %>%  # randomize order of points so that it's not one specific plate plotted on top of all others
  ggplot() +
  geom_hline( yintercept = deltaOD_cutoff_zeroL, col="lightgray" ) +
  geom_vline( xintercept = c( 30, 42 ), col="lightgray" ) +
  geom_point( aes( x=CT, y=absBlue - absYellow, col=plate), size = .8 ) +
  facet_wrap( ~ fct_rev(celsius)) +
  scale_x_reverse(breaks = c(20, 30, 40, 45), labels = c(20, 30, 40, "negative")) +
  scale_color_d3(palette="category10") +
  labs(y = expression("ΔOD"["30 min"]),
       x = "RT-qPCR (CT value)") +
  annotate("text", color = "gray50", x = 50, y = -.06, label = "negative", angle = 90) +
  annotate("text", color = "gray50", x = 50, y = .44, label = "positive", angle = 90)


# Fig 8b
# Sensitivity + Specificity of Zero Lysis plates
# Treating replicates as independent samples

# classify samples into false positives, TP, FN, TN...
lamp_cls <- tbl %>%
  mutate(delta_abs = absBlue - absYellow,
         lamp_result = if_else(delta_abs > deltaOD_cutoff_zeroL, "positive", "negative"),
         qpcr_result = if_else(CT == Inf, "negative", "positive")) %>%
  mutate(result = case_when(lamp_result == "negative" & qpcr_result == "negative" ~ "TN",
                            lamp_result == "positive" & qpcr_result == "positive" ~ "TP",
                            lamp_result == "positive" & qpcr_result == "negative" ~ "FP",
                            lamp_result == "negative" & qpcr_result == "positive" ~ "FN"))

lamp_cls %>%
  group_by(result, heat95) %>%
  tally() %>% ungroup() %>%
  pivot_wider(names_from = result, values_from = n) %>%
  mutate(sensitivity = TP / (TP + FN),
         sensitivity_ci_upper = binom.confint(TP, (TP + FN), method="wilson")$upper,
         sensitivity_ci_lower = binom.confint(TP, (TP + FN), method="wilson")$lower,
         specificity = TN / (TN + FP),
         specificity_ci_upper = binom.confint(TN, (TN + FP), method="wilson")$upper,
         specificity_ci_lower = binom.confint(TN, (TN + FP), method="wilson")$lower)

ct_breaks <- c(0, 25, 30, 35, 40, Inf)

ss_binned <- lamp_cls %>%
  mutate(ct_bin = cut(CT, ct_breaks)) %>%
  group_by(result, heat95, ct_bin) %>%
  tally() %>% ungroup() %>%
  pivot_wider(names_from = result, values_from = n, values_fill = list(n=0)) %>%
  mutate(sensitivity = TP / (TP + FN),
         sensitivity_ci_upper = binom.confint(TP, (TP + FN), method="wilson")$upper,
         sensitivity_ci_lower = binom.confint(TP, (TP + FN), method="wilson")$lower,
         specificity = TN / (TN + FP),
         specificity_ci_upper = binom.confint(TN, (TN + FP), method="wilson")$upper,
         specificity_ci_lower = binom.confint(TN, (TN + FP), method="wilson")$lower)

p_pos_zl <- ss_binned %>%
  filter(!is.na(sensitivity)) %>%
  mutate(celsius = fct_rev(if_else(!heat95, "direct testing", "95°C treatment"))) %>%
  ggplot(aes(x = ct_bin, y = sensitivity, ymin = sensitivity_ci_lower, ymax = sensitivity_ci_upper, color = celsius, group = celsius)) +
  geom_crossbar(fill="white", position = position_dodge(width=.6), width=.5) +
  scale_x_discrete(labels = c("0-25", "26-30", "31-35", "36-40")) +
  labs(x = "RT-qPCR CT value",
       color = "") +
  theme(legend.position = c(0.8, 0.9), legend.background = element_blank()) +
  guides(color = guide_legend(label.position = "left", label.hjust = 1))

p_neg_zl <- ss_binned %>%
  filter(is.na(sensitivity)) %>%
  mutate(celsius = fct_rev(if_else(!heat95, "direct testing", "95°C treatment"))) %>%
  ggplot(aes(x = ct_bin, y = specificity, ymin = specificity_ci_lower, ymax = specificity_ci_upper, color = celsius, group = celsius)) +
  geom_crossbar(fill="white", position = position_dodge(width=.6), width=.5) +
  scale_x_discrete(labels = c("negative")) +
  labs(x="") +
  theme(legend.position = "none")

fig8b <- (p_pos_zl + p_neg_zl + plot_layout(widths = c(4, 1))) &
  coord_cartesian(ylim = c(0, 1)) &
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1), breaks=0:5/5) &
  scale_color_brewer(palette = "Set1", direction = -1L)


fig8a / fig8b +
  plot_layout(heights = c(3, 2)) +
  plot_annotation(tag_levels = "a")

ggsave("Figure_8.png", width=20, height=14.5, units="cm", dpi=400)
ggsave("SVGs/Figure_8.svg", width=20, height=14.5, units="cm")
