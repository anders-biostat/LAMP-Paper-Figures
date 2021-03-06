library(tidyverse)
library(patchwork)
library(glue)
library(ggsci)
library(binom)
library(ggbeeswarm)

source("misc.R")

read_tsv( "data/tecan_values.tsv", col_types = "ccldccccdd" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00003", "CP00005", "CP00006", "CP00008", 
   "CP00009", "CP00010", "CP00011", "CP00012", "CP00013", "CP00016" )

lamp_thresh <- c(.3)

## panel a
tecan %>% 
  filter( plate %in% plates_to_use ) %>%
  filter( minutes==30, gene=="N" ) %>% 
  filter( ! ( plate == "CP00005" & well >= "G" ) ) %>%
  left_join( tblCT ) %>% 
  filter( !is.na(CT) ) %>%
  filter( is.na(wellRemark) ) -> tbl

set.seed(2020)
panel_a <- tbl %>%
  mutate(CT = ifelse(CT > 40, runif(n(), 43, 47), CT)) %>%
  mutate(plate = str_extract(plate, "\\d{2}$")) %>%
  mutate(extraction = recode(extraction, "Qiasymphonie" = "QS", "Qiacube" = "QC")) %>%
  sample_frac() %>%  # shuffle points to randomize overplotting at least
  ggplot() +
  geom_hline(yintercept = lamp_thresh, color = "lightgray" ) +
  #geom_hline(yintercept = lamp_thresh[[1]], color = "lightgray", linetype = "dashed" ) +
  geom_vline(xintercept = c(30, 41.5), color = "lightgray" ) +
  geom_point( aes( x = CT, y = absBlue - absYellow, fill = plate, shape = extraction), colour = "black", alpha = .6, size = 1.2 ) + 
  scale_x_reverse(breaks = c(20, 30, 40, 45), labels = c(20, 30, 40, "negative")) +
  scale_shape_manual(values = c("QS" = 21, "QC" = 24)) +
  scale_fill_d3(palette="category20") +
  guides(fill = guide_legend(override.aes = list(shape = 21), label = FALSE, ncol = 2, title.hjust = 0.5)) +
  #annotate("text", color = "gray50", x = 49, y = -.26, label = glue("negative"), angle = 90) +
  annotate("text", color = "gray50", x = 49, y = 0, label = glue("negative"), angle = 90) +
  #annotate("text", color = "gray50", x = 49, y = .125, label = glue("inconclusive"), angle = 90) +
  annotate("text", color = "gray50", x = 49, y = .425, label = glue("positive"), angle = 90) +
  labs(subtitle = glue("30 min at 65°C\n{nrow(tbl)} samples on {tbl%>%distinct(plate)%>%nrow} plates"),
       x = glue("RT-qPCR (CT value)"),
       y = "RT-LAMP (ΔOD)")
panel_a

tbl %>%
  ggplot() +
  geom_quasirandom(aes(x = plate, y = absBlue - absYellow, fill = plate, shape = extraction), varwidth = TRUE,
                   colour = "black", alpha = .6, size = 1.2) +
  scale_fill_d3(palette="category20", labels=rep("", tbl%>%distinct(plate)%>%nrow), guide = FALSE) +
  scale_shape_manual(values = c("Qiasymphonie" = 21, "Qiacube" = 24)) +
  theme(axis.text.x = element_blank())

## panel b
ct_breaks <- c(0, 30, 35, 40, Inf)

# classify samples into false positives, TP, FN, TN...
lamp_cls <- tbl %>%
  mutate(delta_abs = absBlue - absYellow,
         lamp_result = case_when(delta_abs > lamp_thresh ~ "positive", 
                                 delta_abs <= lamp_thresh ~ "negative",
                                 TRUE                       ~ "inconclusive" ),
         qpcr_result = if_else(CT == Inf, "negative", "positive")) %>%
  mutate(result = case_when(lamp_result == "negative" & qpcr_result == "negative" ~ "TN",
                            lamp_result == "positive" & qpcr_result == "positive" ~ "TP",
                            lamp_result == "positive" & qpcr_result == "negative" ~ "FP",
                            lamp_result == "negative" & qpcr_result == "positive" ~ "FN"))

ss_binned <- lamp_cls %>%
  mutate(ct_bin = cut(CT, ct_breaks)) %>%
  group_by(result, ct_bin) %>%
  tally() %>% ungroup() %>%
  pivot_wider(names_from = result, values_from = n, values_fill = list(n=0)) %>%
  mutate(n = FN + FP + TN + TP,
         sensitivity = TP / (TP + FN),
         sensitivity_ci_upper = binom.confint(TP, (TP + FN), method="wilson")$upper,
         sensitivity_ci_lower = binom.confint(TP, (TP + FN), method="wilson")$lower,
         specificity = TN / (TN + FP),
         specificity_ci_upper = binom.confint(TN, (TN + FP), method="wilson")$upper,
         specificity_ci_lower = binom.confint(TN, (TN + FP), method="wilson")$lower) %>%
  mutate(ct_bin = if_else(str_detect(ct_bin, "Inf"), "negative",
                          paste(str_extract(ct_bin, "(?<=,)\\d+"), str_extract(ct_bin, "\\d+"), sep="-")))

## number of LAMP `incl` samples, which will not be accounted for in sensitivity/specificity calculations
sum(ss_binned$`NA`)

panel_b_pos <- ss_binned %>%
  filter(!is.na(sensitivity)) %>%
  ggplot(aes(x = fct_rev(ct_bin), y = sensitivity, ymin = sensitivity_ci_lower, ymax = sensitivity_ci_upper, group = 1)) +
  geom_crossbar(fill="white", width=.7) +
  #geom_text(aes(y = -0.075, label = n), position = position_dodge(width=.7), size = 3) +
  labs(x="RT-qPCR (CT value)") +
  plot_layout(tag_level = "new")  # prevent making this panel c

panel_b_neg <- ss_binned %>%
  filter(is.na(sensitivity)) %>%
  ggplot(aes(x = ct_bin, y = specificity, ymin = specificity_ci_lower, ymax = specificity_ci_upper, group = 1)) +
  geom_crossbar(fill="white", width=.7) +
  #geom_text(aes(y = -0.075, label = n), position = position_dodge(width=.7), size = 3) +
  scale_x_discrete(labels = c("negative")) +
  labs(x="")

panel_b <- (panel_b_neg + panel_b_pos + plot_layout(widths = c(2, 9))) &
  coord_cartesian(ylim = c(-.08, 1)) &
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-.08,1), breaks=0:5/5) &
  theme(panel.grid.major.y = element_line(colour = "lightgrey"), panel.grid.minor.y = element_line(colour = "lightgrey"))

panel_a + panel_b +
  plot_layout(widths = c(10, 9)) +
  plot_annotation(tag_levels = c("A", NULL))

# Export figures
ggsave("SVGs/Figure_3.svg", width=20, height=10, units="cm")
ggsave("Figure_3.png", width=20, height=10, units="cm", dpi=300)

## Confusion matrix
## Alernative code, for manual comparison

tbl %>%
mutate( CTbin = cut( CT, ct_breaks ) ) %>% 
mutate_at( "CTbin", recode, 	"neg" = "(40,Inf]" ) %>%
mutate( LAMPres = 
   cut( absBlue-absYellow, c( -Inf, lamp_thresh, Inf ) ) %>%
   as.integer %>%
   { c( "neg", "pos" )[.] } ) %>%
group_by( CTbin, LAMPres ) %>%
count() %>%
pivot_wider( names_from = LAMPres, values_from = n, values_fill = c(n=0) )  ->
  confusion_matrix

bind_cols( confusion_matrix, ss_binned ) %>% write_tsv( "LAMP_confMatrix.tsv" )

