library( tidyverse )
source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "ccliccccdd" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT
read_tsv( "data/qPCR_CP00001_E_and_N.tsv" ) -> tblCT_extra

tblCT %>% 
filter( plate == "CP00001" ) %>%
left_join( tblCT_extra ) %>%
mutate_at( c( "vNCOVE", "vNCOVN" ), 
   ~ ifelse( .=="negativ", Inf, suppressWarnings(as.numeric(.)) ) ) %>%  
assertr::verify( ( !is.finite(CT) & !is.finite(vNCOVE) ) | abs( CT - vNCOVE ) < .2 ) %>%
select(-CT) -> qPCRpl1  

tecan %>%
  filter( plate == "CP00001", gene == "N", minutes==30 ) %>%
left_join( qPCRpl1 ) -> tbl

panel_a <- tbl %>%
  mutate_at( c( "vNCOVE", "vNCOVN" ), 
             ~ ifelse( . > 40 | is.na(.), runif( n(), 43, 47 ), suppressWarnings(as.numeric(.)) ) ) %>%
ggplot() +
  geom_vline( xintercept = 42, col="lightgray" ) +
  geom_point( aes( x = vNCOVN, y = absBlue - absYellow ), 
    fill = "black", colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
  scale_x_reverse( breaks = c( 10, 20, 30, 40, 45 ), labels = c( 10, 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value for gene 'N')" ) + 
  ylab( "RT-LAMP (ΔOD)" ) + theme_bw() + theme(panel.grid = element_blank())
panel_a

## calculation of mean difference (assumption of a slope of 1)
tbl_panel_b_diff <- tbl_panel_b %>% mutate(diff = vNCOVN - vNCOVE) %>% summarise(n = n(), mean = mean(diff), se = sd(diff) / sqrt(n))
tbl_panel_b_diff
# tbl_panel_b_lm <- summary(lm(vNCOVN ~ vNCOVE, tbl_panel_b))
# tbl_panel_b_lm$coefficients[1, 1:2]
# tbl_panel_b_lm$coefficients[1, 2] / tbl_panel_b_diff$n

tbl_panel_b <- tbl %>%
  filter_at( c( "vNCOVE", "vNCOVN" ), all_vars(!(. > 40 | is.na(.))) )
panel_b <- tbl_panel_b %>%
  ggplot() + 
  geom_point( aes(x=vNCOVE, y=vNCOVN, fill=vNCOVE), 
              colour = "black", alpha = .8, shape = 21, size = 1.6 ) +
  coord_fixed() +
  geom_abline() +
  scale_fill_ct(  name="CT_E" ) +
  xlab("RT-qPCR (CT value for E gene)" ) +
  ylab("RT-qPCR (CT value for N gene)" ) 
panel_b

panel_a + panel_b +
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(face = "bold"))

# Export figures
ggsave("SVGs/Figure_S1.svg", width=20, height=10, units="cm")
ggsave("Figure_S1.png", width=20, height=10, units="cm", dpi=300)
