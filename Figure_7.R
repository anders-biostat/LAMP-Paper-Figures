library( tidyverse )
library( patchwork )
library( glue )
source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "ccldccccdd" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00035", "CP00036", "CP00037" ) 

tecan %>% 
filter( plate %in% plates_to_use ) %>%  
filter( exclude != "yes" ) %>%
select( -plateRemark ) %>%
left_join( tblCT ) %>%
filter( !is.na(CT) ) %>%
mutate( CT = ifelse( CT>40, 46.5, CT ) ) %>%
filter( !is.na(CT) ) %>%
mutate( facet = ifelse( heat95, "5 min 95 °C prior to testing at 65 °C", "direct testing at 65 °C" ) ) -> tbl

n_65 <- tbl %>% filter( ! heat95) %>% distinct(well, plate) %>% nrow()
n_65_dstnct <- tbl %>% filter( ! heat95) %>% distinct(barcode) %>% nrow()
n_plates_65 <- tbl %>% filter( ! heat95) %>% distinct(plate) %>% nrow()
n_95 <- tbl %>% filter(heat95) %>% distinct(well, plate) %>% nrow()
n_95_dstnct <- tbl %>% filter(heat95) %>% distinct(barcode) %>% nrow()
n_plates_95 <- tbl %>% filter(heat95) %>% distinct(plate) %>% nrow()

ylimits <- range( tbl$absBlue - tbl$absYellow ) 

tbl %>%
  mutate( group = str_c( plate, well, heat95 ) ) %>%
  mutate( celsius = if_else(heat95,
                            glue("5 min 95°C prior to testing, 30 min at 65°C\n",
                                 "{n_95} aliquots from {n_95_dstnct} samples on {n_plates_95} plates"),
                            glue("direct testing, 30 min at 65°C\n",
                                 "{n_65} aliquots from {n_65_dstnct} samples on {n_plates_65} plates"))) %>%
  mutate(group = fct_reorder(group, desc(CT))) %>%
  ggplot() + 
    geom_line( aes( x=minutes, y=absBlue-absYellow, group=group, col=CT ), alpha=.4 ) +
    facet_grid( . ~ fct_rev(celsius), scales = "free_x", space = "free_x" ) +
    scale_color_ct( name="RT-qPCR\nCT value") +
    theme( legend.position = "none" ) +
    ylim( ylimits ) +
    labs(x = "minutes at 65 °C",
         y = expression( "RT-LAMP (ΔOD"["30 min"]~")" ) ) -> plot9a

## color-scale comparisons
# (tbl %>%
#   mutate( group = str_c( plate, well, heat95 ) ) %>%
#   mutate( celsius = if_else(heat95,
#                             glue("5 min 95°C prior to testing, 30 min at 65°C\n",
#                                  "{n_95} aliquots from {n_95_dstnct} samples on {n_plates_95} plates"),
#                             glue("direct testing, 30 min at 65°C\n",
#                                  "{n_65} aliquots from {n_65_dstnct} samples on {n_plates_65} plates"))) %>%
#   mutate(group = fct_reorder(group, desc(CT))) %>%
#   ggplot() + 
#   geom_line( aes( x=minutes, y=absBlue-absYellow, group=group, col=CT ), alpha=.4 ) +
#   facet_grid( . ~ fct_rev(celsius), scales = "free_x", space = "free_x" ) +
# scale_color_gradientn(
#   name="RT-qPCR\nCT value",
#         guide = guide_colourbar(reverse = TRUE),
#         breaks = rev(c(10, 15, 20, 25, 30, 35, 40, 45)),
#         labels = rev(c(10, 15, 20, 25, 30, 35, 40, "neg")),
#         colors = c(viridis::magma(100, direction=1L)[1:75],         # range of CP values
#                    rep("#FFFFFF", 5),                               # white spacer
#                    viridis::magma(100, direction=1L)[rep(90, 10)])  # negative
#       ) +
#   ylim( ylimits ) +
#   labs(x = "minutes at 65 °C",
#        y = expression( "RT-LAMP (ΔOD"["30 min"]~")" ) )) /
# tbl %>%
#   mutate( group = str_c( plate, well, heat95 ) ) %>%
#   mutate( celsius = if_else(heat95,
#                             glue("5 min 95°C prior to testing, 30 min at 65°C\n",
#                                  "{n_95} aliquots from {n_95_dstnct} samples on {n_plates_95} plates"),
#                             glue("direct testing, 30 min at 65°C\n",
#                                  "{n_65} aliquots from {n_65_dstnct} samples on {n_plates_65} plates"))) %>%
#   mutate(group = fct_reorder(group, desc(CT))) %>%
#   ggplot() + 
#   geom_line( aes( x=minutes, y=absBlue-absYellow, group=group, col=CT ), alpha=.4 ) +
#   facet_grid( . ~ fct_rev(celsius), scales = "free_x", space = "free_x" ) +
#   scale_color_ct( name="RT-qPCR\nCT value") +
#   ylim( ylimits ) +
#   labs(x = "minutes at 65 °C",
#        y = expression( "RT-LAMP (ΔOD"["30 min"]~")" ) )


tbl %>%
  mutate( diff = absBlue - absYellow ) %>%
  select( -absBlue, -absYellow ) %>%
  pivot_wider( names_from = minutes, values_from = diff ) %>%
  filter(!is.na(`10`) & !is.na(`30`)) -> tbl_wide

n_65 <- tbl_wide %>% filter( ! heat95) %>% distinct(well, plate) %>% nrow()
n_65_dstnct <- tbl_wide %>% filter( ! heat95) %>% distinct(barcode) %>% nrow()
n_plates_65 <- tbl_wide %>% filter( ! heat95) %>% distinct(plate) %>% nrow()
n_95 <- tbl_wide %>% filter(heat95) %>% distinct(well, plate) %>% nrow()
n_95_dstnct <- tbl_wide %>% filter(heat95) %>% distinct(barcode) %>% nrow()
n_plates_95 <- tbl_wide %>% filter(heat95) %>% distinct(plate) %>% nrow()

tbl_wide %>%
  mutate( celsius = if_else(heat95,
                            glue("5 min 95°C prior to testing, 30 min at 65°C\n",
                                 "{n_95} aliquots from {n_95_dstnct} samples on {n_plates_95} plates"),
                            glue("direct testing, 30 min at 65°C\n",
                                 "{n_65} aliquots from {n_65_dstnct} samples on {n_plates_65} plates"))) %>%
  ggplot() +
     geom_point( aes( `30`-`10`, y=`30`, fill = CT), colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
    facet_grid( . ~ fct_rev(celsius) ) +
    scale_fill_ct(name="RT-qPCR\nCT value") +
    #theme( legend.position = "none" ) +
    ylim( ylimits ) +
    labs(x = expression( "RT-LAMP ( ΔOD"["30 min"] - "ΔOD"["10 min"]~")" ),
         y = expression( "RT-LAMP ( ΔOD"["30 min"]~")" )) -> plot9b

fig9layout <- '
AAAAA
BBBBC
'
plot9a + plot9b + guide_area() +
  plot_layout(design = fig9layout, guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave("Figure_7.png", width=20, height=20, units="cm", dpi=300)
ggsave("SVGs/Figure_7.svg", width=20, height=20, units="cm")
