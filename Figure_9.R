library( tidyverse )
library( patchwork )
library( ggtext )
source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "cclicccdd" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT

plates_to_use <- c( "CP00035", "CP00036", "CP00037" ) 

tecan %>% 
filter( plate %in% plates_to_use ) %>% 
filter( !( plate == "CP00036" & minutes==30 & plateRemark != "30 mins 2nd scan" ) ) %>%
left_join( tblCT ) %>%
filter( !is.na(CT) ) %>%
mutate( CT = ifelse( CT>40, 45, CT ) ) %>%
filter( !is.na(CT) ) %>%
mutate( facet = ifelse( heat95, "5 min 95 °C prior to testing at 65 °C", "direct testing at 65 °C" ) ) -> tbl

tbl %>%
mutate( group = str_c( plate, well, heat95 ) ) %>%
arrange( -CT ) %>%
mutate_at( "group", fct_inorder ) %>% 
ggplot + 
  geom_line( aes( x=minutes, y=absBlue-absYellow, group=group, col=CT ), alpha=.4 ) +
  facet_grid( . ~ fct_rev(facet), scales = "free_x" ) +
  scale_color_ct( name="RT-qPCR\nCT value") +
  labs(x = "minutes at 65 °C",
       y = expression( "LAMP (ΔOD"["30 min"]~")" ) ) -> plot9a

tbl %>%
mutate( diff = absBlue - absYellow ) %>%
select( -absBlue, -absYellow ) %>%
pivot_wider( names_from = minutes, values_from = diff ) %>%
ggplot +
  geom_point( aes( x=`30`-`10`, y=`30`, col=CT ), size=.8 ) +
  facet_grid( . ~ fct_rev(facet) ) +
  scale_color_ct() +
  theme( legend.position = "none" ) +
  coord_fixed() +
  labs(x = expression( "LAMP ( ΔOD"["30 min"] - "ΔOD"["10 min"]~")" ),
       y = expression( "LAMP ( ΔOD"["30 min"]~")" )) -> plot9b


plot9a / plot9b +
  plot_annotation(tag_levels = "a")


ggsave("SVGs/Figure_9.svg", width=20, height=20, units="cm")
ggsave("Figure_9.png", width=15, height=15, units="cm", dpi=300)
