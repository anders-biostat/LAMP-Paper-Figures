library( tidyverse )
library( patchwork )
library( jpeg )
library( grid )
source( "misc.R" )

read_tsv( "data/tecan_values.tsv", col_types = "ccliccccdd" ) -> tecan
read_tsv( "data/plates_with_CTs.tsv" ) -> tblCT


## Figure 3b 
tecan %>% 
left_join( tblCT ) %>%
mutate( CT = ifelse( CT>40, 46.5, CT ) ) %>%
filter( !is.na(CT) ) %>%
filter( plate %in% "CP00001", gene=="N", plateRemark=="2" ) %>% 
mutate(well = fct_reorder(well, desc(CT))) %>%
ggplot +
  geom_line( aes( x=minutes, y=absBlue-absYellow, group=well, col=CT ), size=.5 ) +
  scale_color_ct() +
  xlab( "incubation time [minutes]" ) + 
  ylab( expression( "RT-LAMP (ΔOD)" ) ) -> plot3b
plot3b

## Figure 3c
tbl_3c <- tecan %>% 
  filter( plate == "CP00001" ) %>%
  filter( ( gene=="N" & minutes == 30 ) | ( gene=="1a" & minutes==40 ) ) %>%
  left_join( tblCT ) %>%
  filter( !is.na(CT) ) %>%
  mutate( CT = ifelse( CT>40, runif( n(), 43, 47 ), CT ) ) %>%
  mutate( facet = str_c( gene, " gene, ", minutes, " min" ) )

tbl_3c %>%
  filter(gene == "N") %>%
  ggplot() +
  geom_vline( xintercept = 42, col="darkgray" ) +
  geom_point( aes( x = CT, y = absBlue - absYellow ), fill = "black", colour = "black", alpha = .6, shape = 21, size = 1.2 ) + 
  scale_x_reverse( breaks = c( 10, 20, 30, 40, 45 ), labels = c( 10, 20, 30, 40, "neg" ) ) +
  xlab( "RT-qPCT (CT value)" ) + 
  ylab( "RT-LAMP (ΔOD)" ) -> plot3c

fig3a <- readJPEG("SVGs/Figure_3a.jpeg")
#svg("SVGs/Figure_3.svg", width=20/2.54, height=14/2.54) ## gives strange light borders arround CT colorscale gradient
wrap_elements(panel = rasterGrob(fig3a, x = 0, y = 0.5, just = "left", interpolate = TRUE)) / (plot3b | plot3c) +
  plot_annotation(tag_levels = "a")
#dev.off()

ggsave("SVGs/Figure_3.svg", width=20, height=14, units="cm", dpi = 300) # results in low resolution raster image (https://github.com/r-lib/svglite/issues/86)
ggsave("Figure_3.png", width=20, height=14, units="cm", dpi=300)
