library( tidyverse )

trans_logap <- function( o=0 ) 
  scales::trans_new( "logap", 
    function(x) log10( x + 10^o), 
    function(x) 10^x - 10^o )

mybreaks <- c( 0, 10^(0:5) )
mylabels <- c( 0, 1, 10, expression(10^2), expression(10^4), expression(10^5), expression(10^5) )
  
load( "data/counts.Rda" )

counts %>% filter( matched, cpm>3000 ) %>% select( Plate, Well )-> pos

plates_to_use <- c( "CP00003", "CP00005", "CP00006", "CP00008", 
                    "CP00009", "CP00010", "CP00011", "CP00012", "CP00013", "CP00016" )

counts %>% 
select( -cpm ) %>% 
pivot_wider( names_from = matched, values_from = count, values_fill = c(count=0) ) %>%
left_join( pos %>% add_column( call="pos" ) ) %>%
mutate_at( "call", replace_na, "neg" ) %>%
rename( plate=Plate, well=Well ) %>%
mutate( plate = str_replace( plate, "Plt", "CP000" ) ) %>%
inner_join( tecan %>% filter( minutes==30 ) , by=c("plate","well") ) %>%
filter( plate %in% plates_to_use ) %>%
filter( !(plate == "CP00003" & plateRemark != "2")) %>%
filter( gene=="N" ) %>% 
filter( ! ( plate == "CP00005" & well >= "G" ) ) %>%
left_join( tblCT ) %>%
filter( !is.na(CT) ) -> tbl1



tbl1 %>%
mutate( facet = ifelse( absBlue-absYellow > .3, "ΔOD [30min] > .3", "ΔOD [30min] < .3" ) ) %>%
ggplot( aes( x = `TRUE`+`FALSE`,  y=`TRUE`/(`TRUE`+`FALSE`) ) ) + 
  geom_point( aes( col=absBlue-absYellow ), size=.5, alpha=.8 ) + 
  scale_x_continuous(trans=trans_logap(1), breaks=mybreaks) + 
  theme_bw() + theme( panel.grid = element_blank(), axis.text.x = element_text(angle = 90) ) +
  xlab( "UMI count for well" ) +
  ylab( "fraction of matching UMIs" ) +
  facet_wrap( ~ facet ) +
  scale_color_gradient2( midpoint=.2, low="purple", mid="darkgray", high="brown", 
    limits = c( -.5, .5), oob=scales::squish, name="ΔOD [30min]" )

tbl1 %>%
mutate( `RT-LAMP` = ifelse( absBlue-absYellow > .3, "ΔOD [30min] > .3", "ΔOD [30min] < .3" ) ) %>%
ggplot( aes( x = `TRUE`+`FALSE`,  y=`TRUE`/(`TRUE`+`FALSE`) ) ) + 
  geom_point( aes( col=`RT-LAMP` ), size=.5, alpha=.8 ) + 
  scale_x_continuous(trans=mytrans, breaks=mybreaks) + 
  theme_bw() + theme( panel.grid = element_blank(), axis.text.x = element_text(angle = 90) ) +
  xlab( "UMI count for well" ) +
  ylab( "fraction of matching UMIs" ) 

## This one!
tbl1 %>%
ggplot( aes( x = `TRUE`+`FALSE`,  y=`TRUE` ) ) + 
  geom_point( aes( col=CT ), size=.5, alpha=.8 ) + 
  scale_x_continuous(trans=mytrans, breaks=mybreaks, labels=mylabels ) + 
  scale_y_continuous(trans=mytrans, breaks=mybreaks, labels=mylabels) + 
  theme_bw() + theme( panel.grid = element_blank(), axis.text.x = element_text(angle = 90) ) +
  xlab( "reads total" ) +
  ylab( "reads matching viral amplicon" ) +
  scale_color_ct() +
  coord_fixed()

tbl1 %>%
ggplot( aes( x = `TRUE`+`FALSE`,  y=`TRUE` ) ) + 
  geom_point( aes( col=absBlue-absYellow ), size=.5, alpha=.8 ) + 
  scale_x_continuous(trans=mytrans, breaks=mybreaks, labels=mylabels ) + 
  scale_y_continuous(trans=mytrans, breaks=mybreaks, labels=mylabels) + 
  theme_bw() + theme( panel.grid = element_blank() ) +
  xlab( "reads total" ) +
  ylab( "reads matching viral amplicon" ) +
  scale_color_gradient2( midpoint=.2, low="purple", mid="darkgray", high="brown", 
    limits = c( -.5, .5), oob=scales::squish, name="ΔOD [30min]" ) +
  coord_fixed()

tbl1 %>%
mutate( facet = ifelse( absBlue-absYellow > .3, "ΔOD [30min] > .3", "ΔOD [30min] < .3" ) ) %>%
ggplot( aes( x = `TRUE`+`FALSE`,  y=`TRUE`/(`TRUE`+`FALSE`) ) ) + 
  geom_point( aes( col=call ), size=.5, alpha=.8 ) + 
  scale_x_continuous(trans=mytrans, breaks=mybreaks) + 
  theme_bw() + theme( panel.grid = element_blank(), axis.text.x = element_text(angle = 90) ) +
  xlab( "UMI count for well" ) +
  ylab( "fraction of matching UMIs" ) +
  facet_wrap( ~ facet ) 

tbl1 %>%
ggplot( aes( x = `TRUE`+`FALSE`,  y=`TRUE`/(`TRUE`+`FALSE`) ) ) + 
  geom_point( aes( col=absBlue-absYellow ), size=.3 ) + 
  geom_point( size=.7, col="orange",
    data = ( tbl1 %>% filter( plate=="Plt12", well %in% c( "A5", "A6" ) ) ) ) + 
  scale_x_continuous( trans=mytrans, breaks=c(0,1,10,1e2,2e3,1e4,1e5) ) + 
  theme_bw() + theme( panel.grid = element_blank() ) +
  xlab( "UMIs not mapping to virus" ) +
  ylab( "UMIs mapping to virus" ) +
  facet_wrap( ~ absBlue-absYellow > .3 ) 

tbl1 %>% 
mutate( total = `TRUE`+`FALSE` ) %>% 
ggplot( aes( x=plate, y=total ) ) + 
  ggbeeswarm::geom_beeswarm( aes( col=absBlue-absYellow, label=well ), size=.1 ) + 
  scale_y_continuous(trans=mytrans, breaks=mybreaks) + 
  theme( panel.grid = element_blank(), axis.text.x = element_text(angle = 90) ) + 
  ylab("total UMI count") +
  scale_color_gradient2( midpoint=.2, low="purple", mid="darkgray", high="brown", 
    limits = c( -.5, .5), oob=scales::squish, name="ΔOD [30min]" )
plotly::ggplotly()

tbl1 %>% 
mutate( total = `TRUE`+`FALSE` ) %>% 
ggplot( aes( x=plate, y=total ) ) + 
  ggbeeswarm::geom_beeswarm( aes( col=(WellRow=="A" & WellCol<=6 & WellCol>1) ), size=.1 ) + 
  scale_y_continuous(trans=mytrans, breaks=mybreaks) + 
  theme( panel.grid = element_blank(), axis.text.x = element_text(angle = 90) ) + 
  ylab("total UMI count")
