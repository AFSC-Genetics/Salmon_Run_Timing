library(ggmap)
library(tidyverse)
library(sf)
library(patchwork)

#Load Metadata
SalmonCollections<-read_csv("./Data/Raw/CollectionInfo.csv")

#Get some shapes of US and Canada from rNaturalEarth
BC_coast_sf <- rnaturalearth::ne_countries(country="Canada", scale=10,returnclass = "sf")

BC_coast_st <- st_transform(BC_coast_sf,
                            st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))%>%
  
  {st_crop(x=., y = c(xmin = -179, xmax= -110,ymin = 35, ymax= 70))}

BC_coast_sf <- st_as_sf(BC_coast_st)

US_coast_sf <- rnaturalearth::ne_countries(country="United States of America", scale=10,returnclass = "sf") 

US_coast_st <- st_transform(US_coast_sf,
                            st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")) %>%
  
  {st_crop(x=., y = c(xmin = -179, xmax= -110,ymin = 35, ymax= 70))}

US_coast_sf <- st_as_sf(US_coast_st)


#Get some rivers and creeks
#These are available from ADF&G Anadromous Waters Cat. Geodatabase
#https://www.adfg.alaska.gov/sf/SARR/AWC/index.cfm?ADFG=maps.dataFiles

st_layers(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip")

fc <- sf::st_read(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip", 
                  layer = "AWC_stream")

fc2 <- st_transform(fc,
                    st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

fc3<- fortify(fc2)


fc_AWC_pointCrp<-fc3%>%
  filter(NAME %in% c("Kenai River","Skilak River","Wood River","*Wood River distributary",
                     "Auke Nu Creek",
                     "Lake Creek",
                     "*Lake Two Creek",
                     "Waydelich Creek",
                     "Bay Creek", "Cali Creek",
                     "Auke Lake","Auke Creek"))%>%
  filter(AWC_CODE != "101-75-10300-2014")

rivers10 <- rnaturalearth::ne_download(scale = 10, type = 'rivers_lake_centerlines', 
                                       category = 'physical',
                                       returnclass = "sf")


#Yukon River
YukonKoyR <- rivers10 %>% 
  filter(grepl(pattern="Yukon|Koyukuk",name))

pChum <- ggplot() +
  geom_sf(data = BC_coast_sf, inherit.aes = FALSE,color="black",fill="#cccccc")+
  geom_sf(data = US_coast_sf, inherit.aes = FALSE,color="black",fill="#cccccc")+
  geom_sf(data = fc_AWC_pointCrp, inherit.aes = FALSE,color="blue",alpha=0.25)+
  geom_sf(data = YukonKoyR, inherit.aes = FALSE,color="blue",alpha=0.25)+
  coord_sf(xlim = c(-179, -125),  ylim = c(46.5, 73),expand = FALSE)+
  geom_point(data=SalmonCollections %>%
               mutate(Species = factor(Species,levels=c("Chum","Coho","Sockeye","Pink"))),aes(x=Long,y=Lat,fill=Species),pch=21,size=4)+
  xlab("")+
  ylab("")+
  theme(legend.position="right",
        text = element_text(size=30),
        axis.text=element_text(size=10),
        plot.margin=grid::unit(c(1,1,0.5,0.5), "mm"))+#trbl
  scale_y_continuous(breaks = seq(45, 70, by = 10)) + 
  scale_x_continuous(breaks = seq(-170, -125, by = 20))+
  ggspatial::annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(0.01, "in"), pad_y = unit(0.1, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20"
    )
  )+
  ggspatial::annotation_scale(
    location = "bl",
    bar_cols = c("grey60", "white")
  )+
  annotate("text",y=65,x= -150.0,label="Yukon R.",
           fontface = 'italic',angle = 33,size=3,color='grey50')+
  geom_text(data=SalmonCollections %>% 
              mutate(Species = factor(Species,levels=c("Chum","Coho","Sockeye","Pink")))%>%
              filter(Species == "Chum"),aes(x=Long+2,y=Lat+2,label=Label),size=3,fontface = "italic",angle=45)

#### Coho Skilak R ####
st_layers(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip")

fc <- sf::st_read(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip", 
                  layer = "AWC_stream")

fc2 <- st_transform(fc,
                    st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

fc3<- fortify(fc2)


Skilak_R<-fc3%>%
  filter(NAME %in% c("Kenai River","Skilak River","Soldotna River",
                     "Killey River", "Moose River","Slikok River",
                     "Beaver Creek","Funny River"))

st_layers(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip")

fc <- sf::st_read(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip", 
                  layer = "AWC_lake")

Skilak_L <- st_transform(fc,
                         st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")) 


pSkilak <- ggplot() +
  geom_sf(data = US_coast_sf, inherit.aes = FALSE,color="black",fill="#cccccc")+
  geom_sf(data = Skilak_R, inherit.aes = FALSE,color="blue",alpha=0.25)+
  geom_sf(data = Skilak_L, inherit.aes = FALSE,color="blue",fill="lightblue")+
  coord_sf(xlim = c(-151.45, -149.95),  ylim = c(60.0, 61.0),expand = FALSE)+
  geom_point(data=SalmonCollections %>%
               mutate(Species = factor(Species,levels=c("Chum","Coho","Sockeye","Pink"))),aes(x=Long,y=Lat,fill=Species),pch=21,size=4)+
  xlab("")+
  ylab("")+
  theme(legend.position="right",
        text = element_text(size=30),
        axis.text=element_text(size=10),
        plot.margin=grid::unit(c(1,1,0.5,0.5), "mm"))+#trbl
  scale_y_continuous(breaks = seq(60.3, 60.6, by = .3)) + 
  scale_x_continuous(breaks = seq(-152, -149, by = 1))+
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey60", "white")
  ) +
  annotate("text",y=60.422,x= -150.35,label="Skilak L.",
           fontface = 'italic',angle = -35,size=3,color='grey50')+
  annotate("text",y=60.5,x= -151.35,label="Cook Inlet",
           fontface = 'italic',angle = 90,size=3,color='grey50')+
  annotate("text",y=60.45,x= -151.05,label="Kenai R.",
           fontface = 'italic',angle = -12,size=3,color='grey50')+
  geom_text(data=SalmonCollections %>% 
              mutate(Species = factor(Species,levels=c("Chum","Coho","Sockeye","Pink")))%>%
              filter(Species == "Coho"),aes(x=Long+0.065,y=Lat+0.025,label=Label),size=3,fontface = "italic")

#Sockeye 
st_layers(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip")

fc <- sf::st_read(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip", 
                  layer = "AWC_stream")

fc2 <- st_transform(fc,
                    st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

fc3<- fortify(fc2)


Wood_R<-fc3%>%
  filter(grepl(pattern = "Wood|Aleknagik|Agulowak|Whitefish|Anvil|Teal|Pick|Eagel|Hansen|Happy|Ice|Bear|Yako|Fenno|Lynx|Elva|Hidden",NAME))

st_layers(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip")

fc <- sf::st_read(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip", 
                  layer = "AWC_lake")

Nerka_L <- st_transform(fc,
                        st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")) 


pNerka <- ggplot() +
  geom_sf(data = US_coast_sf, inherit.aes = FALSE,color="black",fill="#cccccc")+
  geom_sf(data = Wood_R, inherit.aes = FALSE,color="blue",alpha=0.25)+
  geom_sf(data = Nerka_L, inherit.aes = FALSE,color="blue",fill="lightblue")+
  coord_sf(xlim = c(-158.5, -159.25),  ylim = c(59.2, 59.625),expand = FALSE)+
  geom_point(data=SalmonCollections %>%
               mutate(Species = factor(Species,levels=c("Chum","Coho","Sockeye","Pink"))),aes(x=Long,y=Lat,fill=Species),pch=21,size=4)+
  xlab("")+
  ylab("")+
  theme(legend.position="right",
        text = element_text(size=30),
        axis.text=element_text(size=10),
        plot.margin=grid::unit(c(1,1,0.5,0.5), "mm"))+#trbl
  scale_y_continuous(breaks = seq(59.2, 59.6, by = .2)) + 
  scale_x_continuous(breaks = seq(-159.1,-158.6, by = 0.4))+
  ggspatial::annotation_scale(
    location = "bl",
    bar_cols = c("grey60", "white")
  )+
  annotate("text",y=59.4325,x= -158.65,label="L. Nerka",
           fontface = 'italic',angle = -35,size=3,color='grey50')+
  annotate("text",y=59.338,x= -158.8,label="L. Aleknagik",
           fontface = 'italic',angle = -35,size=3,color='grey50')+
  annotate("text",y=59.25,x= -158.585,label="Wood R.",
           fontface = 'italic',angle = -80,size=3,color='grey50')+
  geom_text(data=SalmonCollections %>% 
              mutate(Species = factor(Species,levels=c("Chum","Coho","Sockeye","Pink")))%>%
              filter(Species == "Sockeye"),aes(x=Long-0.025,y=Lat-0.012,label=Label),size=3,fontface = "italic")

# Pink Salmon
#load("./Data/Mapping/AWC_data.RData") #fc_AWC_point lives here

### New 
st_layers(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip")

fc <- sf::st_read(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip", 
                  layer = "AWC_Point")

fc2 <- st_transform(fc,
                    st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

fc3<- fortify(fc2)

Auke_Cr<-fc3%>%
  filter(NAME %in% c("Auke Nu Creek",
                     "Lake Creek",
                     "*Lake Two Creek",
                     "Waydelich Creek",
                     "Bay Creek", "Cali Creek",
                     "Auke Lake","Auke Creek"))%>%
  filter(NAME == "Auke Nu Creek" & TYPE== "UPPER" | NAME == "Auke Creek" & TYPE== "LOWER" |
           NAME == "Auke Lake" & TYPE== "LAKE" | NAME == "Bay Creek" & TYPE== "UPPER" |
           NAME == "Lake Creek" & TYPE== "MIDB" | NAME == "*Lake Two Creek" & TYPE== "UPPER" |
           NAME == "Waydelich Creek" & TYPE== "UPPER")%>%
  mutate(NewLab = gsub(pattern="Creek",replacement="Cr.",NAME))%>%
  filter(AWC_CODE != "101-75-10300-2014") %>%
  filter(!grepl(pattern="ANC",QUAD))

st_layers(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip")

fc <- sf::st_read(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip", 
                  layer = "AWC_lake")

Auke_L <- st_transform(fc,
                        st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")) 

fc_AWC_pointCrp2 <- Auke_Cr %>%
  mutate(my_nudge_x = c(0,0.002,0.0025,0,0.008,0.005,-0.0035),
         my_nudge_y = c(0,0.0005,0.0005,0.0005,0.002,0,0))%>%
  filter(NAME %in% c("Lake Creek","Auke Lake"))

fc <- sf::st_read(dsn = "./Data/Mapping/AWC_geodatabase/AWC2024.gdb.zip", 
                  layer = "AWC_stream")

fc2 <- st_transform(fc,
                    st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

fc3<- fortify(fc2)

Auke_Cr<-fc3%>%
  filter(NAME %in% c("Auke Nu Creek",
                     "Lake Creek",
                     "*Lake Two Creek",
                     "Waydelich Creek",
                     "Bay Creek", "Cali Creek",
                     "Auke Lake","Auke Creek"))

##Stat areas for white space
fc <- sf::st_read(dsn = "./Data/Mapping/SEAK_SalmonStatAreas.gdb.zip",
                  layer = "PVB_Southeast_2018_WGS1984")

fc2 <- st_transform(fc,
                    st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))


pAuke <- ggplot() +
  geom_sf(data = (Auke_Cr%>%
                    st_transform(.,
                                 st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))), inherit.aes = FALSE,color="blue",alpha=0.75)+
  geom_sf(data = (Auke_L%>%
                    st_transform(.,
                                 st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))), inherit.aes = FALSE,color="blue",fill="lightblue")+
  geom_sf(data = fc2,color='black', fill="white", inherit.aes = FALSE)+
  geom_sf_text(data = fc_AWC_pointCrp2 ,aes(label=NewLab), 
               inherit.aes = FALSE,color="black",alpha=0.5,size=3,
               nudge_x=fc_AWC_pointCrp2$my_nudge_x,nudge_y=fc_AWC_pointCrp2$my_nudge_y)+
  geom_point(data=SalmonCollections %>%
               mutate(Species = factor(Species,levels=c("Chum","Coho","Sockeye","Pink"))),aes(x=Long,y=Lat,fill=Species),pch=21,size=4)+
  coord_sf(xlim = c(-134.67, -134.62),  ylim = c(58.370, 58.4),expand = FALSE)+
  annotate("text",x=-134.66,y=58.3755,label="Auke\nBay",color="grey",fontface='italic')+
  xlab("")+
  ylab("")+
  theme(legend.position="none",
        text = element_text(size=30),
        axis.text=element_text(size=10))+
  scale_y_continuous(breaks = seq(58.37, 58.4, by = .02)) + 
  scale_x_continuous(breaks = seq(-134.7, -134.6, by = .02))+
  ggspatial::annotation_scale(
    location = "bl",
    bar_cols = c("grey60", "white")
  ) +
  geom_point(aes(y=58.380462,x= -134.642068), pch=25,fill="black",color="black")+#Add a point
  annotate("text",y=58.3798,x= -134.641, colour="black",label="Weir",size=3)+
  theme(panel.background = element_rect(fill = "#cccccc",
                                        colour = "#cccccc",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "#cccccc"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "#cccccc"))+
  geom_text(data=SalmonCollections %>% 
              mutate(Species = factor(Species,levels=c("Chum","Coho","Sockeye","Pink")))%>%
              filter(Species == "Pink") %>% slice_head(n=1),aes(x=Long+0.0009,y=Lat+0.00125,label=Label),size=3,fontface = "italic",angle=36)

#Combine Plots
p1<- pChum + theme(
  panel.background = element_rect(fill = "white", color = NA), # The plot area
  plot.background = element_rect(fill = "white", color = NA),
  axis.title = element_blank(),
  plot.margin = margin(10, 0, 10, 0, "pt"))
p2 <- pAuke +theme(legend.position = "none",
                   axis.title = element_blank(),
                   plot.margin = margin(10, 0, 10, 0, "pt"))+ 
  guides(fill = "none")
p3<- pNerka+theme(legend.position = "none",
                  axis.title = element_blank(),
                  plot.margin = margin(0, 0, 0, 0, "pt"))
p4<- pSkilak +theme(legend.position = "none",
                    axis.title = element_blank(),
                    plot.margin = margin(0, 0, 0, 0, "pt"),
                    panel.background = element_rect(fill = "white", color = NA), # The plot area
                    plot.background = element_rect(fill = "white", color = NA))

PTdf <- SalmonCollections%>%
  mutate(DateTmp = as.Date(paste(`Peak run (est)`,"/24",sep=""),format="%m/%d/%y"),
         PT = format(DateTmp, "%j"),
         PTnum = as.numeric(as.character(PT)))


### Try code Jan as 280 + JD

PTdf <- PTdf %>%
  mutate(PTnumPAD = ifelse(PTnum < 100, PTnum+280,PTnum))%>%
  dplyr::filter(!grepl(pattern="Auke Creek \\(odd",Location))


Prt<-PTdf %>%
  ggplot(.,aes(y=RunTimeFct,x=PTnumPAD,fill=Species))+
  geom_point(pch=21,size=3,position = position_jitter(width = 0, height = 0.15, seed = 123),aes(alpha=0.5))+
  xlab("Day of Year")+
  ylab("Run Timing")+
  ggrepel::geom_text_repel(data=PTdf %>%
              mutate(Species = factor(Species,levels=c("Chum","Coho","Sockeye","Pink"))),aes(x=PTnumPAD,y=RunTimeFct,label=Label),
            size=3,
            fontface = "italic",
            angle=45,
            hjust = 0,
            vjust = -1.2,
            segment.color = NA)+
  scale_fill_brewer(palette = 'Dark2')+
  geom_rect(aes(xmin = 183, xmax = 213, ymin = 0.2, ymax = 0.5),fill="grey50")+#Jul
  annotate("text",y=0.35,x= 198, colour="black",label="July")+
  geom_rect(aes(xmin = 214, xmax = 244, ymin = 0.2, ymax = 0.5),fill="grey50")+#Jul
  annotate("text",y=0.35,x= 229, colour="black",label="August")+
  geom_rect(aes(xmin = 245, xmax = 274, ymin = 0.2, ymax = 0.5),fill="grey50")+#Jul
  annotate("text",y=0.35,x= 259.5, colour="black",label="September")+
  coord_fixed(ratio=10)+
  annotate("segment", x = c(279, 280), xend = c(280, 281), y = c(-.25, -.25), yend = c(.25, .25))+
  annotate("segment", x = c(180, 281), xend = c(279, 281+30), y = c(0.05, 0.05), yend = c(0.05, 0.05))+
  coord_cartesian(clip = "off", ylim = c(0.5, 2)) +
  theme_classic()+
  theme(axis.line.x = element_blank())+
  geom_rect(aes(xmin = 282, xmax = 282+30, ymin = 0.2, ymax = 0.5),fill="grey50")+#Jul
  annotate("text",y=0.35,x= 295, colour="black",label="January")+
  scale_x_continuous(breaks = c(190,230,270,281+15),
                     labels = c("190","230","270","15"))

p5<-Prt +theme(legend.position = "none",
               plot.margin = margin(10, 0, 10, 0, "pt"))+ guides(fill = "none")

my_settings <- list(
  scale_fill_brewer(palette = 'Dark2'),
  guides(size = "none"),
  theme(plot.tag = element_text(size = 14, face = "bold")) # Hardcode the tag size here
)

p1 <- p1 + my_settings
p2 <- p2 + my_settings
p3 <- p3 + my_settings
p4 <- p4 + my_settings
p5 <- p5 + my_settings

pdf("./figures/FigureS1.pdf",width=12,height=12)
(p1 + p2 + p3 + p4)/p5 + plot_layout(guides = "collect",widths = c(6, 3),heights = c(6, 1))+
  plot_annotation(tag_levels = 'A') &
  scale_fill_brewer(palette = 'Dark2') &
  guides(size='none') 
dev.off()
