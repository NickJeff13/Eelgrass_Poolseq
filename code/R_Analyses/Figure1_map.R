##remake Figure 1 

#load libraries
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(viridis)
library(tidyr)

sf_use_s2(FALSE) #makes some of the sf functions run a bit less buggy

#Projections ------------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#colours for the points
colour_df <- data.frame(code=c("NAH","GRB","PORT","HEB","PRJ","SAC","MASI","SAM","L3F","TAYH","POUL","EBAY","NRIV",
                               "MELM","SUM","PETI","POK","RIM","SEPT","BUCK","JB38","JB33","TSW"),
                        col=viridis(n=23))

colour_df <- data.frame(code=c("NAH","GRB","PORT","HEB","PRJ","SAC","MASI","SAM","L3F","TAYH","POUL","EBAY","NRIV",
                                "MELM","SUM","PETI","POK","RIM","SEPT","BUCK","JB38","JB33","TSW"),
                        old_code = c("NAH","GRB","PORT","HEB","PRJ","SAC","MASI","SAM","L3F","TH","POUL","Ebay","NRIV",
                                     "MELM","SUM","PETITE","POK","RIM","SEPT","BUCK","CH38","CH33","TSW"),
                        col=c("dodgerblue2", "#E31A1C","green4","#6A3D9A","#FF7F00","black", "gold1","#FB9A99","#CAB2D6",
                               "#FDBF6F","gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
                                "darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown" )) #from https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
#load the coordinates
eelgrass_coords <- read.csv("Data/allsites_coords.csv")%>%
  st_as_sf(coords=c("long","lat"),crs=latlong,remove=FALSE)%>%
  st_transform(CanProj)%>%
  mutate(code=case_when(code=="TH"~"TAYH", #fix these to match the plot
                        code=="Ebay"~"EBAY",
                        code=="PETITE"~"PETI",
                        code=="CH33"~"JB33",
                        code=="CH38"~"JB38",
                        TRUE ~ code))%>%
  left_join(.,colour_df)


#save the colour dataframe for other plots
write.csv(colour_df,"Data/Colour_codes.csv",row.names=FALSE)  

#Shapefiles
can_eez <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/Canada_EEZ.shp")%>%st_transform(CanProj)

ns_coast <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/Coastline/NS_coastline_project_Erase1.shp")%>%st_transform(CanProj)

#create a map boundary for the analysis based on the Canadian EEZ
can_bounding <- can_eez%>%
              st_transform(latlong)%>%
              st_buffer(2)%>% #2 degree buffer on the Canadian EEZ extent
              st_transform(CanProj)%>%
              st_bbox() #bounding box used for plot limits #will give some warnings but it doesn't distort the bounding box by much and it's only a bounding box anyway

point_bound <- rbind(
  
                eelgrass_coords%>%
                  st_bbox()%>%
                  st_as_sfc()%>%
                  st_as_sf()%>%
                  st_buffer(25*1000),
                
                can_eez%>%
                  st_transform(CanProj)%>%
                  st_bbox()%>%
                  st_as_sfc()%>%
                  st_as_sf()
                
                  )%>%
  
                st_bbox()%>%
                st_as_sfc()%>%
                st_as_sf()%>%
                st_buffer(200*1000)%>%
                st_bbox()%>%
                st_as_sfc()%>%
                st_as_sf()
                
                          

#create a basemap for plotting MPAs
basemap <- rbind(ne_states(country = "Canada",returnclass = "sf")%>%
                   dplyr::select(name_en,geometry)%>%
                   st_as_sf()%>%
                   st_union()%>%
                   st_transform(latlong)%>%
                   st_as_sf()%>%
                   mutate(country="Canada"),
                 ne_states(country = "United States of America",returnclass = "sf")%>%
                   dplyr::select(name_en,geometry)%>%
                   st_as_sf()%>%
                   st_union()%>%
                   st_transform(latlong)%>%
                   st_as_sf()%>%
                   mutate(country="USA"),
                 ne_states(country = "Greenland",returnclass = "sf")%>%
                   dplyr::select(name_en,geometry)%>%
                   st_as_sf()%>%
                   st_union()%>%
                   st_transform(latlong)%>%
                   st_as_sf()%>%
                   mutate(country="USA"))%>%
  st_union()%>%
  st_as_sf()%>%
  st_transform(CanProj)%>%
  st_intersection(.,point_bound)

Canada <-  ne_states(country = "Canada",returnclass = "sf")%>%
  dplyr::select(name_en,geometry)%>%
  st_as_sf()%>%
  st_union()%>%
  st_transform(CanProj)%>%
  st_as_sf()%>%
  mutate(country="Canada")

Atlantic_bound <- eelgrass_coords%>%
                  filter(grouping=="atlantic",code!="BUCK")%>%
                  st_bbox()%>%
                  st_as_sfc()%>%
                  st_as_sf()%>%
                  st_buffer(100*1000)%>%
                  st_bbox()
                  

Atlantic_bound

#Atlantic map

atlantic <- ggplot()+
  geom_sf(data=basemap,lwd=0,col=NA,fill="grey60")+
  geom_sf(data=Canada,lwd=0,col=NA,fill="grey50")+
  geom_sf(data=eelgrass_coords,aes(fill=col),shape=21,size=4)+
  coord_sf(expand=0,xlim=Atlantic_bound%>%st_bbox()%>%.[c(1,3)],ylim=Atlantic_bound%>%st_bbox()%>%.[c(2,4)])+
  theme_bw()+
  annotation_scale(location="br")+
  annotation_north_arrow(location="tl",height=unit(0.75, "cm"),width=unit(0.75, "cm"))+
  scale_fill_identity()+
  theme(legend.position = "none")

#Inset map
inset <- ggplot()+
  geom_sf(data=basemap,lwd=0,col=NA,fill="grey60")+
  geom_sf(data=Canada,lwd=0,col=NA,fill="grey50")+
  geom_sf(data=eelgrass_coords,aes(fill=col),shape=21)+
  geom_sf(data=Atlantic_bound%>%st_as_sfc()%>%st_as_sf(),fill=NA,col="black")+
  coord_sf(expand=0,xlim=point_bound%>%st_bbox()%>%.[c(1,3)],ylim=point_bound%>%st_bbox()%>%.[c(2,4)])+
  theme_bw()+
  scale_fill_identity()+
  theme(legend.position = "none")

ggsave("Figures/atlantic_point_map.png",atlantic,height=8,width=6,units="in",dpi=300)
ggsave("Figures/eelgrass_point_map_inset.png",inset,height=8,width=8,units="in",dpi=300)

#micro inset map

#MASI-SAC
masi_inset_bounds <- eelgrass_coords%>%
                     filter(code %in% c("MASI","SAC"))%>%
                     st_bbox()%>%
                     st_as_sfc()%>%
                     st_as_sf()%>%
                     st_buffer(5*1000)%>%
                     st_bbox()
masi_sac <- ggplot()+
            geom_sf(data=ns_coast,lwd=0,col=NA,fill="grey60")+
            geom_sf(data=eelgrass_coords%>%filter(code %in% c("MASI","SAC")),aes(fill=col),shape=21,size=4)+
            coord_sf(expand=0,xlim=masi_inset_bounds[c(1,3)],ylim=masi_inset_bounds[c(2,4)])+
            theme_bw()+
            scale_fill_identity()+
            theme(legend.position = "none")+
            annotation_scale(location="br")


#get high resolution coastline data for James bay - https://maps-cartes.ec.gc.ca/arcgis/rest/services/ShorelineMappingNWRC_CartographieDuLittoralCNRF/MapServer

jb_link <- "https://maps-cartes.ec.gc.ca/arcgis/rest/services/ShorelineMappingNWRC_CartographieDuLittoralCNRF/MapServer/2"
jb_link <- "https://maps-cartes.ec.gc.ca/arcgis/rest/services/ShorelineMappingNWRC_CartographieDuLittoralCNRF/MapServer/3"

jb_sf <- esri2sf::esri2sf(paste0(jb_link, "0"))
st_write(jb_sf,"Data/James_Bay.shp")

jb_inset_bounds <- eelgrass_coords%>%
                  filter(code %in% c("JB33","JB38"))%>%
                  st_bbox()%>%
                  st_as_sfc()%>%
                  st_as_sf()%>%
                  st_buffer(15*1000)%>%
                  st_bbox()

jb_inset <- ggplot()+
  geom_sf(data=Canada,lwd=0,col=NA,fill="grey60")+
  geom_sf(data=eelgrass_coords%>%filter(code %in% c("JB33","JB38")),aes(fill=col),shape=21,size=4)+
  coord_sf(expand=0,xlim=jb_inset_bounds[c(1,3)],ylim=jb_inset_bounds[c(2,4)])+
  theme_bw()+
  scale_fill_identity()+
  theme(legend.position = "none")+
  annotation_scale(location="br")
