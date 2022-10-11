#Buffered land polygon

#load libraries ------------
library(R.matlab)
library(dplyr)
library(raster)
library(sf)
library(tidyr)
library(rnaturalearth)
library(ggplot2)

sf_use_s2 = FALSE

#projections --------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=19N +datum=WGS84 +units=m +no_defs"

#load coordinates for the sites ----------
eelgrass_coords <- read.csv("Data/EelgrassCoords.csv")%>%
  st_as_sf(coords=c("Long","Lat"),crs=latlong,remove=FALSE)

#read in the current range for extractions ----------- # going for a larger exent for plotting and the new raster-based analysis
eelgrass_region <- eelgrass_coords%>%
  filter(Region != "Pacific")%>%
  st_bbox()

eelgrass_region[1] <- -87 #slightly larger than the extent of the rasters from BNAM
eelgrass_region[2] <- 30
eelgrass_region[3] <- -52 
eelgrass_region[4] <- 68

eelgrass_region <- eelgrass_region%>%
  st_as_sfc()%>%
  st_as_sf()

plot_extent <- eelgrass_region%>%st_bbox()

plot_extent[1] <- -85 #for plotting make it a bit smaller than this region
plot_extent[2] <- 41
plot_extent[3] <- -52 
plot_extent[4] <- 64.5


#set up a dynamic azimutal equidistance projection for the study regions
#https://gis.stackexchange.com/questions/121489/1km-circles-around-lat-long-points-in-many-places-in-world/121539#121539

region_centre <- st_centroid(eelgrass_region)

aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                st_coordinates(region_centre)[2], st_coordinates(region_centre)[1])

#read in basemap to the large-scale extent (see environmental_extractions_largescale.R)
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
  st_transform(latlong)%>%
  st_intersection(.,eelgrass_region)%>%
  st_transform(aeqd) #equal area projection

#do a buffered using the new equal areas projection
basemap_buffered <- basemap%>%st_buffer(20*1000)%>%st_simplify() #20km buffer

#the eelgrass 'zone' is the sea within 20km of land (buffered - basemap)
eelgrass_zone <- basemap_buffered%>%st_difference(.,basemap)

#write the shapefile
st_write(obj = eelgrass_zone%>%st_transform(latlong),dsn = "Data/eelgrass_zone.shp")

#plot it ** note the extent is larger than the plot which is trimmed using the plot extent parameters. 
ggplot()+
  geom_sf(data=basemap%>%st_transform(latlong),col="black",fill=NA)+
  geom_sf(data=eelgrass_zone%>%st_transform(latlong),fill="darkgreen")+
  geom_sf(data=eelgrass_coords%>%filter(Region != "Pacific"),size=2)+
  theme_bw()+
  coord_sf(expand=0,xlim=c(plot_extent[c(1,3)]),ylim=c(plot_extent[c(2,4)]))


