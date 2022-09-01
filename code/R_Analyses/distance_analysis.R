## code for calculating distance metrics between sample sites of the eel grass genomics program 

#load libraries --------
library(rnaturalearth)
library(rnaturalearthhires)
library(sf)
library(tidyr)
library(dplyr)
library(raster)
library(ggplot2)
library(gdistance)
library(nngeo)
library(fasterize)
library(reshape2)
library(patchwork)


#Projections -----
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#Source functions for analysis ------------
source("code/coord_bump.R")
source("code/lcp_function.R")

##Load in sample site coordinates ------------
site_coords <- read.csv("Data/site_coordinates.csv")

##Bras d'Or sites cause issues, so we have to manually bump onto water or the coord_bump function (internal to lcp_function) will 'bump' it into the Bras d'Or where it becomes an issue
north_river_org <- site_coords%>%filter(code == "NRIV")
eastbay_org <- site_coords%>%filter(code=="Ebay")

site_coords[site_coords$code == "NRIV","lat"] <- 46.319114 # this is 12.65 km from North River (as the fish would swim) this will need to be added on later -- Google Earth Path
site_coords[site_coords$code == "NRIV","long"] <- -60.520537

site_coords[site_coords$code == "Ebay","lat"] <-  46.314814 # this is 92.56 km from North River (as the fish would swim) this will need to be added on later -- Google Earth Path
site_coords[site_coords$code == "Ebay","long"] <- -60.398887

NRIV_offset <- 12.65
Ebay_offset <- 92.56

site_coords <-site_coords%>%
              st_as_sf(coords=c("long","lat"),crs=latlong,remove=FALSE)%>%
              mutate(site_id = location) # to match code in lcp_function that was developed for the biodiversity_edna project


#create bounding box for the coordinates
site_bounds_all <- site_coords%>%
                   st_bbox()%>%
                   st_as_sfc()%>%
                   st_as_sf()

site_bounds_east <- site_coords%>%
                    filter(grouping != "pacific")%>%
                    st_bbox()%>%
                    st_as_sfc()%>%
                    st_as_sf()

site_bounds_atlantic <- site_coords%>%
                        filter(grouping == "atlantic")%>%
                        st_bbox()%>%
                        st_as_sfc()%>%
                        st_as_sf()

#Source basemap -------------- 
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
                   mutate(country="USA"))%>%st_union()%>%st_as_sf()%>%suppressMessages()

#Create a bounding box for the analysis. The way the LCP analysis works is that it will seek a way around a trimmed barrier. So if we 
#cut the land too short it will come around the artificial bottom and work its way around the backend (e.g., go from the GOF through the US and to James bay)
#so we manually create a bounding box here that can be applied int he code. 
bound_box_large <- st_bbox(c(xmin = -100, xmax = -45, ymax = 80, ymin = 25), crs = latlong)%>%
                   st_as_sfc()%>%
                   st_as_sf()

bound_box_small <- st_bbox(c(xmin = -75, xmax = -45, ymax = 57, ymin = 35), crs = latlong)%>%
  st_as_sfc()%>%
  st_as_sf()

      ggplot()+
         geom_sf(data=basemap)+
         geom_sf(data=bound_box_large,fill=NA)+
         geom_sf(data=bound_box_small,fill=NA)+
         geom_sf(data=site_coords)+
         theme_bw()+
         coord_sf(xlim=c(-105,-40),ylim=c(20,85))


#Do a large scale (coarse scale - 10 km) analysis for the largest extent 
large_scale_lcp <- lcp_function(x=site_coords%>%
                                  filter(grouping != "pacific"),
                                basemap=basemap,
                                bound_box = bound_box_large,
                                resolution = 10,
                                transition_name = "atlantic_arcitc",
                                rad=42,
                                lines=TRUE,
                                dirs=8,
                                matrix_out = TRUE,
                                recalculate = TRUE) #this will recalculate the transition layer which you don't need to do, but it does make sure everything is in line.

small_scale_lcp <- lcp_function(x=site_coords%>%
                                  filter(grouping == "atlantic"),
                                basemap=basemap,
                                bound_box = bound_box_small,
                                resolution = 5,
                                transition_name = "atlantic",
                                rad=25,
                                lines=TRUE,
                                dirs=8,
                                matrix_out = FALSE,
                                recalculate = TRUE) #this will recalculate the transition layer which you don't need to do, but it does make sure everything is in line. 

#plot results to make sure no major errors occurred (e.g., land skipping)
lines_ls <- large_scale_lcp[[2]]
lines_sm <- small_scale_lcp[[2]]

basemap_large <- basemap%>%st_intersection(bound_box_large)
basemap_small <- basemap%>%st_intersection(bound_box_small)

plotbounds_lg <- lines_ls%>%
                st_bbox()%>%
                st_as_sfc()%>%
                st_as_sf()%>%
                st_buffer(1)%>%
                st_bbox()

plotbounds_sm <- lines_sm%>%
                st_bbox()%>%
                st_as_sfc()%>%
                st_as_sf()%>%
                st_buffer(1)%>%
                st_bbox()

plot_large <- ggplot()+
            geom_sf(data=basemap_large)+
            geom_sf(data=site_coords%>%filter(grouping != "pacific"),size=2)+
            geom_sf(data=lines_ls%>%filter(len>1))+
            theme_bw()+
            coord_sf(xlim = plotbounds_lg[c(1,3)],ylim = plotbounds_lg[c(2,4)])

plot_small <- ggplot()+
              geom_sf(data=basemap_small)+
              geom_sf(data=site_coords%>%filter(grouping == "atlantic"),size=2)+
              geom_sf(data=lines_sm%>%filter(len>1))+
              theme_bw()+
              coord_sf(xlim = plotbounds_sm[c(1,3)],ylim = plotbounds_sm[c(2,4)])

least_cost_plots <- plot_large + plot_small + plot_layout(ncol=2)

least_cost_plots

ggsave("output/least_cost_paths_maps.png",least_cost_plots,width=7,units="in",dpi=600)

#Complete distance matrix --------

    #get the distance data
    dists_ls <- large_scale_lcp[[1]]
    dists_sm <- small_scale_lcp[[1]]
    
    #get the distance data and fix the offsets
    dists_ls <- large_scale_lcp[[1]]%>%
                rowwise()%>%
                mutate(dist_adj = ifelse(start == "North River" & end !="North River",dist + NRIV_offset,dist), #Bras d'Or offsets
                       dist_adj = ifelse(start == "East Bay" & end !="East Bay",dist_adj + Ebay_offset,dist_adj),
                       dist_adj = ifelse(end == "North River" & start !="North River",dist_adj + NRIV_offset,dist_adj),
                       dist_adj = ifelse(end == "East Bay" & start !="East Bay",dist_adj + Ebay_offset,dist_adj))%>%
                data.frame()%>%
                dplyr::select(start,end,dist_adj)%>%
                rename(dist=dist_adj)
    
    dists_sm <- small_scale_lcp[[1]]%>%
                rowwise()%>%
                mutate(dist_adj = ifelse(start == "North River" & end !="North River",dist + NRIV_offset,dist), #Bras d'Or offsets
                       dist_adj = ifelse(start == "East Bay" & end !="East Bay",dist_adj + Ebay_offset,dist_adj),
                       dist_adj = ifelse(end == "North River" & start !="North River",dist_adj + NRIV_offset,dist_adj),
                       dist_adj = ifelse(end == "East Bay" & start !="East Bay",dist_adj + Ebay_offset,dist_adj))%>%
                  data.frame()%>%
                  dplyr::select(start,end,dist_adj)%>%
                  rename(dist=dist_adj)
    
    #now the two James Bay sites need to be added into the smaller scaled analysis for the 'atlantic' sites and save as a distance matrix
    dists_combo <- rbind(dists_sm,dists_ls%>%filter(start %in% c("James Bay CH33","James Bay CH38") | end %in% c("James Bay CH33","James Bay CH38")))
    
    lcp_dist <- dists_combo%>%
                  pivot_wider(names_from = end,values_from = dist)%>%
                  dplyr::select(-start)%>%
                  as.matrix()
    
    row.names(lcp_dist) <-  colnames(lcp_dist) 
    
#Now get the 'as the bird would fly distances'    

    geo_dist <- site_coords%>%
               filter(grouping != "pacific")%>%
               mutate(location = factor(location,levels=rownames(lcp_dist)))%>% #so the output matches the matrix from the lcp analysis
               arrange(location)%>%
               st_distance()%>%
               as.matrix()/1000
    
    geo_dist <- geo_dist%>%units::drop_units()
    
    rownames(geo_dist) <- site_coords%>%
                          filter(grouping != "pacific")%>% # this is overkill because it should just be the factor levels - rownames(lcp_dist) -  but this makes sure. 
                          mutate(location = factor(location,levels=rownames(lcp_dist)))%>%
                          arrange(location)%>%
                          pull(location)
    
    colnames(geo_dist) <- site_coords%>%
                          filter(grouping != "pacific")%>%
                          mutate(location = factor(location,levels=rownames(lcp_dist)))%>%
                          arrange(location)%>%
                          pull(location)

#save the outputs
    save(lcp_dist,file="output/lcp_dissimilarity.RData")
    save(geo_dist,file="output/geographic_dissimilarity.RData")
    
