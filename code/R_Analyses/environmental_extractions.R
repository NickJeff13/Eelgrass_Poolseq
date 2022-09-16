#load libraries ------------
library(R.matlab)
library(dplyr)
library(raster)
library(sf)
library(tidyr)

sf_use_s2 = FALSE

#source functions ----------
source("code/R_Analyses/raster_gen.R")
source("code/R_Analyses/coord_bump.R")

#projections --------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=19N +datum=WGS84 +units=m +no_defs"

#load coordinates for the sites ----------
eelgrass_coords <- read.csv("Data/EelgrassCoords.csv")%>%
                   st_as_sf(coords=c("Long","Lat"),crs=latlong,remove=FALSE)

#read in the current range for extractions -----------
eelgrass_region <- eelgrass_coords%>%
                   filter(Region != "Pacific")%>%
                   st_bbox()%>%
                   st_as_sfc()%>%
                   st_as_sf()%>%
                   st_buffer(3) # three degree buffer

#generate seasonal and annual rasters ------------

  #Present data
  
  # bnam_present <- c("Data/bnam/PRESENT_CLIMATE/SSS_PC.mat",
  #                  "Data/bnam/PRESENT_CLIMATE/Sbtm_PC.mat",
  #                  "Data/bnam/PRESENT_CLIMATE/SST_PC.mat",
  #                  "Data/bnam/PRESENT_CLIMATE/Tbtm_PC.mat")
  # 
  #       for(i in bnam_present){
  #           bnam_raster_gen_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = min,funname="min",output_dir = "output/bnam_rasters/present_climate/")
  #           bnam_raster_gen_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = max,funname="max",output_dir = "output/bnam_rasters/present_climate/")
  #           bnam_raster_gen_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = mean,funname="mean",output_dir = "output/bnam_rasters/present_climate/")
  #       }
  # 
  # #projections RCP 4.5 and 8.5 for 2075
  # 
  # bnam_future_2075 <- c("Data/bnam/2075/RCP 4.5/dSST_dSSS_F_R45_2066-2085.mat",
  #                        "Data/bnam/2075/RCP 4.5/dTbtm_dSbtm_F_R45_2066-2085.mat",
  #                        "Data/bnam/2075/RCP 8.5/dSST_dSSS_F_R85_2066-2085.mat",
  #                        "Data/bnam/2075/RCP 8.5/dTbtm_dSbtm_F_R85_2066-2085.mat")
  # 
  #       for(i in bnam_future_2075){
  #         
  #         bnam_raster_gen_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = min,funname="min",output_dir = "output/bnam_rasters/future_2075/")
  #         bnam_raster_gen_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = max,funname="max",output_dir = "output/bnam_rasters/future_2075/")
  #         bnam_raster_gen_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = mean,funname="mean",output_dir = "output/bnam_rasters/future_2075/")
  #         
  #       }
  # 
  # #projections RCP 4.5 and 8.5 for 2055
  # 
  # bnam_future_2055 <- c("Data/bnam/2055/RCP 4.5/dSST_dSSS_F_R45_2046-2065.mat",
  #                      "Data/bnam/2055/RCP 4.5/dTbtm_dSbtm_F_R45_2046-2065.mat",
  #                      "Data/bnam/2055/RCP 8.5/dSST_dSSS_F_R85_2046-2065.mat",
  #                      "Data/bnam/2055/RCP 8.5/dTbtm_dSbtm_F_R85_2046-2065.mat")
  # 
  #     for(i in bnam_future_2055){
  #       
  #       bnam_raster_gen_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = min,funname="min",output_dir = "output/bnam_rasters/future_2055/")
  #       bnam_raster_gen_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = max,funname="max",output_dir = "output/bnam_rasters/future_2055/")
  #       bnam_raster_gen_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = mean,funname="mean",output_dir = "output/bnam_rasters/future_2055/")
  #       
  #     }
  
  
#extract the values for each site ---------------
  
  
  #present climate extractions --
  
  #get the coordinates
  eelgrass_present_coords <- eelgrass_coords%>%
                             filter(Region != "Pacific",!Code == "Ebay")%>% # no values for the Bras d'Or lakes and though they can be bumped to the adjacent ocean that might not work for the Bras d'or given how unique they are
                             st_transform("+proj=longlat +datum=NAD83 +no_defs")%>% #bnam projection
                             data.frame()%>%
                             rename(lat=Lat,long=Long)%>%
                             dplyr::select(lat,long)
  
  #do a check of the coordinates as some are still on land. 
  test_present_ras <- raster("output/bnam_rasters/present_climate/annual_max_Sbtm.mon.tif")
    
  #run check to adjust coordinates
  bump_check <- NULL
  for(i in 1:nrow(eelgrass_present_coords)){bump_check <- rbind(bump_check,coord_bump(eelgrass_present_coords[i,],test_present_ras,radius = 20))}

  eelgrass_present_coords$long_a <- bump_check$long
  eelgrass_present_coords$lat_a <- bump_check$lat
  
  #secondary check because the present salinity seems to have a slightly different coverage that requires a different bump
  test_ras_2 <- raster("output/bnam_rasters/present_climate/annual_min_SSS.mon.tif")
  
  eelgrass_present_coords2 <- eelgrass_present_coords[,c("long_a","lat_a")]%>%rename(long=long_a,lat=lat_a)
  
  #run check to adjust coordinates - need to do this a second time since the salinity present rasters seem to be slightly different resolutions
  bump_check2 <- NULL
  for(i in 1:nrow(eelgrass_present_coords)){bump_check2 <- rbind(bump_check2,coord_bump(eelgrass_present_coords2[i,],test_ras_2,radius = 20))}
  
  eelgrass_present_coords$long_a <- bump_check2$long
  eelgrass_present_coords$lat_a <- bump_check2$lat
  
  #manually bump NRIV which seems to cause a few issues
  eelgrass_present_coords[which(eelgrass_present_coords$lat==46.2969),"lat_a"] <- 46.363760
  eelgrass_present_coords[which(eelgrass_present_coords$lat==46.2969),"long_a"] <- -60.467881
  
  #spatial 'sp' with adjusted coordinates
  eelgrass_sp <- eelgrass_present_coords%>%
    st_as_sf(coords=c("long_a","lat_a"),crs=proj4string(test_present_ras))%>%
    as_Spatial()
    
    # plot to view the bumped coordinates
    # ind <- eelgrass_present_coords%>%
    #        st_as_sf(coords=c("long","lat"),crs=proj4string(test_present_ras))%>%
    #        as_Spatial()%>%
    #        extract(test_present_ras,.)%>%
    #        is.na(.)
    # 
    # demo_df <- eelgrass_present_coords%>%
    #            filter(ind)%>%
    #            st_as_sf(coords=c("long","lat"),crs=proj4string(test_present_ras))%>%
    #            as_Spatial()
    # 
    # demo_df_adj <- eelgrass_present_coords%>%
    #                 filter(ind)%>%
    #                 st_as_sf(coords=c("long_a","lat_a"),crs=proj4string(test_present_ras))%>%
    #                 as_Spatial()
    # 
    # plot(test_present_ras)
    # points(demo_df,pch=19,col="red")
    # points(demo_df_adj,pch=19,col="blue")

  
 #create data.frames that can be used to assemble the data later 
  present_df <- data.frame(filename=dir("output/bnam_rasters/present_climate/"))%>%
                rowwise()%>%
                mutate(content=gsub(".mon.tif","",filename),
                       season = strsplit(content,"_",content)%>%unlist()%>%.[1],
                       variable = strsplit(content,"_",content)%>%unlist()%>%.[2],
                       parameter = strsplit(content,"_",content)%>%unlist()%>%.[3],
                       rcp = NA,
                       period="present",
                       rasterpath = paste0("output/bnam_rasters/present_climate/",filename))%>%
                data.frame()%>%
                dplyr::select("filename","content","season","variable","period","rcp","parameter","rasterpath")
  
  future_2055_df <- data.frame(filename=dir("output/bnam_rasters/future_2055/"))%>%
                    rowwise()%>%
                    mutate(content=gsub(".tif","",filename),
                           season = strsplit(content,"_",content)%>%unlist()%>%.[1],
                           variable = strsplit(content,"_",content)%>%unlist()%>%.[2],
                           content = strsplit(content,"_",content)%>%unlist()%>%.[3],
                           period = strsplit(content,"\\.",content)%>%unlist()%>%.[2],
                           rcp = strsplit(content,"\\.",content)%>%unlist()%>%.[1],
                           parameter = strsplit(content,"\\.",content)%>%unlist()%>%.[3],
                           content = paste(season,variable,gsub("\\.","_",content),sep="_"),
                           rasterpath = paste0("output/bnam_rasters/future_2055/",filename))%>%
                    data.frame%>%
                    dplyr::select("filename","content","season","variable","period","rcp","parameter","rasterpath")
  
  future_2075_df <- data.frame(filename=dir("output/bnam_rasters/future_2075/"))%>%
                    rowwise()%>%
                    mutate(content=gsub(".tif","",filename),
                           season = strsplit(content,"_",content)%>%unlist()%>%.[1],
                           variable = strsplit(content,"_",content)%>%unlist()%>%.[2],
                           content = strsplit(content,"_",content)%>%unlist()%>%.[3],
                           period = strsplit(content,"\\.",content)%>%unlist()%>%.[2],
                           rcp = strsplit(content,"\\.",content)%>%unlist()%>%.[1],
                           parameter = strsplit(content,"\\.",content)%>%unlist()%>%.[3],
                           content = paste(season,variable,gsub("\\.","_",content),sep="_"),
                           rasterpath = paste0("output/bnam_rasters/future_2075/",filename))%>%
                    data.frame%>%
                    dplyr::select("filename","content","season","variable","period","rcp","parameter","rasterpath")
  
  bnam_raster_df <- rbind(present_df,future_2055_df,future_2075_df)
  
  #start a data.frame that will be used to associated extractions with names 'code'
  eelgrass_extracts <- eelgrass_present_coords%>%
                      dplyr::select(lat)%>%
                      left_join(.,eelgrass_coords%>%data.frame()%>%dplyr::select(Lat,Code)%>%rename(lat=Lat))%>%
                      dplyr::select(Code)
                      
  #do the extractions ~~ takes < 5 min ------------- 
  for(i in 1:nrow(bnam_raster_df)){
    
    #crude progress bar for tracking errors 
    message(paste0(i," of ",nrow(bnam_raster_df)))
    
    #extract raster
    ras <- raster(bnam_raster_df[i,"rasterpath"])
    
    #values for each new variable assigned as columns from extracts
    eelgrass_extracts[[bnam_raster_df[i,"content"]]] <- with(eelgrass_extracts,raster::extract(ras,eelgrass_sp))
  }
  
#final data.frame of extracts form ------
  bnam_extracts <- eelgrass_extracts%>%
                    gather(content,value,-Code)%>%
                    left_join(.,bnam_raster_df%>%dplyr::select(-c(filename,rasterpath)))%>%
                    dplyr::select(Code,rcp,period,season,variable,parameter,value)

  
  #format from anomalies to real values
  bnam_formatted <- bnam_extracts 
  
  parameters <- c("Sbtm","Tbtm","SSS","SST")
  
  for(i in unique(bnam_extracts$Code)){
      for(s in unique(bnam_extracts$season)){
        for(v in c("max","min","mean")){
          for(p in parameters){
            
            ind <- which(bnam_formatted$Code == i &  #index for the values (1 for the present day, 2 for each of the RCP values for 2055 and 2075)
                         bnam_formatted$season == s &
                         bnam_formatted$variable == v &
                         grepl(p,bnam_formatted$parameter))
          
            #the first indexed value is the present day value which has an 'NA' as the rcp.  The check will warn you if this isn't the case
            if(which(is.na(bnam_formatted[ind,"rcp"])) != 1){message(paste("Something wrong with indexing",i,s,v,p,sep=" "))}
            
            bnam_formatted[ind[2:5],"value"] <- bnam_formatted[ind[1],"value"]+bnam_formatted[ind[2:5],"value"]
          
          }#parameter loop
        }#variable loop
      }#Season loop
   } #Code loop
  
#now clean up the naming because the 'd' is not required as these are no long anomalies. 
  bnam_formatted$parameter <- gsub("d","",bnam_formatted$parameter)
  
#save outputs ------------------
  save(bnam_extracts,file="output/bnam_extracts_anomalies.RData")
  save(bnam_formatted,file="output/bnam_extracts_values.RData")
  
  eelgrass_coords_formatted <- eelgrass_present_coords%>%left_join(.,eelgrass_coords%>%data.frame()%>%dplyr::select(Lat,Code)%>%rename(lat=Lat))
  write.csv(eelgrass_coords_formatted,"output/bnam_extracts_bumped_coords.csv",row.names=FALSE)
  
  