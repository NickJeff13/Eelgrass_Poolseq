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

#read in the current range for extractions ----------- # going for a larger exent for plotting and the new raster-based analysis
eelgrass_region <- eelgrass_coords%>%
                   filter(Region != "Pacific")%>%
                   st_bbox()

eelgrass_region[1] <- -85
eelgrass_region[2] <- 41
eelgrass_region[3] <- -52 
eelgrass_region[4] <- 64.5

eelgrass_region <- eelgrass_region%>%
                   st_as_sfc()%>%
                   st_as_sf()

#generate seasonal and annual rasters ------------

  #Present data
   bnam_present <- c("Data/bnam/PRESENT_CLIMATE/SSS_PC.mat",
                    "Data/bnam/PRESENT_CLIMATE/Sbtm_PC.mat",
                    "Data/bnam/PRESENT_CLIMATE/SST_PC.mat",
                    "Data/bnam/PRESENT_CLIMATE/Tbtm_PC.mat")
  
         for(i in bnam_present){
             bnam_raster_gen_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = min,funname="min",output_dir = "output/bnam_rasters/large_scale/present_climate/")
             bnam_raster_gen_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = max,funname="max",output_dir = "output/bnam_rasters/large_scale/present_climate/")
             bnam_raster_gen_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = mean,funname="mean",output_dir = "output/bnam_rasters/large_scale/present_climate/")
         }
  
  #projections RCP 4.5 and 8.5 for 2075
   bnam_future_2075 <- c("Data/bnam/2075/RCP 4.5/dSST_dSSS_F_R45_2066-2085.mat",
                          "Data/bnam/2075/RCP 4.5/dTbtm_dSbtm_F_R45_2066-2085.mat",
                          "Data/bnam/2075/RCP 8.5/dSST_dSSS_F_R85_2066-2085.mat",
                          "Data/bnam/2075/RCP 8.5/dTbtm_dSbtm_F_R85_2066-2085.mat")
  
         for(i in bnam_future_2075){
  
           bnam_raster_gen_rcp_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = min,funname="min",output_dir = "output/bnam_rasters/large_scale/future_2075/")
           bnam_raster_gen_rcp_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = max,funname="max",output_dir = "output/bnam_rasters/large_scale/future_2075/")
           bnam_raster_gen_rcp_season(x=i,region=eelgrass_region,crs_val = latlong,rasterfun = mean,funname="mean",output_dir = "output/bnam_rasters/large_scale/future_2075/")
           
         }
   
#Get the projected conditions in the future (present + anomaly) ---------------
   
   #get the files lists for current (cf) and future (ff)
   
   cf_dir <- "output/bnam_rasters/large_scale/present_climate/"
   ff_dir <- "output/bnam_rasters/large_scale/future_2075/"
   
   cf <- dir(cf_dir,full.names = T)%>%
     data.frame(file=.)%>%
     mutate(name=gsub(".tif","",file),
            name=gsub(cf_dir,"",name),
            name=gsub("_PC","",name))%>%
     arrange(name)%>%
     rename(file_cf = file)
   
   #file names (full paths) of the large scale 2075 rasters
   ff <- dir(ff_dir,full.names = T)
   
   #RCP 4.5
   ff_45 <- ff[grepl("rcp45",ff)]%>%
     data.frame(file=.)%>%
     mutate(name=gsub("_rcp45.tif","",file),
            name=gsub(ff_dir,"",name))%>%
     arrange(name)%>%
     rename(file_45 = file)
   
   #RCP 8.5
   ff_85 <- ff[grepl("rcp85",ff)]%>%
     data.frame(file=.)%>%
     mutate(name=gsub("_rcp85.tif","",file),
            name=gsub(ff_dir,"",name))%>%
     arrange(name)%>%
     rename(file_85 = file)
   
   #merge it together in a master file match
   master_files <- cf%>%
     left_join(.,ff_45)%>%
     left_join(.,ff_85)
   
   #create directory to store the combined rasters
   future_condition_dir <- "output/bnam_rasters/large_scale/future_conditions/"
   if(!dir.exists(future_condition_dir )){dir.create(future_condition_dir)}
   
   #Function to combine the rasters --
   raster_combined <- function(x){
     
     #directory where the rasters will be written to
     out_dir <- "output/bnam_rasters/large_scale/future_conditions/"
     
     #read in the rasters (this is done on a row by row basis)
     ras_current <- raster(x$file_cf)
     ras_45 <- raster(x$file_45)
     ras_85 <- raster(x$file_85)
     
     #sum the rasters
     out_45 <- sum(ras_current,ras_45)
     out_85 <- sum(ras_current,ras_85)
     
     #write the rasters
     writeRaster(out_45,filename = paste0(out_dir,x$name,"_rcp45",".tif"),overwrite=TRUE)
     writeRaster(out_85,filename = paste0(out_dir,x$name,"_rcp85",".tif"),overwrite=TRUE)
     
   }
   
#Run the raster summation. Note this is a for loop but it is as fast as apply because it is writing each step. Added a progress message
   
   for(i in 1:nrow(master_files)){
     
     message(paste0("working on ",i," of ",nrow(master_files)))
     
     raster_combined(master_files[i,])
     
   }
   
   
#check to make sure the math is worked
   
   i = sample(1:60,1) #get a random file segment and check
   
   message("Running a raster compilation check on ",master_files[i,"name"],".")
   
   ras1 <- raster(master_files[i,"file_cf"])
   ras2 <- raster(master_files[i,"file_45"])
   ras3 <- raster(master_files[i,"file_85"])
   ras4 <- raster(paste0("output/bnam_rasters/large_scale/future_conditions/",master_files[i,"name"],"_rcp45.tif"))
   ras5 <- raster(paste0("output/bnam_rasters/large_scale/future_conditions/",master_files[i,"name"],"_rcp85.tif"))
   
   focuspoint <- data.frame(lon= -78.29496, lat= 59.14498)%>%
     st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
     as_Spatial()
   
   test1 <- round(raster::extract(ras1,focuspoint) + raster::extract(ras2,focuspoint),5) == round(raster::extract(ras4,focuspoint),5);test1
   test2 <- round(raster::extract(ras1,focuspoint) + raster::extract(ras3,focuspoint),5) == round(raster::extract(ras5,focuspoint),5);test2
   
   if(test1 & test2){message("Raster validation complete. It worked.")}
   

  
  