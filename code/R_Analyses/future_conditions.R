#Code to create 'future' layers from BNAM

#Load libraries
library(dplyr)
library(tidyr)
library(raster)
library(sf)
library(sp)
library(R.matlab)

#set projection to be used with the BNAM files
latlong <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#load in extent of focus (Maritimes bioreigon)
region <- read_sf("data/MaritimesPlanningArea.shp")%>%st_transform(latlong)

#load functions for creating the raster geotiffs
source("code/raster_gen.R")

#Generate the rasters extractions from BNAM (only need to do this once!). i have commented out because the rasters on on github
 
   matlab_files <- dir("data/BNAM data/",full.names = TRUE)
  
  
    # 
    # for(f in anom_files){
    #   
    #   bnam_raster_gen(f,region=region,crs_val=latlong)
    #   
    # }

#Combine the anomalies and the contemporary rasters to create a 'future' conditions

   #RCP 4.5 
   #Sea Surface Salinity
   sum_rasters(file1 = "output/bnam_rasters/RCP45_2075_dSSS_ann.tif",
               file2 = "output/bnam_rasters/SSS_ann.tif",
               outputfile = "output/future_climate/RCP45_2075_future_SSS.tif")
   
   #Sea Surface Temperature
   sum_rasters(file1 = "output/bnam_rasters/RCP45_2075_dSST_ann.tif",
               file2 = "output/bnam_rasters/SST_ann.tif",
               outputfile = "output/future_climate/RCP45_2075_future_SST.tif")
   
   #Bottom Salinity
   sum_rasters(file1 = "output/bnam_rasters/RCP45_2075_dSbtm_ann.tif",
               file2 = "output/bnam_rasters/Sbtm_ann.tif",
               outputfile = "output/future_climate/RCP45_2075_future_BS.tif")
   
   #Bottom Temperature
   sum_rasters(file1 = "output/bnam_rasters/RCP45_2075_dTbtm_ann.tif",
               file2 = "output/bnam_rasters/Tbtm_ann.tif",
               outputfile = "output/future_climate/RCP45_2075_future_BT.tif")
   
   #RCP 8.5 
   #Sea Surface Salinity
   sum_rasters(file1 = "output/bnam_rasters/RCP85_2075_dSSS_ann.tif",
               file2 = "output/bnam_rasters/SSS_ann.tif",
               outputfile = "output/future_climate/RCP85_2075_future_SSS.tif")
   
   #Sea Surface Temperature
   sum_rasters(file1 = "output/bnam_rasters/RCP85_2075_dSST_ann.tif",
               file2 = "output/bnam_rasters/SST_ann.tif",
               outputfile = "output/future_climate/RCP85_2075_future_SST.tif")
   
   #Bottom Salinity
   sum_rasters(file1 = "output/bnam_rasters/RCP85_2075_dSbtm_ann.tif",
               file2 = "output/bnam_rasters/Sbtm_ann.tif",
               outputfile = "output/future_climate/RCP85_2075_future_BS.tif")
   
   #Bottom Temperature
   sum_rasters(file1 = "output/bnam_rasters/RCP85_2075_dTbtm_ann.tif",
               file2 = "output/bnam_rasters/Tbtm_ann.tif",
               outputfile = "output/future_climate/RCP85_2075_future_BT.tif")

#do a check to make sure it worked
  sample_station <- c(44,-62.5)
  
  ras1 <- raster("output/bnam_rasters/SST_ann.tif") #now
  ras2 <- raster("output/bnam_rasters/RCP85_2075_dSST_ann.tif") #anomaly
  ras3 <- raster("output/future_climate/RCP85_2075_future_SST.tif") #future
  
  now <- extract(ras1,sample_station)[1]
  future <- extract(ras3,sample_station)[1]
  difference <- future-now
  
  #if the difference is the same as the anomaly the addition worked (TRUE)
  round(difference,4) == round(extract(ras2,sample_station)[1],4) #there is some trimming of significant figures so better to round


