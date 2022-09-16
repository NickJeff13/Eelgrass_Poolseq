##First order adjustment function function 'coord_bump' ---------------------
coord_bump <- function(x,ras,radius = 10){
  
  require(dplyr)
  require(sf)
  require(raster)
  require(nngeo)
  
  # x = dataframe with a column for lat and long that has the same projection as ras
  # ras = bathymetry raster
  # radius = search radius in km (default is 10 or 10000 m)
 
  rasproj <- proj4string(ras)
  
  x <- st_as_sf(x,coords=c("long","lat"),crs=rasproj) #convert to sf
  x2 <- x%>%st_transform("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
  
  #set up a dynamic azimutal equidistance projection for each point 
  #https://gis.stackexchange.com/questions/121489/1km-circles-around-lat-long-points-in-many-places-in-world/121539#121539
  
    aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                  st_coordinates(x2)[2], st_coordinates(x2)[1])
  
  cr <- x%>%
        st_transform(aeqd)%>%
        st_buffer(radius*1000)%>% # km search radius for water
        st_as_sf()%>%
        st_transform(proj4string(ras))%>% #Back to the projection of the raster
        extent()

  bathy <- crop(ras,cr)
  
  #this is a logical to catch scenarios in the present climate where the rasters have 0 instead of NAs for land
  if(sum(is.na(values(bathy))) == 0){values(bathy)[values(bathy) == 0] <- NA}
  
  xy <- st_coordinates(x)%>%SpatialPoints()
  
  depth <- raster::extract(bathy,xy)

  if(is.na(depth)){
    
    if(sum(!is.na(values(bathy)))>1){
         bathypoints <- coordinates(bathy)[!is.na(values(bathy)),]%>%
                        data.frame()%>%
                        st_as_sf(coords=c("x","y"),crs=rasproj)
    }else{bathypoints <- data.frame(x=coordinates(bathy)[!is.na(values(bathy)),][1],
                                    y=coordinates(bathy)[!is.na(values(bathy)),][2])%>%
                         st_as_sf(coords=c("x","y"),crs=rasproj)}
                  
  nearest_point <- bathypoints[suppressMessages(st_nn(x,bathypoints,progress = FALSE))[[1]],]
  
  return(st_coordinates(nearest_point)%>%data.frame()%>%rename(long=X,lat=Y))
  
  } else {return(st_coordinates(x)%>%data.frame()%>%rename(long=X,lat=Y))}
      
}
