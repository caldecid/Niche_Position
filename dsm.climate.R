# Function for cleaning and formating dsm and climate ---------------------



##package
library(ape)
library(sf)
library(tidyr)
library(dplyr)
library(raster)
library(maptools)
library(rgdal)
library(stringr)
library(rgeos)
library(cleangeo)
library(CENFA)
library(usdm)
library(svMisc)



# dsm and climate function ------------------------------------------------

dsm.climate.vectorized <- function(realm, dsm, w){
  #realm = realm from wwf terrestrial ecosystems (see metadata)
  #dsm = SpatialPolygondataframe of the distribution species model 
  #w = weather raster brick
  
  ########################Dealing with weather###############################
  
  ##empty list for storing the climatic variables for each biome
  
  clim.l <- list()
  
  ##empty list for storing the highly correlated climatic variables
  
  clim.vif <- list()
  
  ##loop for cropping and storing the weather of each biome
  
  for(i in 1:length(realm)){
    
    clim.l[[i]] <- crop(x = w, y = extent(realm[i, ]))
    
    ##mask for cropping the exact biome extension
    
    clim.l[[i]] <- mask(clim.l[[i]], realm[i, ])
    
    ##detecting highly correlated variables
    
    clim.vif[[i]] <- vifstep(clim.l[[i]])
    
    ##excluding highly correlated variables
    clim.l[[i]] <- exclude(clim.l[[i]], clim.vif[[i]])
  }
  
  
  ###############Dealing with dsm###############################
  
  #eliminate island species
  
  not.island <- which(is.na(dsm@data[["ISLAND"]])) ##vector with nonisland sp
  
  dsm.continental <- dsm[not.island, ]
  
  ###exclude the species which range <0.11 due to not running in ENFA
  
  dsm.min.area <- dsm.continental[gArea(dsm.continental, byid = TRUE)
                                  > 0.11, ]
  
  ###coercing coordinate reference system 
  crs(dsm.min.area) <- crs(clim.l[[1]])
  
  ##we are going to consider just extant, native, and resident species
  ##Presence
  dsm.min.area <- dsm.min.area[dsm.min.area@data$PRESENCE == 1, ]
  ##origin
  dsm.min.area <- dsm.min.area[dsm.min.area@data$ORIGIN == 1, ]
  ##seasonal
  dsm.min.area <- dsm.min.area[dsm.min.area@data$SEASONAL == 1, ] 
  
  ###obtaining dsm data
  dsm.data <- dsm.min.area@data
  
  ###simplyfying the dsm
  dsm.simple <- gSimplify(dsm.min.area, tol = 0.00001)
  
  ###gSimplify transforms the dsm into a spatialpolygon, so we need the data
  dsm.simple.data <- SpatialPolygonsDataFrame(dsm.simple, data = dsm.data)
  
  ###aggregating for dissolving the polygons of repeated species
  dsm.ag <- aggregate(dsm.simple.data, by = "BINOMIAL")
  
  ##saving the output
  save(list(clim.l, dsm.ag), 
       file = paste0(as.character(dsm), "climate", ".RData"))
  
  ##returning the climate and dsm
  return(list(clim.l, dsm.ag))
}

