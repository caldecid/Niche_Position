#packages
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
library(stringr)
library(progress)


# ENFA within biomes function ---------------------------------------------


enfa.within.biomes <- function(realm, dsm, w){
  #realm = one realm from wwf terrestrial ecosystems divided by its biomes
  #dsm = SpatialPolygondataframe of the distribution species model 
  #w = weather raster brick
  
  ######first section of the function is to clean and prepare data ##########
  
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
  
  ###naming climates with the biome names
  names(clim.l) <- realm@data$BIOME_D
  
  ###############Dealing with dsm###############################
  
  #eliminate island species
  
  not.island <- which(is.na(dsm@data[["ISLAND"]])) ##vector with nonisland sp
  
  dsm.continental <- dsm[not.island, ]
  
  ###exclude the species which range <0.11 due to not running in ENFA
  
  dsm.min.area <- dsm.continental[gArea(dsm.continental, byid = TRUE)
                                  > 0.11, ]
  
  ###coercing coordinate reference system 
  crs(dsm.min.area) <- crs(clim.l)
  
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
  print("aggregating dsm")
  
  dsm.ag <- try(aggregate(dsm.simple.data, by = "BINOMIAL"))
  
  print("finishing of aggregating dsm")
  ###################data prepared for ENFA####################################
  
  ######PAM for determining the presence of species within biomes
  
  #empty matrix
  pam.biomes <- matrix(NA, nrow = length(dsm.ag), ncol = length(realm)+1)
  
  #naming columns and first row
  colnames(pam.biomes) <- c("species", realm@data$BIOME_D)
  pam.biomes[,1] <- dsm.ag@data$BINOMIAL
  
  ###matching projections
  crs(dsm.ag) <- crs(realm)
  ###Loop for building the PAM of the biomes
  for(i in 1:length(realm)) {
    pam.biomes[,i+1] <- gIntersects(realm[i,], dsm.ag, byid = TRUE)
  }
  
  ###list for keeping the intersection polygons of each species in each biome
  sp.polygon.list <- vector("list", length = length(realm))
  
  ##loop for creating a spatialpolygon for the species dsm occupyied within 
  ##each biome
  for(i in seq_along(sp.polygon.list)){
    sp.polygon.list[[i]] <- try(gIntersection(realm[i,],
                                dsm.ag[which(pam.biomes[,i+1] == TRUE),],
                                              byid = TRUE))
  }
  
  ##loop for obtaining the IDs and building a dataframe for the polygons
  data.biome.list <- vector("list", length = length(realm))
  for(i in seq_along(sp.polygon.list)){
   data.biome.list[[i]] <- try(data.frame(BINOMIAL = dsm.ag$BINOMIAL[which(pam.biomes[,i+1]
                         == TRUE)],row.names = sapply(slot(sp.polygon.list[[i]],
                                      "polygons"),function(x) slot(x, "ID"))))
  }
  
  ##loop for creating a spatialpolygondataframe
  sp.biome.list <- vector("list", length = length(realm))
  for(i in seq_along(sp.biome.list)){
    sp.biome.list[[i]] <- try(SpatialPolygonsDataFrame(sp.polygon.list[[i]],
                                             data = data.biome.list[[i]]))
  }
  
  ##naming the lists with the biome names
  names(sp.biome.list) <- realm@data$BIOME_D
  
  ##removing elements (biomes) without species
  sp.biome.list <- sp.biome.list[sapply(sp.biome.list,
                            function(x)class(x)== "SpatialPolygonsDataFrame")]
  
  ##In climate, removing biomes that are not represented by species
  clim.l <- clim.l[which(names(clim.l) %in% names(sp.biome.list))]
  
  ##empty list for storing the ENFA results for each species in each biome
  enfa.bio <- vector("list", length(sp.biome.list))##one main list for each biome
  
  ##naming enfa.bio
  names(enfa.bio) <- names(sp.biome.list)
  
  for(i in 1:length(sp.biome.list)){
    enfa.bio[[i]] <- vector("list", length(sp.biome.list[[i]])) 
    #a nested list for each species present in each biome
  }
  
  ###print progress
  print("starting ENFA")
  
  ##loop for performing the ENFA analysis
  
   pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)",
                         total = length(sp.biome.list))
  
  for(j in 1:length(sp.biome.list)){
      pb$tick()
      for(i in 1:length(sp.biome.list[[j]])){
        enfa.bio[[j]][[i]] <- suppressWarnings(try(enfa(x = clim.l[[j]], 
                                            s.dat = sp.biome.list[[j]][i,],
                                            parallel = TRUE, n = 3)))
      }
      
    }

  
  ########################ENFA results###################################
  
  # list storing the principal results of ENFA ------------------------------
  
  ##empty list for storing the results
  output.l <-  vector("list", length(enfa.bio))
  
  for(j in 1:length(enfa.bio)){
    print(j)
    output.l[[j]] <- data.frame(matrix(NA, nrow = length(enfa.bio[[j]]),
                                       ncol = 4))
    colnames(output.l[[j]]) <- c("biome", "species", "marginality", 
                                 "specialization")
    for(i in 1:length(enfa.bio[[j]])){
      output.l[[j]][i,1] <- names(sp.biome.list)[j] ##name biome
      output.l[[j]][i,2] <- sp.biome.list[[j]]@data$BINOMIAL[i] ##name sp
      if(class(enfa.bio[[j]][[i]]) != "enfa"){
        output.l[[j]][i,3] <- NA
        output.l[[j]][i,4] <- NA
      } else{
        output.l[[j]][i,3] <- enfa.bio[[j]][[i]]@marginality
        output.l[[j]][i,4] <- enfa.bio[[j]][[i]]@specialization
      }
      
    }
    
  }  
  
  results <- do.call("rbind", output.l)##results without weighted means
  
  ##dropping NAs
  results <- drop_na(results)
  
  print("ENFA ended")
  
  ############weighting marginality and specialization according to area######
  
  ##data frame with the species which are present in only one biome
  results.unique <- results[ave(seq_along(results$species), 
                                results$species, FUN = length) == 1, ]%>%
                            arrange(species)
  
  ###species which are  present in more than one biome nested by biome
  results.biome <- results[-which(results$species %in% 
                          results.unique$species), ] %>%
                      group_by(biome) %>%
                      nest()
  
  print("starting measuring Area")
  
  ##excluding biomes which have unique species
  sp.biome.list.1 <- sp.biome.list[which(names(sp.biome.list)
                                         %in% results.biome$biome)]
  ##calculating area of each species withing each biome
  for(i in seq_along(sp.biome.list.1)){
    
    results.biome[[2]][[i]]$area <- gArea(sp.biome.list.1[[i]]
                               [which(sp.biome.list.1[[i]]@data$BINOMIAL %in%
                               results.biome[[2]][[i]]$species),], byid = TRUE)
  }
  
  print("ending measuring Area")
  
  ###unnesting by biome and nesting by species
  results.biome <- results.biome %>% unnest(data) %>% 
                    group_by(species) %>% nest()
  
  
  print("starting weighting")
  
  ##weighting marginality and specialization by area
  
  ##creating columns for storing the weighting variables
  results.biome$marginality <- rep(NA, nrow(results.biome))
  
  results.biome$specialization <- rep(NA, nrow(results.biome))
  
  
#weighted marginality for each species
  
  for(j in 1:nrow(results.biome)){
    results.biome$marginality[j] <- 
      stats::weighted.mean(results.biome[[2]][[j]]$marginality,
             (results.biome[[2]][[j]]$area)/sum(results.biome[[2]][[j]]$area))
    
    #weighted specialization
    results.biome$specialization[j] <- 
      stats::weighted.mean(results.biome[[2]][[j]]$specialization,
           (results.biome[[2]][[j]]$area)/sum(results.biome[[2]][[j]]$area))
  }
 
  print("unifying dataframes")
  ##renaming columns of results.unique
  results.unique <- results.unique %>%
    rename(data = biome) %>%
    mutate(data = as.list(data))
  
  ##binding data sets
  results.total <- bind_rows(results.biome, results.unique)
  ##saving outputs
  save(enfa.bio, results.total, sp.biome.list, file = "enfa.results.RData")
  
  return(list(enfa.bio, results.total, sp.biome.list))
  
}

  