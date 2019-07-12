#####################################
## set working directory
#####################################

setwd("~/CARD-FISH/CARD-FISH_water_project/")

#set working directory in windows
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
###################################
## Load required libraries
###################################
library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(splines)
library(gridExtra)

###################################
## import environmental data and scale to center
###################################
#import data
metadata <- read.csv("PS99_samples_meta_EF_MC_nutrient_corrected.csv", sep = ",", dec = ".", header = TRUE)

#remove station SV2 from the dataset
metadata <- subset(metadata, !StationName == "SV2")

#list environmental parameters that need to be scaled
env.par= c("Temperature","Salinity", 
           "Chla_fluor", "NO3.NO2", 
           "NO2", "NO3", 
           "SiO3", "PO4",
           "NH4")

#scale parameters
scale1 <- function(x) scale(x, center = FALSE, scale = TRUE)[,1]

metadata %>% 
  mutate_at(env.par, scale1) -> env

# #Data standarization
# env_raw <- metadata
# 
# centre_scale2 <- function(x) {
#   scale(x, scale = FALSE)
# }
# env.par= c("Temperature","Salinity", "Chla_fluor", "NO3.NO2", "NO2", "NO3", "SiO3", "PO4","NH4")
# env=env_raw
# env[,env.par]<- apply(env[,env.par], MARGIN=2, FUN=centre_scale2)
# env <- subset(env, !StationName == "SV2")

###################################
## Import counts and calculate cell concentrations
###################################
#raw counts
raw.counts.SH <- read.csv("FOV_all_groups_SH.csv", sep = ",", dec = ".", header = TRUE)

#split sample name
raw.counts.SH <- raw.counts.SH %>% 
  separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-")

#remove SV2 station and ABYSS depth
raw.counts.SH <- subset(raw.counts.SH, !StationName == "SV2")
raw.counts.SH <- subset(raw.counts.SH, !Depth == "ABYSS")

###calculate cell concentration for each field of view
#microscope calculation factor
calc.factor <- 99515.5458411807
#DAPI
raw.counts.SH$DAPI_conc <- (raw.counts.SH$DAPI_Nr_Set*calc.factor)/raw.counts.SH$Volume
#FISH
raw.counts.SH$FISH_conc <- (raw.counts.SH$DAPI_Nr_SubSet*calc.factor)/raw.counts.SH$Volume

#calculate DAPI concetration per sample
raw.counts.SH %>% 
  group_by(StationName, Depth, Domain) %>% 
  summarise (DAPI.conc.mn = mean(DAPI_conc),
             DAPI.conc.md =  median(DAPI_conc), 
             DAPI.conc.sd = sd(DAPI_conc),
             FISH.conc.mn = mean(FISH_conc),
             FISH.conc.md =  median(FISH_conc), 
             FISH.conc.sd = sd(FISH_conc),
             n = length(DAPI_conc)) -> counts_all

# 
# counts_DAPI <- as.data.frame(as.list((aggregate(DAPI_Nr_Set~SAMPLE_NAME, 
#                                                 data=raw.counts.SH, 
#                                                 FUN = function(x) c(mn = mean(x), md =  median(x), sd = sd(x), n = length(x))))))
# 
# counts_DAPI.vol <- as.data.frame(as.list((aggregate(Volume~SAMPLE_NAME, 
#                                                     data=raw.counts.SH, 
#                                                     FUN = mean))))
# 
# counts_DAPI <- cbind(counts_DAPI,counts_DAPI.vol[,2])
# colnames(counts_DAPI)[6] <- "volume"
# 
# counts_DAPI$conc.mn <- (counts_DAPI$DAPI_Nr_Set.mn*calc.factor)/counts_DAPI$volume
# counts_DAPI$conc.md <- (counts_DAPI$DAPI_Nr_Set.md*calc.factor)/counts_DAPI$volume
# counts_DAPI$conc.sd <- (counts_DAPI$DAPI_Nr_Set.sd*calc.factor)/counts_DAPI$volume
# 
# counts_DAPI <- counts_DAPI %>% 
#   separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-")
# 
# counts_DAPI <- subset(counts_DAPI, !StationName == "SV2")


###FISH



# counts_FISH <- as.data.frame(as.list(aggregate(DAPI_Nr_SubSet ~SAMPLE_NAME, data=raw.counts.SH,
#                                                           FUN = function(x) c(mn = mean(x), md =  median(x), sd = sd(x), n = length(x)))))
# counts_FISH.vol <- as.data.frame(as.list((aggregate(Volume ~SAMPLE_NAME, data=raw.counts.SH,
#                                                                   FUN = mean))))
# 
# counts_FISH <- cbind(counts_FISH,counts_FISH.vol[,2])
# colnames(counts_FISH)[6] <- "volume"
# 
# counts_FISH$conc.mn <- (counts_FISH$DAPI_Nr_SubSet.mn*calc.factor)/counts_FISH$volume
# counts_FISH$conc.md <- (counts_FISH$DAPI_Nr_SubSet.md*calc.factor)/counts_FISH$volume
# counts_FISH$conc.sd <- (counts_FISH$DAPI_Nr_SubSet.sd*calc.factor)/counts_FISH$volume
# 
# counts_FISH <- counts_FISH %>%
#   separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-")
# counts_FISH <- subset(counts_FISH, !StationName == "SV2")

##################################
## call water layers and add regions and water masses
##################################
#Regions as in ice covered:EGC ice-free:WSC (Fadeev et al., 2019)
EGC <- c("EG1","EG4")
N <- c("N3","N4","N5")
WSC <- c("HG9","HG7","HG5","HG4","HG2","HG1")

#water masses
PSWw <- c("EG1","EG4", "HG4", "HG5", "HG7", "HG9", "N3", "N4", "N5")
AW <- c("HG2","HG1", "SV2")

counts_all$Region <- "EGC"
counts_all$Region[counts_all$StationName %in% N] <- "N"
counts_all$Region[counts_all$StationName %in% WSC] <- "WSC"


#add the longitude of each station from the metadata
for (station in metadata$StationName){
  counts_all$long <- NA
  counts_all$long[counts_all$StationName == station] <- metadata[metadata$StationName == station,]$Longitude..degrees_east.
    }

# counts_FISH$long[counts_FISH$StationName == "EG1"] <- -5.418
# counts_FISH$long[counts_FISH$StationName == "EG4"] <- -2.729
# counts_FISH$long[counts_FISH$StationName == "HG1"] <- 6.088
# counts_FISH$long[counts_FISH$StationName == "HG2"] <- 4.907
# counts_FISH$long[counts_FISH$StationName == "HG4"] <- 4.185
# counts_FISH$long[counts_FISH$StationName == "HG5"] <- 3.8
# counts_FISH$long[counts_FISH$StationName == "HG7"] <- 3.5
# counts_FISH$long[counts_FISH$StationName == "HG9"] <- 2.841
# counts_FISH$long[counts_FISH$StationName == "N4"] <- 4.508
# counts_FISH$long[counts_FISH$StationName == "N3"] <- 5.166
# counts_FISH$long[counts_FISH$StationName == "N5"] <- 3.062
# #counts_FISH$long[counts_FISH$StationName == "SV2"] <- 9.514
# 
# counts_DAPI$long[counts_DAPI$StationName == "EG1"] <- -5.418
# counts_DAPI$long[counts_DAPI$StationName == "EG4"] <- -2.729
# counts_DAPI$long[counts_DAPI$StationName == "HG1"] <- 6.088
# counts_DAPI$long[counts_DAPI$StationName == "HG2"] <- 4.907
# counts_DAPI$long[counts_DAPI$StationName == "HG4"] <- 4.185
# counts_DAPI$long[counts_DAPI$StationName == "HG5"] <- 3.8
# counts_DAPI$long[counts_DAPI$StationName == "HG7"] <- 3.5
# counts_DAPI$long[counts_DAPI$StationName == "HG9"] <- 2.841
# counts_DAPI$long[counts_DAPI$StationName == "N4"] <- 4.508
# counts_DAPI$long[counts_DAPI$StationName == "N3"] <- 5.166
# counts_DAPI$long[counts_DAPI$StationName == "N5"] <- 3.062
# #counts_DAPI$long[counts_DAPI$StationName == "SV2"] <- 9.514
# 

##################################
#save image for further analysis
##################################
save.image("card_fish_water_column_ps99.Rdata")

