require(dplyr)
require(tidyr)

##########################################
# Import and preprocess cell counts tables
##########################################
#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

#microscope calculation factor
calc.factor <- 99515.5458411807

## Import raw counts
raw.counts <- read.csv("./paper_scripts/Data/FOV_all_groups_SH.csv",
                       sep = ",", dec = ".", header = TRUE)

#Split sample name and remove SV2 station and ABYSS depth and calculate counts per FOV
raw.counts.SH <- raw.counts %>% 
  separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-") %>% 
  filter(StationName != "SV2" & Depth != "ABYSS") %>% 
  mutate(DAPI_conc= (DAPI_Nr_Set*calc.factor)/Volume,
         FISH_conc = (DAPI_Nr_SubSet*calc.factor)/Volume)

#Calculate abundances per sample and add regions
counts_all <-  raw.counts.SH %>% 
  filter(FISH_conc > 0) %>% 
  mutate(Region = ifelse(StationName %in% c("EG1","EG4"), "EGC", 
                         ifelse(StationName %in% c("N3","N4","N5"), "N", 
                                ifelse(StationName %in% c("N3","N4","N5"), "N", "WSC"))),
         Depth = factor(Depth, levels = c("SRF","EPI","MESO","BATHY"))) %>% 
  group_by(Region, StationName, Depth, Domain) %>% 
  summarise(DAPI.conc.mn = mean(DAPI_conc),
            DAPI.conc.md =  median(DAPI_conc), 
            DAPI.conc.sd = sd(DAPI_conc),
            FISH.conc.mn = mean(FISH_conc),
            FISH.conc.md =  median(FISH_conc), 
            FISH.conc.sd = sd(FISH_conc),
            n = length(DAPI_conc))%>%
            ungroup %>%
            as.data.frame()

