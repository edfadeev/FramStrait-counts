###################################
## Load required libraries
###################################
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")


###################################
## defined functions
###################################
#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

#microscope calculation factor
calc.factor <- 99515.5458411807

###################################
## Import counts and calculate cell concentrations
###################################
#raw counts
raw.counts <- read.csv("./Rwork/FOV_all_groups_SH.csv", sep = ",", dec = ".", header = TRUE)

#split sample name and remove SV2 station and ABYSS depth and calculate concentrations per FOV
raw.counts.SH <- raw.counts %>% 
  separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-") %>% 
  filter(StationName != "SV2" & Depth != "ABYSS") %>% 
  mutate(DAPI_conc= (DAPI_Nr_Set*calc.factor)/Volume,
         FISH_conc = (DAPI_Nr_SubSet*calc.factor)/Volume)


#calculate concetrations per sample and add regions
counts_all <-  raw.counts.SH %>% 
  filter(FISH_conc > 0) %>% 
  mutate(Region = ifelse(StationName %in% c("EG1","EG4"), "EGC", 
                         ifelse(StationName %in% c("N3","N4","N5"), "N", "WSC")),
         Depth = factor(Depth, levels = c("SRF","EPI","MESO","BATHY"))) %>% 
  group_by(Region, StationName, Depth, Domain) %>% 
  summarise(DAPI.conc.mn = mean(DAPI_conc),
            DAPI.conc.md =  median(DAPI_conc), 
            DAPI.conc.sd = sd(DAPI_conc),
            FISH.conc.mn = mean(FISH_conc),
            FISH.conc.md =  median(FISH_conc), 
            FISH.conc.sd = sd(FISH_conc),
            n = length(DAPI_conc))

##################################
# Depth profiles by taxa 
##################################
taxa.depth <- c("SAR11", "BACT", "GAM", "CREN",  "VER", "DELTA", "SAR324",   "SAR202", "SAR406")


#calculate mean conc., adjust the depth
depth_FISH_by_regions <- counts_all%>% ungroup() %>% 
  filter(Domain %in% taxa.depth) %>% 
  group_by(Region, Depth, Domain) %>% 
  summarise (mean.abund = mean(FISH.conc.mn),
             se.abund = se(FISH.conc.mn))%>% ungroup() %>% 
  mutate_if(is.factor, as.character) %>% 
  mutate(Depth=replace(Depth, Depth=="SRF", 20),
         Depth=replace(Depth, Depth=="EPI", 100),
         Depth=replace(Depth, Depth=="MESO", 1000),
         Depth=replace(Depth, Depth=="BATHY", 2000)) %>% 
  mutate(Depth = as.numeric(Depth),
         Domain = factor(Domain, levels = c("SAR11","BACT", "GAM", "CREN",  "VER", "OPI","DELTA", "SAR324", "CFX",  "SAR202", "SAR406")))

Depth_profiles_abs_abundance.p <- ggplot(depth_FISH_by_regions, aes(y = mean.abund, x = Depth, color = Region, group = Region))+
  geom_point()+
  geom_line()+
  #stat_smooth(method = "glm", formula = y ~ poly(x, 3), se = FALSE)+
  #geom_line(data = spline_int, aes(x = x, y = y))+
  #geom_smooth(aes(colour = Region), method = "auto")+
  scale_y_log10()+
  scale_x_reverse()+
  coord_flip()+
  facet_wrap(~Domain)+
  labs(y = "Cell concentration (cells/mL)]", x = "Depth (m)")+
  scale_color_manual(breaks = c("EGC", "N", "WSC"), values = c("blue", "gray", "red"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")
