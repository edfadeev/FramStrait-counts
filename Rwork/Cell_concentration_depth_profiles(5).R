#####################################
## set working directory
#####################################

setwd("~/CARD-FISH/CARD-FISH_water_project/")

load("card_fish_water_column_ps99.Rdata")
###################################
## Load required libraries
###################################
library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)

##################################
# Depth profiles cell concentration
##################################
deep.data <- counts_FISH
counts.deep <- deep.data[,c("StationName", "Depth", "Domain", "conc.mn", "Region")]

counts.deep.agg <- aggregate(conc.mn~Domain+Depth+Region, counts.deep,mean)
counts.deep.agg$conc.mn <- log10(counts.deep.agg$conc.mn)


counts.deep.agg$Depth<- factor(counts.deep.agg$Depth, 
                           levels = rev(c("DCM", "EPI", "MESO", "BATHY", "ABYSS")))

counts.deep.agg_egc <- subset(counts.deep.agg[counts.deep.agg$Region %in% c("EGC"),])
counts.deep.agg_wsc <- subset(counts.deep.agg[counts.deep.agg$Region %in% c("WSC"),])

# fit text in facets:
swr = function(string, nwrap=20) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)
domain.labs = swr(domain.labs)

domain.labs <- c("Alteromonadaceae Colwelliaceae Pseudoalteromonadaceae", "Archaea", "Bacteroidetes", "Chloroflexi", "Thaumarcheota", "Deltaproteobacteria", "Bacteria", "Gammaproteobacteria", "Opitutae", "Polaribacter", "Rhodobacteraceae", "SAR11", "SAR202", "SAR324", "Marinimicrobia (SAR406)", "Verrucomicrobia")
names(domain.labs) <- c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")
#test <-  aggregate(dapi.mn~Domain+Depth+Region, counts.rel.ab,mean)
#test_2 <- aggregate(fish.mn~Domain+Depth+Region, counts.rel.ab,mean) 
#test_3 <- test_2 %>% left_join(test, by=c("Domain", "Depth", "Region")) 
Plot_depth_profiles_absolute_abundance.p <- ggplot(data = counts.deep.agg, aes(y = conc.mn, x = Depth, color =Region, group=Region))+
  geom_line(aes(),linetype=1)+
  coord_flip()+
  facet_wrap(.~Domain, labeller = labeller(Domain = domain.labs))+
  labs(title = , y = "Cell density [log10(cells/ml)]", x = "", color = "Region")+ 
  scale_x_discrete(labels=c("DCM" = "Surface", "EPI" = "Epipelagic","MESO" = "Mesopelagic", "BATHY" = "Bathypelagic", "ABYSS" = "Abyssopelagic"))+
  scale_color_manual(breaks = c("EGC", "WSC"),labels = c("Ice-covered", "Ice-free"), values = c("blue", "red"))+
  theme(legend.position = "bottom")+
  theme_classic(base_size = 12)
Plot_depth_profiles_absolute_abundance.p 

save.image("card_fish_water_column_ps99.Rdata")

###rel. abundance

counts.deep.rel <- fish.rel %>% left_join(dapi.rel ,by =c("StationName", "Domain", "Depth", "Region"))
counts.deep.rel.f <- counts.deep.rel[, c("StationName", "Depth", "Domain", "fish.mn", "Region")]
counts.deep.rel.d <- counts.deep.rel[, c("StationName", "Depth", "Domain", "dapi.mn", "Region")]

counts.deep.rel.f %>% spread(Domain, fish.mn) -> counts.deep.rel.f.sp
counts.deep.rel.d %>% spread(Domain, dapi.mn) -> counts.deep.rel.d.sp

counts.deep.rel.f.agg <- aggregate(fish.mn~Domain+Depth+Region, counts.deep.rel.f,mean)
counts.deep.rel.d.agg <- aggregate(dapi.mn~Domain+Depth+Region, counts.deep.rel.d,mean)

counts.deep.f.d.agg <- counts.deep.rel.f.agg %>% left_join(counts.deep.rel.d.agg ,by =c("Domain", "Depth", "Region"))
counts.deep.f.d.agg$proportion <- counts.deep.f.d.agg$fish.mn*100/counts.deep.f.d.agg$dapi.mn

counts.deep.f.d.agg$Depth<- factor(counts.deep.f.d.agg$Depth, 
                               levels = rev(c("DCM", "EPI", "MESO", "BATHY", "ABYSS")))

Plot_depth_profiles_re_abund.p <- ggplot(data = counts.deep.f.d.agg, aes(y = proportion, x = Depth, color =Region, group=Region))+
  geom_line(aes(),linetype=1)+
  coord_flip()+
  facet_wrap(.~Domain, labeller = labeller(Domain = domain.labs))+
  labs(title = , y = "Relative abundance (% DAPI)", x = "", color = "Region")+ 
  scale_x_discrete(labels=c("DCM" = "Surface", "EPI" = "Epipelagic","MESO" = "Mesopelagic", "BATHY" = "Bathypelagic", "ABYSS" = "Abyssopelagic"))+
  scale_color_manual(breaks = c("EGC", "WSC"),labels = c("Ice-covered", "Ice-free"), values = c("blue", "red"))+
  theme(legend.position = "bottom")+
  theme_classic()
Plot_depth_profiles_re_abund.p 



