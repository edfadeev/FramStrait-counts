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

Plot_depth_profiles.p <- ggplot()+
  geom_point(data = counts.deep.agg, aes(y = conc.mn, x = Depth, colour = Domain))+
  geom_line(data = counts.deep.agg, linetype=1, aes(y = conc.mn, x = Depth, colour = Domain, group = Domain),size = 0.8)+ 
  coord_flip()+
  facet_grid(Region~Domain)+
  theme_plot

##################################
# Concentration surfaces of groups
##################################

dcm.group <- subset(counts_FISH[counts_FISH$Depth %in% c("DCM"),])
dcm.group.data <- aggregate(conc.mn~Region+Domain, dcm.group,mean)
dcm.group.data$Domain <- factor(dcm.group.data$Domain, levels = rev(c("EUB","SAR11","BACT","GAM","VER", "OPI", "POL", "ROS", "CFX", "ALT", "ARCH", "DELTA", "SAR324", "CREN", "SAR202", "SAR406" )))
dcm.group.data$conc.mn <- log10(dcm.group.data$conc.mn)

big_diff <- dcm.group.data %>% 
  spread(Region, conc.mn) %>% 
  group_by(Domain) %>% 
  mutate(Diff = WSC-EGC,
         Diffx = Diff/EGC) 

right_label <- dcm.group.data %>%
  group_by(Domain) %>%
  top_n(1)

plot_label <- big_diff %>%
  select(Domain, Diffx) %>%
  right_join(right_label)

surface_differences <- ggplot(dcm.group.data, aes(conc.mn, Domain)) +
  geom_line(aes(group = Domain)) +
  geom_point(aes(color = Region)) +
  scale_colour_manual(values=c("EGC" = "blue", "WSC" = "red")) +
  geom_text(data = plot_label, aes(label = paste0(scales::percent(round(Diffx,2)))),
            size = 5, hjust = -.5)+
  scale_x_continuous(limits = c(2, 6.5), breaks = seq(2, 6.5, by = 0.5)) +
  xlab("Log10(cells/ml)") +
  scale_y_discrete(labels=rev(c("Bacteria","SAR11","Bacteroidetes","Gammaproteobacteria","Verrucomicrobia", "Opitutae", "Polaribacter", "Rhodobacteraceae", "Chloroflexi", "Ateromonadales", "Archaea", "Deltaproteobacteria", "SAR324", "Thaumarchaeota", "SAR202", "SAR406"))) +
  labs(y = "") +
  theme_plot

                     
save.image("card_fish_water_column_ps99.Rdata")


