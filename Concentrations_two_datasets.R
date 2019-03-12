###################################
##required libraries
###################################
library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)

#plot theme
theme_plot <- theme(axis.title.x = element_text(size = 18), 
                    axis.title.y = element_text(size = 18),
                    axis.text.x = element_text(size = 15),
                    axis.text.y = element_text(size = 15),
                    strip.text = element_text(size = 15),
                    legend.text= element_text(size = 12),
                    plot.title = element_text(size=18, lineheight=.8, face="bold", hjust = 0.5),
                    panel.background = element_rect(fill = 'white', colour = 'white'),
                    panel.spacing = unit(.05, "lines"),
                    panel.border = element_rect(color = "black", fill = NA, size = 1), 
                    strip.background = element_rect(color = "black", size = 1),
                    panel.grid.minor = element_line(linetype = 'dotted', color = 'grey', size= 0.5),
                    panel.grid.major = element_line(linetype = 'dotted', color = 'grey', size= 0.5))

###################################

###################################
## DAPI and FISH cell concentrations 
###################################
#import data set
raw.counts.SH <- read.csv("~/CARD-FISH/CARD-FISH_water_project/Concentration/single_hibr_counted/FOV_all_groups_SH.csv", sep = ",", dec = ".", header = TRUE)

raw.counts.DH <- read.csv("~/CARD-FISH/CARD-FISH_water_project/Concentration/double_hibr_counted/FOV_all_groups_DH.csv", sep = ",", dec = ".", header = TRUE)

#calculation factor
calc.factor <- 99515.5458411807

## 1.) Concentrations:

# a) DAPI and FISH cell concentrations in cells/mL in SH:
library(magrittr)
Concentration_data_SH <- raw.counts.SH %>%
  separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-")

Concentration_data_SH$conc.DAPI <- (Concentration_data_SH$DAPI_Nr_Set*calc.factor)/Concentration_data_SH$Volume
Concentration_data_SH$conc.FISH <- (Concentration_data_SH$DAPI_Nr_SubSet*calc.factor)/Concentration_data_SH$Volume

#write csv Concentration_data for single counted hibridized:

#write.csv(Concentration_data_SH, file = "~/CARD-FISH/CARD-FISH_water_project/Concentration/single_hibr_counted/concentrations_SH.csv") 

# b) DAPI and FISH cell concentrations in cells/mL in DH:
library()
Concentration_data_DH <- magrittr::raw.counts.SH %>%
  separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-")

Concentration_data_DH$conc.DAPI <- (Concentration_data_DH$DAPI_Nr_Set*calc.factor)/Concentration_data_DH$Volume
Concentration_data_DH$conc.FISH <- (Concentration_data_DH$DAPI_Nr_SubSet*calc.factor)/Concentration_data_DH$Volume

#write csv Concentration_data for single counted hibridized:

#write.csv(Concentration_data_DH, file = "~/CARD-FISH/CARD-FISH_water_project/Concentration/double_hibr_counted/concentrations_DH.csv") 

##boxplot of count ditribution
library(magrittr)
boxplot_data_SH <- raw.counts.SH %>%
  separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-")

boxplot_data_SH_EUB_DAPI <- dplyr::select(filter(boxplot_data_SH, Domain == "EUB"),c(StationName,Depth,Domain,SMP_No,FOV_No,DAPI_Nr_Set,DAPI_Nr_SubSet,Volume))
boxplot_data_SH_EUB_DAPI_DCM <- dplyr::select(filter(boxplot_data_SH_EUB_DAPI, Depth == "DCM"),c(StationName,Depth,Domain,SMP_No,FOV_No,DAPI_Nr_Set,DAPI_Nr_SubSet,Volume))
boxplot_data_SH_ARCH_DAPI <- dplyr::select(filter(boxplot_data_SH, Domain == "ARCH"),c(StationName,Depth,Domain,SMP_No,FOV_No,DAPI_Nr_Set,DAPI_Nr_SubSet,Volume))
boxplot_data_SH_ARCH_DAPI_DCM <- dplyr::select(filter(boxplot_data_SH_ARCH_DAPI, Depth == "DCM"),c(StationName,Depth,Domain,SMP_No,FOV_No,DAPI_Nr_Set,DAPI_Nr_SubSet,Volume))


boxplot_data_SH_EUB_DAPI_DCM$conc.DAPI <- (boxplot_data_SH_EUB_DAPI_DCM$DAPI_Nr_Set*calc.factor)/boxplot_data_SH_EUB_DAPI_DCM$Volume
boxplot_data_SH_EUB_DAPI_DCM$conc.FISH <- (boxplot_data_SH_EUB_DAPI_DCM$DAPI_Nr_SubSet*calc.factor)/boxplot_data_SH_EUB_DAPI_DCM$Volume

boxplot_data_SH_ARCH_DAPI_DCM$conc.DAPI <- (boxplot_data_SH_ARCH_DAPI_DCM$DAPI_Nr_Set*calc.factor)/boxplot_data_SH_ARCH_DAPI_DCM$Volume
boxplot_data_SH_ARCH_DAPI_DCM$conc.FISH <- (boxplot_data_SH_ARCH_DAPI_DCM$DAPI_Nr_SubSet*calc.factor)/boxplot_data_SH_ARCH_DAPI_DCM$Volume



EGC <- c("EG1","EG4")
HG <- c("HG9","HG7", "HG5", "HG4", "HG2","HG1")
N <- c("N3", "N4", "N5")
SV <- c("SV2")

EGCN <- c("EG1","EG4","N3", "N4", "N5")

boxplot_data_SH_EUB_DAPI_DCM$Region[boxplot_data_SH_EUB_DAPI_DCM$StationName %in% EGC] <- "EGC"
boxplot_data_SH_EUB_DAPI_DCM$Region[boxplot_data_SH_EUB_DAPI_DCM$StationName %in% HG] <- "HG"
boxplot_data_SH_EUB_DAPI_DCM$Region[boxplot_data_SH_EUB_DAPI_DCM$StationName %in% N] <- "N"
boxplot_data_SH_EUB_DAPI_DCM$Region[boxplot_data_SH_EUB_DAPI_DCM$StationName %in% SV] <- "SV"
#boxplot_data_SH_EUB_DAPI_DCM$Region[boxplot_data_SH_EUB_DAPI_DCM$StationName %in% EGCN] <- "EGCN"
boxplot_data_SH_ARCH_DAPI_DCM$Region[boxplot_data_SH_ARCH_DAPI_DCM$StationName %in% EGC] <- "EGC"
boxplot_data_SH_ARCH_DAPI_DCM$Region[boxplot_data_SH_ARCH_DAPI_DCM$StationName %in% HG] <- "HG"
boxplot_data_SH_ARCH_DAPI_DCM$Region[boxplot_data_SH_ARCH_DAPI_DCM$StationName %in% N] <- "N"
boxplot_data_SH_ARCH_DAPI_DCM$Region[boxplot_data_SH_ARCH_DAPI_DCM$StationName %in% SV] <- "SV"


EUB_FISH_boxplot <- ggplot(boxplot_data_SH_EUB_DAPI_DCM, aes(Region, conc.FISH))+geom_boxplot()+#facet_grid(Region,.)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "DAPI [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

EUB_DAPI_boxplot <- ggplot(boxplot_data_SH_EUB_DAPI_DCM, aes(Region, conc.DAPI))+geom_boxplot()+#facet_grid(Depth~Domain)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "CARD-FISH [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

FISH_scatterplot <- ggplot(boxplot_data_SH_EUB_DAPI_DCM, aes(x=Region, y=conc.FISH)) + geom_point()+
  geom_text(label=boxplot_data_SH_EUB_DAPI_DCM$StationName)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "CARD-FISH [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))
  
EUB_DAPI_scatterplot <- ggplot(boxplot_data_SH_EUB_DAPI_DCM, aes(x=Region, y=conc.DAPI)) + geom_point()+
  geom_text(label=boxplot_data_SH_EUB_DAPI_DCM$StationName)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "CARD-FISH [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

#for archaea:

ARCH_FISH_boxplot <- ggplot(boxplot_data_SH_ARCH_DAPI_DCM, aes(Region, conc.FISH))+geom_boxplot()+#facet_grid(Region,.)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "DAPI [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

ARCH_DAPI_boxplot <- ggplot(boxplot_data_SH_ARCH_DAPI_DCM, aes(Region, conc.DAPI))+geom_boxplot()+#facet_grid(Depth~Domain)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "CARD-FISH [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

ARCH_FISH_scatterplot <- ggplot(boxplot_data_SH_ARCH_DAPI_DCM, aes(x=Region, y=conc.FISH)) + geom_point()+
  geom_text(label=boxplot_data_SH_ARCH_DAPI_DCM$StationName)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "CARD-FISH [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

ARCH_DAPI_scatterplot <- ggplot(boxplot_data_SH_ARCH_DAPI_DCM, aes(x=Region, y=conc.DAPI)) + geom_point()+
  geom_text(label=boxplot_data_SH_ARCH_DAPI_DCM$StationName)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "CARD-FISH [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

#boxplots with average and error bars


#### 1.) DAPI comparison median and mean

counts.aggregated_DAPI <- as.data.frame(as.list((aggregate(DAPI_Nr_Set ~SAMPLE_NAME, data=raw.counts.SH, 
                                                           FUN = function(x) c(mn = mean(x), md =  median(x), sd = sd(x), n = length(x))))))
counts.aggregated_DAPI.volume <- as.data.frame(as.list((aggregate(Volume ~SAMPLE_NAME, data=raw.counts.SH, 
                                                                  FUN = mean))))

counts.aggregated_DAPI <- cbind(counts.aggregated_DAPI,counts.aggregated_DAPI.volume[,2])
colnames(counts.aggregated_DAPI)[6] <- "volume"

counts.aggregated_DAPI$conc.mn <- (counts.aggregated_DAPI$DAPI_Nr_Set.mn*calc.factor)/counts.aggregated_DAPI$volume
counts.aggregated_DAPI$conc.md <- (counts.aggregated_DAPI$DAPI_Nr_Set.md*calc.factor)/counts.aggregated_DAPI$volume
counts.aggregated_DAPI$conc.sd <- (counts.aggregated_DAPI$DAPI_Nr_Set.sd*calc.factor)/counts.aggregated_DAPI$volume

counts.aggregated_FISH <- as.data.frame(as.list(aggregate(DAPI_Nr_SubSet ~SAMPLE_NAME, data=raw.counts.SH, 
                                                          FUN = function(x) c(mn = mean(x), md =  median(x), sd = sd(x), n = length(x)))))
counts.aggregated_FISH.volume <- as.data.frame(as.list((aggregate(Volume ~SAMPLE_NAME, data=raw.counts.SH, 
                                                                  FUN = mean))))

counts.aggregated_FISH <- cbind(counts.aggregated_FISH,counts.aggregated_FISH.volume[,2])
colnames(counts.aggregated_FISH)[6] <- "volume"

counts.aggregated_FISH$conc.mn <- (counts.aggregated_FISH$DAPI_Nr_SubSet.mn*calc.factor)/counts.aggregated_FISH$volume
counts.aggregated_FISH$conc.md <- (counts.aggregated_FISH$DAPI_Nr_SubSet.md*calc.factor)/counts.aggregated_FISH$volume
counts.aggregated_FISH$conc.sd <- (counts.aggregated_FISH$DAPI_Nr_SubSet.sd*calc.factor)/counts.aggregated_FISH$volume

counts.aggregated_DAPI <- counts.aggregated_DAPI %>% 
  separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-")

counts.aggregated_FISH <- counts.aggregated_FISH %>% 
  separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-")

counts.aggregated_DAPI$Region[counts.aggregated_DAPI$StationName %in% EGC] <- "EGC"
counts.aggregated_DAPI$Region[counts.aggregated_DAPI$StationName %in% HG] <- "HG"
counts.aggregated_DAPI$Region[counts.aggregated_DAPI$StationName %in% N] <- "N"
counts.aggregated_DAPI$Region[counts.aggregated_DAPI$StationName %in% SV] <- "SV"

counts.aggregated_FISH$Region[counts.aggregated_FISH$StationName %in% EGC] <- "EGC"
counts.aggregated_FISH$Region[counts.aggregated_FISH$StationName %in% HG] <- "HG"
counts.aggregated_FISH$Region[counts.aggregated_FISH$StationName %in% N] <- "N"
counts.aggregated_FISH$Region[counts.aggregated_FISH$StationName %in% SV] <- "SV"

counts.aggregated_DAPI_DCM <- dplyr::select(filter(counts.aggregated_DAPI, Domain == "EUB"),c(StationName,Depth,Domain,DAPI_Nr_Set.mn,DAPI_Nr_Set.md,DAPI_Nr_Set.sd,DAPI_Nr_Set.n,volume,conc.mn,conc.md,conc.sd,Region))
counts.aggregated_DAPI_DCM <- dplyr::select(filter(counts.aggregated_DAPI_DCM, Depth == "DCM"),c(StationName,Depth,Domain,DAPI_Nr_Set.mn,DAPI_Nr_Set.md,DAPI_Nr_Set.sd,DAPI_Nr_Set.n,volume,conc.mn,conc.md,conc.sd,Region))

counts.aggregated_EUB_DCM <- dplyr::select(filter(counts.aggregated_FISH, Domain == "EUB"),c(StationName,Depth,Domain,DAPI_Nr_SubSet.mn,DAPI_Nr_SubSet.md,DAPI_Nr_SubSet.sd,DAPI_Nr_SubSet.n,volume,conc.mn,conc.md,conc.sd,Region))
counts.aggregated_EUB_DCM <- dplyr::select(filter(counts.aggregated_EUB_DCM, Depth == "DCM"),c(StationName,Depth,Domain,DAPI_Nr_SubSet.mn,DAPI_Nr_SubSet.md,DAPI_Nr_SubSet.sd,DAPI_Nr_SubSet.n,volume,conc.mn,conc.md,conc.sd,Region))

counts.aggregated_DAPI_ARCH_DCM <- dplyr::select(filter(counts.aggregated_DAPI, Domain == "ARCH"),c(StationName,Depth,Domain,DAPI_Nr_Set.mn,DAPI_Nr_Set.md,DAPI_Nr_Set.sd,DAPI_Nr_Set.n,volume,conc.mn,conc.md,conc.sd,Region))
counts.aggregated_DAPI_ARCH_DCM <- dplyr::select(filter(counts.aggregated_DAPI_ARCH_DCM, Depth == "DCM"),c(StationName,Depth,Domain,DAPI_Nr_Set.mn,DAPI_Nr_Set.md,DAPI_Nr_Set.sd,DAPI_Nr_Set.n,volume,conc.mn,conc.md,conc.sd,Region))

counts.aggregated_ARCH_DCM <- dplyr::select(filter(counts.aggregated_FISH, Domain == "ARCH"),c(StationName,Depth,Domain,DAPI_Nr_SubSet.mn,DAPI_Nr_SubSet.md,DAPI_Nr_SubSet.sd,DAPI_Nr_SubSet.n,volume,conc.mn,conc.md,conc.sd,Region))
counts.aggregated_ARCH_DCM <- dplyr::select(filter(counts.aggregated_ARCH_DCM, Depth == "DCM"),c(StationName,Depth,Domain,DAPI_Nr_SubSet.mn,DAPI_Nr_SubSet.md,DAPI_Nr_SubSet.sd,DAPI_Nr_SubSet.n,volume,conc.mn,conc.md,conc.sd,Region))


#############

#DAPI #add the region, plot the region label by stations 
counts.aggregated_DAPI_DCM_aggr  <- aggregate(conc.mn~StationName+Region,counts.aggregated_DAPI_DCM,mean)
counts.aggregated_DAPI_DCM_aggr1  <- aggregate(conc.sd~StationName+Region,counts.aggregated_DAPI_DCM,mean)
counts.aggregated_DAPI_DCM_aggr2 <- merge(counts.aggregated_DAPI_DCM_aggr, counts.aggregated_DAPI_DCM_aggr1)

EUB_DAPI_scatterplot_withsd2 <- ggplot(counts.aggregated_DAPI_DCM_aggr2, aes(x=Region, y=conc.mn)) + geom_point()+
  geom_text(label=counts.aggregated_DAPI_DCM_aggr2$StationName)+
  geom_errorbar(aes(ymin = conc.mn-conc.sd, ymax = conc.mn+conc.sd), width = 0.15)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "CARD-FISH [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

#EUB
counts.aggregated_EUB_DCM_aggr  <- aggregate(conc.mn~StationName+Region,counts.aggregated_EUB_DCM,mean)
counts.aggregated_EUB_DCM_aggr1  <- aggregate(conc.sd~StationName+Region,counts.aggregated_EUB_DCM,mean)
counts.aggregated_EUB_DCM_aggr2 <- merge(counts.aggregated_EUB_DCM_aggr1, counts.aggregated_EUB_DCM_aggr)

EUB_scatterplot_withsd2 <- ggplot(counts.aggregated_EUB_DCM_aggr2, aes(x=Region, y=conc.mn)) + geom_point()+
  geom_text(label=counts.aggregated_EUB_DCM_aggr2$StationName)+
  geom_errorbar(aes(ymin = conc.mn-conc.sd, ymax = conc.mn+conc.sd), width = 0.15)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "CARD-FISH [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

#ARCH
counts.aggregated_ARCH_DCM_aggr  <- aggregate(conc.mn~StationName+Region,counts.aggregated_ARCH_DCM,mean)
counts.aggregated_ARCH_DCM_aggr1  <- aggregate(conc.sd~StationName+Region,counts.aggregated_ARCH_DCM,mean)
counts.aggregated_ARCH_DCM_aggr2 <- merge(counts.aggregated_ARCH_DCM_aggr1, counts.aggregated_ARCH_DCM_aggr)

ARCH_scatterplot_withsd2 <- ggplot(counts.aggregated_ARCH_DCM_aggr2, aes(x=Region, y=conc.mn)) + geom_point()+
  geom_text(label=counts.aggregated_ARCH_DCM_aggr2$StationName)+
  geom_errorbar(aes(ymin = conc.mn-conc.sd, ymax = conc.mn+conc.sd), width = 0.15)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "CARD-FISH [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

#########
#Patterns of specific groups vertical and horizontally, until 2500 faceting according to region and add error bars, 

#EUB

counts.EUB <- dplyr::select(filter(counts.aggregated_FISH, Domain == "EUB"),c(StationName,Depth,Domain,DAPI_Nr_SubSet.mn,DAPI_Nr_SubSet.md,DAPI_Nr_SubSet.sd,DAPI_Nr_SubSet.n,volume,conc.mn,conc.md,conc.sd,Region))
counts.EUB_aggr  <- aggregate(conc.mn~Region+Depth,counts.EUB,mean)
counts.EUB_aggr1  <- aggregate(conc.mn~Region+Depth,counts.EUB,mean)
counts.EUB_aggr2 <- merge(counts.EUB_aggr1, counts.EUB_aggr)
#counts.EUB_aggr3 <- aggregate(conc.mn~Region+Depth,counts.EUB,mean)
counts.EUB_aggr2$Depth<- factor(counts.EUB_aggr2$Depth, 
                             levels = rev(c("DCM", "EPI", "MESO", "BATHY", "ABYSS")))
counts.EUB_aggr2$Depth<- as.numeric(counts.EUB_aggr2$Depth, 
                                levels = rev(c("DCM", "EPI", "MESO", "BATHY", "ABYSS")))


abund_abs_deep_groups.p <- ggplot()+
  geom_point(data = counts.EUB_aggr2, aes(y = conc.mn, x = Depth, colour = Region))+
  geom_line(data = counts.EUB_aggr2, aes(y = conc.mn, x = Depth, colour = Region),size = 0.8)+
  # geom_line(data = data.deep.agg.melt[data.deep.agg.melt$variable %in% includep,], linetype=1, aes(y = value, x = Depth, group = variable, colour = variable),size = 0.8)+
  # geom_line(data = data.deep.agg.melt[data.deep.agg.melt$variable %in% includep,], linetype=1, aes(y = value, x = Depth, group = variable, colour = variable),size = 0.8)+
  #scale_x_reverse(breaks = c(50,100,1000,2500,4000,5500))+
  #scale_y_continuous(breaks = c(0.1,0.25,0.5,0.75,1))+
  coord_flip()
  #scale_color_viridis(option="plasma", discrete = TRUE)+ #uncommented to get default color for more groups
  #facet_grid(~Region)
#theme_classic(base_size = 12)
