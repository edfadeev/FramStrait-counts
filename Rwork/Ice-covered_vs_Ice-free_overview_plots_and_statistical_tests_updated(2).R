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
library(splines)
library(gridExtra)
library(ggbiplot)

##################################
#plot theme
##################################
theme_plot <- theme(axis.title.x = element_text(size = 14), 
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 14),
                    strip.text = element_text(size = 14),
                    legend.text= element_text(size = 12),
                    plot.title = element_text(size=14, lineheight=.8, face="bold", hjust = 0.5),
                    panel.background = element_rect(fill = 'white', colour = 'white'),
                    panel.spacing = unit(.05, "lines"),
                    panel.border = element_rect(color = "black", fill = NA, size = 1), 
                    strip.background = element_rect(color = "black", size = 1),
                    panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5),
                    panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))


##################################
## Overview plots
##################################
#scatter of cell concentration all groups all depths
counts_FISH$Depth <- factor(counts_FISH$Depth, levels = c("DCM","EPI","MESO","BATHY","ABYSS"))
counts_DAPI$Depth <- factor(counts_DAPI$Depth, levels = c("DCM","EPI","MESO","BATHY","ABYSS"))

#FISH
counts_FISH_aggr  <- aggregate(conc.mn~StationName+Region+Depth+long,counts_FISH,mean)
counts_FISH_aggr1  <- aggregate(conc.sd~StationName+Region+Depth+long,counts_FISH,mean)
counts_FISH_aggr2 <- merge(counts_FISH_aggr1, counts_FISH_aggr)

Scatter_FISH_counts <- ggplot(counts_FISH_aggr2, aes(x=long, y=conc.mn)) + geom_point()+facet_grid(Depth~.)+
  geom_text(label=counts_FISH_aggr2$StationName, nudge_y= 1)+
  geom_errorbar(aes(ymin = conc.mn-conc.sd, ymax = conc.mn+conc.sd), width = 0.15)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "CARD-FISH [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

#DAPI
counts_DAPI_aggr  <- aggregate(conc.mn~StationName+Region+Depth+long,counts_DAPI,mean)
counts_DAPI_aggr1  <- aggregate(conc.sd~StationName+Region+Depth+long,counts_DAPI,mean)
counts_DAPI_aggr2 <- merge(counts_DAPI_aggr1, counts_DAPI_aggr)

Scatter_DAPI_counts <- ggplot(counts_DAPI_aggr2, aes(x=long, y=conc.mn)) + geom_point()+facet_grid(Depth~.)+
  geom_text(label=counts_DAPI_aggr2$StationName, nudge_y= 1)+
  geom_errorbar(aes(ymin = conc.mn-conc.sd, ymax = conc.mn+conc.sd), width = 0.15)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "CARD-FISH [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))
Scatter_DAPI_counts

#Boxplots of cell concentration DCM:
#FISH
Boxplot_FISH_counts <- ggplot(subset(counts_FISH, Depth %in% c("DCM")), aes(Region, conc.mn))+
  geom_boxplot()+facet_grid(Depth~.)+ geom_text(aes(label=Domain), size=3, position = position_jitter(w = 0.3))+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "cell conc. [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))
Boxplot_FISH_counts

Boxplot_FISH_counts_all <- ggplot(counts_FISH, aes(Region, conc.mn))+
  geom_boxplot()+facet_grid(Depth~.)+ geom_text(aes(label=Domain), size=3, position = position_jitter(w = 0.3))+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "cell conc. [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))
Boxplot_FISH_counts_all

#DAPI 
Boxplot_DAPI_counts <- ggplot(subset(counts_DAPI, Depth %in% c("DCM")), aes(Region, conc.mn))+
  geom_boxplot()+facet_grid(Depth~.)+ geom_text(aes(label=Domain), size=3, position = position_jitter(w = 0.3))+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "cell conc. [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

Boxplot_DAPI_counts_all <- ggplot(counts_DAPI, aes(Region, conc.mn))+
  geom_boxplot()+facet_grid(Depth~.)+ geom_text(aes(label=Domain), size=3, position = position_jitter(w = 0.3))+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "cell conc. [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

##################################
## East vs West 
##################################
#Plot of bacterial (EUB) cell concentration all groups Ice-covered vs Ice-free region 

counts_FISH_EUB <- subset(counts_FISH, Domain %in% c("EUB"))
counts_FISH_EUB  <- aggregate(conc.mn~StationName+Region+Depth+long,counts_FISH_EUB,mean)

 data.fish <-counts_FISH_aggr2 
#data.fish$conc.mn <- log10(data.fish$conc.mn)
#data.fish$smooth <- 1

counts_FISH_EUB$conc.mn <- log10(counts_FISH_EUB$conc.mn)
counts_FISH_EUB$smooth <- 1

Plot_eub_east_west.p <- ggplot(na.omit(subset(counts_FISH_EUB, Depth %in% c("DCM", "EPI", "MESO", "BATHY"))), aes(x = long, y = conc.mn, shape = Depth, group = Depth))+scale_shape_manual(values=c(19, 17, 15, 3), labels= c("Surface", "Epipelagic", "Mesopelagic", "Bathypelagic"))+  
  geom_point(aes(colour = Region), size =4)+
  xlab("Longitude [°East]")+
  ylab("Cell density [log10(cells/ml)]")+
  geom_text(aes(label = StationName), nudge_y= 0.08)+
  geom_smooth(aes(linetype = Depth), method = "gam", formula = y ~ ns(x,2), 
              show.legend = FALSE, se=FALSE, colour = "black")+
  ggpmisc::stat_poly_eq(formula = y ~ ns(x,2), parse = TRUE)
Plot_eub_east_west.p + labs(title = "Bacteria density in the water layers\n", color = "Region\n") +
  scale_colour_manual(values=c("EGC"="blue","WSC"="red"), labels= c("Ice covered", "Ice free"))+
  theme_plot

#Plot of total cells (DAPI)
data.dapi <-counts_DAPI_aggr2 
data.dapi$conc.mn <- log10(data.dapi$conc.mn)
data.dapi$smooth <- 1

Plot_dapi_east_west.p <- ggplot(na.omit(subset(data.dapi, Depth %in% c("DCM", "EPI", "MESO", "BATHY"))), aes(x = long, y = conc.mn, shape = Depth, group = Depth))+scale_shape_manual(values=c(19, 17, 15, 3), labels= c("Surface", "Epipelagic", "Mesopelagic", "Bathypelagic"))+  
  geom_point(aes(colour = Region), size =4)+
  xlab("Longitude [°East]")+
  ylab("Cell density [log10(cells/ml)]")+
  geom_text(aes(label = StationName), nudge_y= 0.08)+
  geom_smooth(aes(linetype = Depth), method = "gam", formula = y ~ ns(x,2), 
              show.legend = FALSE, se=FALSE, colour = "black")+
  ggpmisc::stat_poly_eq(formula = y ~ ns(x,2), parse = TRUE)
Plot_dapi_east_west.p + labs(title = "Total cell density in the water layers\n", color = "Region\n") +
  scale_colour_manual(values=c("EGC"="blue","WSC"="red"), labels= c("Ice covered", "Ice free"))+
  theme_plot

##################################
# Relative abundance
##################################
## Get proportions as proportions of DAPI
dapi.rel <- counts_DAPI
colnames(dapi.rel)[9] <- 'dapi.mn'
colnames(dapi.rel)[11] <- 'dapi.sd' 

fish.rel <- counts_FISH
colnames(fish.rel)[9] <- 'fish.mn'
colnames(fish.rel)[11] <- 'fish.sd' 

counts.rel.ab <- dapi.rel %>% left_join(fish.rel, by=c("StationName", "Domain", "Depth"))
counts.rel.ab <- counts.rel.ab[,c("StationName", "Depth", "Domain", "dapi.mn", "fish.mn")]
counts.rel.ab$proportion <- counts.rel.ab$fish.mn*100/counts.rel.ab$dapi.mn
counts.rel.ab$Region[counts.rel.ab$StationName %in% WSC] <- "WSC"
counts.rel.ab$Region[counts.rel.ab$StationName %in% EGC] <- "EGC"

##################################
# Wilcoxon test dcm FISH
##################################
w_fish_test_dcm <- counts_FISH[,c("conc.mn", "Domain", "Depth", "Region")]
w_fish_test_dcm <- subset(w_fish_test_dcm, Depth %in% c("DCM"))

w.EGC.dcm <- subset(w_fish_test_dcm, Region %in% c("EGC"))
w.EGC.dcm <-w.EGC.dcm[,c("Domain", "conc.mn")]
w.EGC.dcm.agg <- aggregate(.~Domain, data=w.EGC.dcm, FUN = mean)
w.EGC.dcm.agg <- w.EGC.dcm.agg[,c("conc.mn")]

w.WSC.dcm <- subset(w_fish_test_dcm, Region %in% c("WSC"))
w.WSC.dcm <-w.WSC.dcm[,c("Domain", "conc.mn")]
w.WSC.dcm.agg <- aggregate(.~Domain, data=w.WSC.dcm, FUN = mean)
w.WSC.dcm.agg <- w.WSC.dcm.agg[,c("conc.mn")]

w.test <- wilcox.test(w.EGC.dcm.agg, w.WSC.dcm.agg, paired = TRUE)
w.test

# Wilcoxon signed rank test
# 
# data:  w.EGC.dcm.agg and w.WSC.dcm.agg
# V = 31, p-value = 0.05768
# alternative hypothesis: true location shift is not equal to 0

##################################
# Wilcoxon test dcm dapi
##################################

w_dapi_test_dcm <- counts_DAPI[,c("conc.mn", "Domain", "Depth", "Region")]
w_dapi_test_dcm <- subset(w_dapi_test_dcm, Depth %in% c("DCM"))

w.d.EGC.dcm <- subset(w_dapi_test_dcm, Region %in% c("EGC"))
w.d.EGC.dcm <-w.d.EGC.dcm[,c("Domain", "conc.mn")]
w.d.EGC.dcm.agg <- aggregate(.~Domain, data=w.d.EGC.dcm, FUN = mean)
w.d.EGC.dcm.agg <- w.d.EGC.dcm.agg[,c("conc.mn")]

w.d.WSC.dcm <- subset(w_dapi_test_dcm, Region %in% c("WSC"))
w.d.WSC.dcm <-w.d.WSC.dcm[,c("Domain", "conc.mn")]
w.d.WSC.dcm.agg <- aggregate(.~Domain, data=w.d.WSC.dcm, FUN = mean)
w.d.WSC.dcm.agg <- w.d.WSC.dcm.agg[,c("conc.mn")]

w.d.test <- wilcox.test(w.d.EGC.dcm.agg, w.d.WSC.dcm.agg, paired = TRUE)
w.d.test


# Wilcoxon signed rank test with continuity correction
# 
# data:  w.d.EGC.dcm.agg and w.d.WSC.dcm.agg
# V = 0, p-value = 0.0004814
# alternative hypothesis: true location shift is not equal to 0


save.image("card_fish_water_column_ps99.Rdata")

