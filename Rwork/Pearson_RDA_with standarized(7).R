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
library(psych)
library(corrplot)
library(RColorBrewer)
library(PerformanceAnalytics)

##################################
# Env. parameter coorelation test with groups
##################################

#Pearson’s product moment correlation

#omitNA measurements

cor.test.2 <- counts.rda.meta[ , c(9, 10, 12, 13, 15, 17, 18, 20, 22, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39)]
                                      
cor.test.2 <-na.omit(cor.test.2[ , c(9, 10, 12, 13, 15, 17, 18, 20, 22, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39)])

#Calculate the Pearson’s correlation coefficient in a matrix
cor(cor.test.2)

pairs.panels(cor.test.2)

#Correlogram
x <- cor(cor.test.2)
corrplot(x, type="upper", order="hclust", col=brewer.pal(n=8, name="PuOr"), method="pie", addCoef.col = "black",srt = 45, number.cex= 0.75)

#####

counts_FISH_DCM <- subset(counts_FISH, Depth %in% c("DCM"))
counts_FISH_DCM<- counts_FISH_DCM[,c("Domain", "conc.mn", "StationName")]
counts_FISH_DCM %>% spread(Domain, conc.mn) -> counts_FISH_DCM_agg

env.DCM.sp<- subset(env, Depth %in% c("DCM"))

env.DCM.sp <- env.DCM.sp[,c("StationName","Temperature", "Salinity", "Chla_fluor", "NO3.NO2", "NO2", "NO3", "SiO3", "PO4", "NH4")]

counst_dcm_for_sp <- env.DCM.sp %>% left_join(counts_FISH_DCM_agg, by=c("StationName"))

counst_dcm_for_sp.par = c("ALT", "ARCH", "BACT","CFX","CREN","DELTA" ,"EUB" ,"GAM","OPI", "POL","ROS", "SAR11","SAR202","SAR324","SAR406", "VER")
counst_dcm_for_sp[,counst_dcm_for_sp.par]<- apply(counst_dcm_for_sp[,counst_dcm_for_sp.par], MARGIN=2, FUN=centre_scale2)


cor.test.3 <- subset(counst_dcm_for_sp, StationName!="HG1")
cor.test.3 <-cor.test.3[colnames(cor.test.3) != "StationName"]

#Calculate the Pearson’s correlation coefficient in a matrix
cor(cor.test.3)

pairs.panels(cor.test.3)

table <- as.data.frame(cor(cor.test.3))

write.csv(table, file="~/CARD-FISH/CARD-FISH_water_project/peason.csv")

#Correlogram
x <- cor(cor.test.3)
corrplot(x, type="upper", order="hclust", col=brewer.pal(n=8, name="PuOr"), method="pie", addCoef.col = "black",srt = 45, number.cex= 0.75)

##################################
#RDA with standarized nutrients
##################################
counts.rel.ab.surf <- subset(counts.rel.ab[counts.rel.ab$Depth %in% c("DCM", "EPI"),])
counts.rel.ab.surf.par = c("ALT", "ARCH", "BACT","CFX","CREN","DELTA" ,"EUB" ,"GAM","OPI", "POL","ROS", "SAR11","SAR202","SAR324","SAR406", "VER")
counts.rel.ab.surf[,counts.rel.ab.surf.par]<- apply(counts.rel.ab.surf[,counts.rel.ab.surf.par], MARGIN=2, FUN=centre_scale2)
counts.rda.surf <- counts.rel.ab.surf[,c("Domain", "StationName", "proportion", "Depth")]
counts.rda.surf %>% spread(Domain, proportion) -> counts.sp.rda
counts.rda <- counts.sp.rda[,c("StationName", "Depth", "ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")]
env.DCM <- subset(env[env$Depth %in% c("DCM", "EPI"),])
counts.rda.meta <- env.DCM %>% left_join(counts.rda, by=c("StationName", "Depth"))
rownames(counts.rda.meta) <- env.DCM$sample_ID
rownames(counts.rda) <- env.DCM$sample_ID 
rownames(env.DCM) <- env.DCM$sample_ID
counts.rda.meta$sample_ID <- as.factor(counts.rda.meta$sample_ID)
counts.rda.meta <- counts.rda.meta[colnames(counts.rda.meta) != "X"]
counts.rda.o <- counts.rda[colnames(counts.rda) != "StationName"]
counts.rda.ord <-counts.rda.o[colnames(counts.rda.o) != "Depth"]
counts.rda.meta$Region[counts.rda.meta$StationName %in% EGC] <- "EGC"
counts.rda.meta$Region[counts.rda.meta$StationName %in% WSC] <- "WSC"
counts.rda.meta[is.na(counts.rda.meta)] <- 0

#redundancy test
ord <- rda(counts.rda.ord ~ Temperature + Salinity + Chla_fluor + NH4 + SiO3 + PO4 + NO3, data=counts.rda.meta, scale = TRUE, center = FALSE)
plot(ord)
#generate dataframe for the samples
ps_scores <- vegan::scores(ord)
sites <- data.frame(ps_scores$sites)
sites$sample_ID <- as.factor(rownames(sites))
species <- data.frame(ps_scores$species)
species$sample_ID <- as.factor(rownames(species))
species$species <- rownames(species)


sites <- sites %>%
  left_join(counts.rda.meta)

# Adding environmental variables as arrows
arrowmatFL <- vegan::scores(ord, display = "bp")
# Add labels, make a data.frame
arrowdfFL <- data.frame(labels = rownames(arrowmatFL), arrowmatFL)
# Define the arrow aesthetic mapping
arrow_mapFL <- aes(xend = 1 *RDA1, 
                   yend =1 * RDA2, 
                   x = 0, 
                   y = 0, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)
label_mapFL <- aes(x = 1.1 * RDA1, 
                   y = 1.1 * RDA2, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))

# plot
evals_prop <- 100 * ord$CCA$eig[1:2] # NOT CORRECT VALUE!!!
#use summary(ord) to add manually in the axis the "proportion explained" value
rda_plot <- ggplot() +
  geom_point(data = sites, aes(x = RDA1, y = RDA2, colour = Region, shape=Depth),  size = 4) + scale_shape_manual(values =c(17,16),labels=c("Surface", "Epipelagic"))+
  geom_text(data=species,aes(x=RDA1,y=RDA2,label=species),alpha=0.5, size = 4)+
  geom_text(data = sites,aes(x = RDA1, y = RDA2,label = StationName), vjust = -1.2,size=4)+
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("\nAxis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]\n", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ord$CCA$eig[2] / ord$CCA$eig[1])) +
  theme_plot+
  geom_segment(mapping = arrow_mapFL, size = .5, data = arrowdfFL, color = "black", arrow = arrowhead) +
  geom_text(mapping = label_mapFL, size = 4, data = arrowdfFL, show.legend = FALSE, vjust = -0.9)+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'black', size= 0.5))
rda_plot + scale_colour_manual(values=c("EGC"="blue","WSC"="red"), labels = c("Ice-covered", "Ice-free"))+
  theme_plot
