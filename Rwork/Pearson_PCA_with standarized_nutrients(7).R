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
library(factoextra)

##################################
# Pearson’s product moment correlation with rela. abundance
##################################
#DCM and EPI:

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
#DCM only and standarizing cell counts and nutrients:

#Data preparation:
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
#PCA environmental variables and taxonomic groups (DCM)
##################################
#rel. abund. 
test1_pca <- counts.rda.meta[,c( "StationName", "Depth", "Temperature", "Salinity", "Chla_fluor", "NO3.NO2", "NO2", "NO3", "SiO3", "PO4", "NH4", "ALT",  "ARCH",  "BACT",  "CFX",  "CREN",  "DELTA","EUB", "GAM" ,"OPI","POL","ROS","SAR11","SAR202","SAR324","SAR406", "VER", "Region")]
test1_pca <- subset(test1_pca, Depth == "DCM")
test1_pca <- subset(test1_pca, !StationName == "HG1")
rownames(test1_pca) <- test1_pca$StationName

test2_pca <- test1_pca[,c("Temperature", "Salinity", "Chla_fluor", "NO3.NO2", "NO2", "NO3", "SiO3", "PO4", "NH4","ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")]

names(test2_pca)[names(test2_pca) == "Chla_fluor"] <- "Chl a"
names(test2_pca)[names(test2_pca) == "NO3.NO2"] <- "Total N"

test2_pca.prcomp <- prcomp(test2_pca, center = TRUE,scale = TRUE)
summary(test2_pca.prcomp)
Region <- as.factor(test2_abs_pca$Region)

autoplot(test2_pca.prcomp, labels=rownames(test2_pca))

#plot rel. abun 
PCA_DCM_rel <- fviz_pca_biplot(test2_pca.prcomp,
                repel = TRUE,
                col.var = c("#696969"),
                col.ind = Region,
                palette = c("blue",  "red"), habillage = test2_abs_pca$Region, title = "")

####Absolute counts
test1_abs_pca <- subset(counts_FISH, Depth == "DCM")
env.DCM.pca <- subset(env.DCM, Depth == "DCM")
test1_abs_pca <- subset(test1_abs_pca, !StationName == "HG1")
env.DCM.pca<- subset(env.DCM.pca, !StationName == "HG1")
env.DCM.pca <- env.DCM.pca[,c("StationName", "Temperature", "Salinity", "Chla_fluor", "NO3.NO2", "NO2", "NO3", "SiO3", "PO4", "NH4")]

test1_abs_pca <-test1_abs_pca[,c("StationName", "Domain", "Region",  "conc.mn")]
test1_abs_pca %>% spread(Domain, conc.mn) -> test1_abs_pca_sp

test2_abs_pca <- test1_abs_pca_sp %>% left_join(env.DCM.pca, by=c("StationName")) 


rownames(test2_abs_pca) <- test2_abs_pca$StationName

test3_abs_pca <- test2_abs_pca[,c("Temperature", "Salinity", "Chla_fluor", "NO3.NO2", "NO2", "NO3", "SiO3", "PO4", "NH4","ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")]

names(test3_abs_pca)[names(test3_abs_pca) == "Chla_fluor"] <- "Chl a"
names(test3_abs_pca)[names(test3_abs_pca) == "NO3.NO2"] <- "Total N"

test2_abs_pca.prcomp <- prcomp(test3_abs_pca, center = TRUE,scale = TRUE)
summary(test2_abs_pca.prcomp)

Region <- as.factor(test2_abs_pca$Region)
test2_abs_pca$Region[10] <- "Ice-covered"

autoplot(test2_abs_pca.prcomp, labels=rownames(test2_abs_pca))

#plot 
PCA_DCM_abs <- fviz_pca_biplot(test2_abs_pca.prcomp,
                               repel = TRUE,
                               col.var = c("#696969"),
                               col.ind = Region,
                               palette = c("blue",  "red"),habillage = test2_abs_pca$Region,  title = "")

pca_both<- cowplot::plot_grid(PCA_DCM_abs, PCA_DCM_rel, 
                       labels = c("A", "B"),
                       ncol = 1, nrow = 2)

