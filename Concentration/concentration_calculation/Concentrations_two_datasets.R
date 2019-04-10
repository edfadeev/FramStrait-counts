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
# WARNING: make sure to clean the environment, silence and re-run whatever is necesary for the statistical analyses
###################################
## import data set
###################################
#raw counts
raw.counts.SH <- read.csv("~/CARD-FISH/CARD-FISH_water_project/Concentration/single_hibr_counted/FOV_all_groups_SH.csv", sep = ",", dec = ".", header = TRUE)
#environmental data

metadata <- read.csv("~/CARD-FISH/PS99_samples_meta_EF_MC_nutrient_corrected.csv", sep = ",", dec = ".", header = TRUE)
###################################
#calculate cell concentration
###################################

#calculation factor
calc.factor <- 99515.5458411807

#cell concentration for DAPI
counts.aggregated_DAPI <- as.data.frame(as.list((aggregate(DAPI_Nr_Set ~SAMPLE_NAME, data=raw.counts.SH, 
                                                           FUN = function(x) c(mn = mean(x), md =  median(x), sd = sd(x), n = length(x))))))
counts.aggregated_DAPI.volume <- as.data.frame(as.list((aggregate(Volume ~SAMPLE_NAME, data=raw.counts.SH, 
                                                                  FUN = mean))))

counts.aggregated_DAPI <- cbind(counts.aggregated_DAPI,counts.aggregated_DAPI.volume[,2])
colnames(counts.aggregated_DAPI)[6] <- "volume"

counts.aggregated_DAPI$conc.mn <- (counts.aggregated_DAPI$DAPI_Nr_Set.mn*calc.factor)/counts.aggregated_DAPI$volume
counts.aggregated_DAPI$conc.md <- (counts.aggregated_DAPI$DAPI_Nr_Set.md*calc.factor)/counts.aggregated_DAPI$volume
counts.aggregated_DAPI$conc.sd <- (counts.aggregated_DAPI$DAPI_Nr_Set.sd*calc.factor)/counts.aggregated_DAPI$volume

counts.aggregated_DAPI <- counts.aggregated_DAPI %>% 
  separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-")

#cell concentration for FISH
counts.aggregated_FISH <- as.data.frame(as.list(aggregate(DAPI_Nr_SubSet ~SAMPLE_NAME, data=raw.counts.SH, 
                                                          FUN = function(x) c(mn = mean(x), md =  median(x), sd = sd(x), n = length(x)))))
counts.aggregated_FISH.volume <- as.data.frame(as.list((aggregate(Volume ~SAMPLE_NAME, data=raw.counts.SH, 
                                                                  FUN = mean))))

counts.aggregated_FISH <- cbind(counts.aggregated_FISH,counts.aggregated_FISH.volume[,2])
colnames(counts.aggregated_FISH)[6] <- "volume"

counts.aggregated_FISH$conc.mn <- (counts.aggregated_FISH$DAPI_Nr_SubSet.mn*calc.factor)/counts.aggregated_FISH$volume
counts.aggregated_FISH$conc.md <- (counts.aggregated_FISH$DAPI_Nr_SubSet.md*calc.factor)/counts.aggregated_FISH$volume
counts.aggregated_FISH$conc.sd <- (counts.aggregated_FISH$DAPI_Nr_SubSet.sd*calc.factor)/counts.aggregated_FISH$volume

counts.aggregated_FISH <- counts.aggregated_FISH %>% 
  separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-")

#call surface water layer and add regions

surf <- c("DCM")
surf2 <- c("DCM", "EPI")
depth <- c("MESO", "BATHY", "ABYSS")
EGC <- c("EG1","EG4")
HG <- c("HG9","HG7", "HG5", "HG4", "HG2","HG1")
N <- c("N3", "N4", "N5")
SV <- c("SV2")
WSC <- c("HG9","HG7","HG5","HG4","HG2","HG1","N3","N4","N5", "SV2")

counts.aggregated_DAPI$Region[counts.aggregated_DAPI$StationName %in% EGC] <- "EGC"
counts.aggregated_DAPI$Region[counts.aggregated_DAPI$StationName %in% HG] <- "HG"
counts.aggregated_DAPI$Region[counts.aggregated_DAPI$StationName %in% N] <- "N"
counts.aggregated_DAPI$Region[counts.aggregated_DAPI$StationName %in% SV] <- "SV"
counts.aggregated_DAPI$Region[counts.aggregated_DAPI$StationName %in% WSC] <- "WSC"

counts.aggregated_FISH$Region[counts.aggregated_FISH$StationName %in% EGC] <- "EGC"
counts.aggregated_FISH$Region[counts.aggregated_FISH$StationName %in% WSC] <- "WSC" # uncomment if necessary
counts.aggregated_FISH$Region[counts.aggregated_FISH$StationName %in% HG] <- "HG"
counts.aggregated_FISH$Region[counts.aggregated_FISH$StationName %in% N] <- "N"
counts.aggregated_FISH$Region[counts.aggregated_FISH$StationName %in% SV] <- "SV"

#add the longitude

counts.aggregated_FISH$long[counts.aggregated_FISH$StationName == "EG1"] <- -5.418
counts.aggregated_FISH$long[counts.aggregated_FISH$StationName == "EG4"] <- -2.729
counts.aggregated_FISH$long[counts.aggregated_FISH$StationName == "HG1"] <- 6.088
counts.aggregated_FISH$long[counts.aggregated_FISH$StationName == "HG2"] <- 4.907
counts.aggregated_FISH$long[counts.aggregated_FISH$StationName == "HG4"] <- 4.185
counts.aggregated_FISH$long[counts.aggregated_FISH$StationName == "HG5"] <- 3.8
counts.aggregated_FISH$long[counts.aggregated_FISH$StationName == "HG7"] <- 3.5
counts.aggregated_FISH$long[counts.aggregated_FISH$StationName == "HG9"] <- 2.841
counts.aggregated_FISH$long[counts.aggregated_FISH$StationName == "N4"] <- 4.508
counts.aggregated_FISH$long[counts.aggregated_FISH$StationName == "N3"] <- 5.166
counts.aggregated_FISH$long[counts.aggregated_FISH$StationName == "N5"] <- 3.062
counts.aggregated_FISH$long[counts.aggregated_FISH$StationName == "SV2"] <- 9.514


#Boxplots of cell concentration all groups all depths

counts.aggregated_FISH$Depth <- factor(counts.aggregated_FISH$Depth, levels = c("DCM","EPI","MESO","BATHY","ABYSS"))

Boxplot_FISH_counts <- ggplot(counts.aggregated_FISH[counts.aggregated_FISH$Depth %in% surf,], aes(Region, conc.mn))+
  #ggplot(counts.aggregated_FISH[counts.aggregated_FISH$Depth %in% surf,], aes(Region, conc.mn))+
  geom_boxplot()+facet_grid(Depth~.)+ geom_jitter()+ # add counts.aggregated_FISH[counts.aggregated_FISH$Depth %in% surf,] as data to select surface counts
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "cell conc. [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

#Boxplots of cell concentration DAPI all depths
counts.aggregated_DAPI$Depth <- factor(counts.aggregated_DAPI$Depth, levels = c("DCM","EPI","MESO","BATHY","ABYSS"))

Boxplot_DAPI_counts <- ggplot(counts.aggregated_DAPI[counts.aggregated_DAPI$Depth %in% surf,], aes(Region, conc.mn))+
  #ggplot(counts.aggregated_FISH[counts.aggregated_FISH$Depth %in% surf,], aes(Region, conc.mn))+
  geom_boxplot()+facet_grid(Depth~.)+ geom_jitter()+ # add counts.aggregated_FISH[counts.aggregated_FISH$Depth %in% surf,] as data to select surface counts
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "cell conc. [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

#scatter of cell concentration all groups all depths

counts.aggregated_FISH_aggr  <- aggregate(conc.mn~StationName+Region+Depth+long,counts.aggregated_FISH,mean)
counts.aggregated_FISH_aggr1  <- aggregate(conc.sd~StationName+Region+Depth+long,counts.aggregated_FISH,mean)
counts.aggregated_FISH_aggr2 <- merge(counts.aggregated_FISH_aggr1, counts.aggregated_FISH_aggr)

Scatter_FISH_counts <- ggplot(counts.aggregated_FISH_aggr2, aes(x=long, y=conc.mn)) + geom_point()+facet_grid(Depth~.)+
  geom_text(label=counts.aggregated_FISH_aggr2$StationName)+
  geom_errorbar(aes(ymin = conc.mn-conc.sd, ymax = conc.mn+conc.sd), width = 0.15)+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "CARD-FISH [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

#Plot of cell concentration all groups East vs West surface water

data.eub <- subset(counts.aggregated_FISH, Domain == "EUB")  # change EUB for other domains 
data.eub$conc.mn <- log10(data.eub$conc.mn)
data.eub$smooth <- 1

library(splines)
library(ggplot2)
library(gridExtra)
library(viridis)

Plot_east_west_EUB.p <- ggplot(na.omit(subset(data.eub, Depth %in% c("DCM", "MESO", "BATHY"))), aes(x = long, y = conc.mn, colour = Region, shape = Depth, group = Depth))+ 
  geom_point(size=5)+  
  xlab("Longitude [°East]")+
  ylab("log10(cells/ml")+
  geom_text(aes(label = StationName))+
  geom_smooth(aes(linetype = Depth), method = "gam", formula = y ~ ns(x,2), 
              show.legend = FALSE, se=FALSE, colour = "black")+
  scale_fill_manual(values=c("EGC"="blue","HG"="red"))+
  ggpmisc::stat_poly_eq(formula = y ~ ns(x,2), parse = TRUE)+
  theme_classic(base_size = 12)


##########################
# Environmental 
##########################

##Data standarization. By default, this function will standardize the data (mean zero, unit variance). To indicate that we just want to subtract the mean, we need to turn off the argument scale = FALSE.
env_raw <- metadata
#env_raw <- subset(env_raw, StationName!="HG1") #remove HG1 (NA) # silence to run adonis and permanova analyses

centre_scale2 <- function(x) {
   scale(x, scale = FALSE)
 }

env.par= c("Temperature","Salinity", "Chla_fluor")
env=env_raw
env[,env.par]<- apply(env[,env.par], MARGIN=2, FUN=centre_scale2)

## data correlation Temperature salinity by Pearson correlation test

library("ggpubr")
ggscatter(test, x = "Temperature", y = "Salinity", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Temperature", ylab = "Salinity")


res <- cor.test(env$Temperature, env$Salinity, 
                method = "pearson")
res

# normality test: data is not normally distributed p-values should be greater than the significance level 0.05

shapiro.test(env$Temperature) # p-value = 6.466e-06
shapiro.test(env$Salinity) # p-value = 1.135e-08
shapiro.test(env$Chla_fluor) # p-value = 0.0001543
shapiro.test(env$NO3.NO2) # p-value = 0.002956
shapiro.test(env$NO2) # p-value = 7.592e-05
shapiro.test(env$NO3) # p-value = 0.003924
shapiro.test(env$SiO3) # p-value = 0.001508
shapiro.test(env$PO4) # p-value = 0.02344
shapiro.test(env$NH4) # p-value = 2.357e-07

## non-parametric tests:

# Kendall rank correlation test

res2 <- cor.test(env$Temperature, env$Salinity,  method="kendall")
res2 # tau 0.3038223 

res2.1 <- cor.test(env$Temperature, env$Chla_fluor,  method="kendall")
res2.1 # tau -0.2577137

res2.2 <- cor.test(env$Salinity, env$Chla_fluor,  method="kendall")
res2.2 # tau -0.3012709

# Spearman rank correlation coefficient

res3 <-cor.test(env$Temperature, env$Salinity,  method = "spearman")
res3 # rho 0.4401271

res3.1 <-cor.test(env$Temperature, env$Chla_fluor,  method = "spearman")
res3.1 # rho -0.3400739 

res3.2 <-cor.test(env$Salinity, env$Chla_fluor,  method = "spearman")
res3.2 # rho -0.4157426  

# -1 indicates a strong negative correlation : this means that every time x increases, y decreases 
# 0 means that there is no association between the two variables (x and y) 
# 1 indicates a strong positive correlation : this means that y increases with x 

#Temperature and salinity are correlated to each other 
#Temperature and Chlorophyll a are negatively correlated to each other 
#Salinity and Chlorophyll a are negatively correlated to each other 

## non-parametric test for surface (DCM, 100M) without HG1:

surface <- c("DCM", "EPI")
env.surf <- subset(env[env$Depth %in% surface,])
#env.surf <- subset(env.surf, StationName!="HG1") # silence to run RDA and PCA
env.surf$Region[env.surf$StationName %in% EGC] <- "EGC"
env.surf$Region[env.surf$StationName %in% WSC] <- "WSC"
#env.surf$Region[env.surf$StationName %in% HG] <- "HG"
#env.surf$Region[env.surf$StationName %in% N] <- "N"
#env.surf$Region[env.surf$StationName %in% SV] <- "SV"

# Kendall rank correlation test: between physico-chemical parameters: 

res2.surf <- cor.test(env.surf$Temperature, env.surf$Salinity,  method="kendall")
res2.surf # tau 0.5844156 
res2.1.surf <- cor.test(env.surf$Temperature, env.surf$Chla_fluor,  method="kendall")
res2.1.surf # tau -0.2863347 
res2.2.surf <- cor.test(env.surf$Salinity, env.surf$Chla_fluor,  method="kendall")
res2.2.surf # tau -0.3557492

# Kendall rank correlation test: between physico-chemical parameters and nutrients: 

res2.3.surf <- cor.test(env.surf$Temperature, env.surf$NO3.NO2,  method="kendall") # here change only nutrient to see corrleation nutrient temperature
res2.3.surf # tau 0.4285714 
res2.4.surf <- cor.test(env.surf$Temperature, env.surf$NO3,  method="kendall") # here change only nutrient to see corrleation nutrient temperature
res2.4.surf # tau 0.4511941
res2.5.surf <- cor.test(env.surf$Temperature, env.surf$NO2,  method="kendall") # here change only nutrient to see corrleation nutrient temperature
res2.5.surf # tau 0.4289097 
res2.6.surf <- cor.test(env.surf$Temperature, env.surf$PO4,  method="kendall") # here change only nutrient to see corrleation nutrient temperature
res2.6.surf # tau 0.4095948 
res2.7.surf <- cor.test(env.surf$Temperature, env.surf$SiO3,  method="kendall") # here change only nutrient to see corrleation nutrient temperature
res2.7.surf # tau 0.2689811
res2.8.surf <- cor.test(env.surf$Temperature, env.surf$NH4,  method="kendall") # here change only nutrient to see corrleation nutrient temperature
res2.8.surf # tau 0.1413082 

res2.9.surf <- cor.test(env.surf$Chla_fluor, env.surf$NO3.NO2,  method="kendall") 
res2.9.surf # tau -0.4598709  
res2.10.surf <- cor.test(env.surf$Chla_fluor, env.surf$NO3,  method="kendall") 
res2.10.surf # tau -0.473913
res2.11.surf <- cor.test(env.surf$Chla_fluor, env.surf$NO2,  method="kendall")  
res2.11.surf # tau -0.2850987 
res2.12.surf <- cor.test(env.surf$Chla_fluor, env.surf$PO4,  method="kendall") 
res2.12.surf # tau -0.4061174 
res2.13.surf <- cor.test(env.surf$Chla_fluor, env.surf$SiO3,  method="kendall") 
res2.13.surf # tau -0.473913
res2.14.surf <- cor.test(env.surf$Chla_fluor, env.surf$NH4,  method="kendall") 
res2.14.surf # tau 0.1106368

# Spearman rank correlation coefficient between physico-chemical parameters: 

res3.surf <-cor.test(env.surf$Temperature, env.surf$Salinity,  method = "spearman")
res3.surf # rho 0.6871824
res3.1.surf <-cor.test(env.surf$Temperature, env.surf$Chla_fluor,  method = "spearman")
res3.1.surf # rho -0.3981926
res3.2.surf <-cor.test(env.surf$Salinity, env.surf$Chla_fluor,  method = "spearman")
res3.2.surf # rho -0.4975996

# Spearman rank correlation coefficient between physico-chemical parameters and nutrients: 

res2.3.surf <- cor.test(env.surf$Temperature, env.surf$NO3.NO2,  method = "spearman") # here change only nutrient to see corrleation nutrient temperature
res2.3.surf # rho 0.5855449  
res2.4.surf <- cor.test(env.surf$Temperature, env.surf$NO3,  method = "spearman") # here change only nutrient to see corrleation nutrient temperature
res2.4.surf # rho 0.5953121
res2.5.surf <- cor.test(env.surf$Temperature, env.surf$NO2,  method = "spearman") # here change only nutrient to see corrleation nutrient temperature
res2.5.surf # rho 0.5933266 
res2.6.surf <- cor.test(env.surf$Temperature, env.surf$PO4,  method = "spearman") # here change only nutrient to see corrleation nutrient temperature
res2.6.surf # rho 0.5613344 
res2.7.surf <- cor.test(env.surf$Temperature, env.surf$SiO3,  method = "spearman") # here change only nutrient to see corrleation nutrient temperature
res2.7.surf # rho 0.3767298
res2.8.surf <- cor.test(env.surf$Temperature, env.surf$NH4,  method = "spearman") # here change only nutrient to see corrleation nutrient temperature
res2.8.surf # rho 0.2244277 

res2.9.surf <- cor.test(env.surf$Chla_fluor, env.surf$NO3.NO2,  method = "spearman") 
res2.9.surf # rho -0.7031912   
res2.10.surf <- cor.test(env.surf$Chla_fluor, env.surf$NO3,  method = "spearman") 
res2.10.surf # rho -0.710452
res2.11.surf <- cor.test(env.surf$Chla_fluor, env.surf$NO2,  method = "spearman")  
res2.11.surf # rho -0.3719947 
res2.12.surf <- cor.test(env.surf$Chla_fluor, env.surf$PO4,  method = "spearman") 
res2.12.surf # rho -0.5744984 
res2.13.surf <- cor.test(env.surf$Chla_fluor, env.surf$SiO3,  method = "spearman") 
res2.13.surf # rho -0.6576271
res2.14.surf <- cor.test(env.surf$Chla_fluor, env.surf$NH4,  method = "spearman") 
res2.14.surf # rho 0.1496607


#Temperature and salinity are correlated to each other 
#Temperature and Chlorophyll a are negatively correlated to each other 
#Salinity and Chlorophyll a are negatively correlated to each other 
#Temperature and the nutrients are correlated to each other
#Chlorophyll and the nutrients are negatively correlated to each other 
#Temperature and NH4 are weakely correlated to each other
#Chlorophyll and NH4 are weakely correlated to each other

##########################
# Environmental surface and counts 
##########################

#RDA
library(tidyr)
env.surf.DCM.only <- subset(env.surf[env.surf$Depth %in% c("DCM"),])
counts.surf <- subset(counts.aggregated_FISH[counts.aggregated_FISH$Depth %in% c("DCM"),])
counts_data_for_stats.sub <- counts.surf[,c("Depth", "Domain", "StationName", "conc.mn")]
counts_data_for_stats.sub.agg <- aggregate(conc.mn~StationName+Domain+Depth, counts_data_for_stats.sub,mean)
counts_data_for_stats.sub.agg %>% spread(Domain, conc.mn) -> counts_data_for_stats.sub.agg

counts_data_for_stats.sub.agg <- subset(counts_data_for_stats.sub.agg, StationName!="HG1")
counts.surf.metadata <- env.surf.DCM.only %>% left_join(counts_data_for_stats.sub.agg, by=c("StationName"))
counts.surf.metadata$Domain <- NULL

counts.stats <- counts.surf.metadata[,c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")]
rownames(counts.stats) <- counts.surf.metadata$sample_ID
rownames(env.surf.DCM.only) <- env.surf.DCM.only$sample_ID
env.surf.DCM.only$sample_ID <- as.factor(env.surf.DCM.only$sample_ID)

ord <- rda(counts.stats ~ Temperature + Chla_fluor + SiO3 + NO3, data=env.surf.DCM.only, scale = TRUE, center = FALSE) # scale FALSE since it is already scaled

#generate dataframe for the samples
ps_scores <- vegan::scores(ord)
sites <- data.frame(ps_scores$sites)
sites$sample_ID <- as.factor(rownames(sites))

sites <- sites %>%
  left_join(env.surf.DCM.only)

############arrows - no need to touch, just run it########
# Now add the environmental variables as arrows
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
####################################

#plot
evals_prop <- 100 * ord$CCA$eig[1:2]
rda_plot <- ggplot() +
  geom_point(data = sites, aes(x = RDA1, y = RDA2, colour = Region),  size = 4) +
  geom_text(data = sites,aes(x = RDA1, y = RDA2,label = StationName), nudge_y= 0.1,size=3)+
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("\nAxis1 [%s%% variance]", round(evals_prop[1], 2)), 
       y = sprintf("Axis2 [%s%% variance]\n", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ord$CCA$eig[2] / ord$CCA$eig[1])) +
  theme_plot+ 
  geom_segment(mapping = arrow_mapFL, size = .5, data = arrowdfFL, color = "black", arrow = arrowhead) + 
  geom_text(mapping = label_mapFL, size = 5, data = arrowdfFL, show.legend = FALSE, nudge_y= 0.1)+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

# PCA attempt DCM
pca.data <- counts_data_for_stats.sub.agg

rownames(pca.data) <- pca.data$StationName
pca.data <- pca.data[,c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")]

pca.data.prcomp <- prcomp(pca.data, center = TRUE,scale. = TRUE)
summary(pca.data.prcomp)
str(pca.data.prcomp)
library(ggbiplot)

ggbiplot(pca.data.prcomp, labels=rownames(pca.data))

pca.data.Region <- c(rep("EGC", 2), rep("HG", 5), rep("N", 3), rep("SV", 1))

ggbiplot(pca.data.prcomp,ellipse=TRUE,  labels=rownames(pca.data), groups=pca.data.Region)+theme_plot

##########################
# Depth profiles cell concentration
##########################
 
data.deep <- counts.aggregated_FISH[,c("StationName", "Depth", "Domain", "conc.mn", "conc.sd", "Region")]
data.deep.agg <- aggregate(conc.mn~Domain+Depth+Region, data.deep,mean)
data.deep.agg1 <- aggregate(conc.sd~Domain+Depth+Region, data.deep,mean)
data.deep.m <- merge(data.deep.agg, data.deep.agg1, by = c("Domain", "Depth", "Region"))
data.deep.m$conc.mn <- log10(data.deep.m$conc.mn)

data.deep.m$Depth<- factor(data.deep.m$Depth, 
                        levels = rev(c("DCM", "EPI", "MESO", "BATHY", "ABYSS")))

Plot_depth_profiles.p <- ggplot()+
  geom_point(data = data.deep.m, aes(y = conc.mn, x = Depth, colour = Domain))+
  geom_line(data = data.deep.m, linetype=1, aes(y = conc.mn, x = Depth, colour = Domain, group = Domain),size = 0.8)+ # this geom_line is not connecting the depth right 
  coord_flip()+
  facet_grid(Region~Domain)+
  theme_plot

##
counts.odv <- counts.aggregated_FISH
counts.odv.sub <- counts.odv[,c("Depth", "Domain", "StationName", "conc.mn")]

counts_data_for_ODV <- aggregate(conc.mn~StationName+Depth, counts.odv.sub,mean)
env.for.odv <- env


counts_metadata_for_ODV <- env.for.odv %>% left_join(counts_data_for_ODV, by=c("StationName", "Depth"))
counts_metadata_for_ODV$Layer[counts_metadata_for_ODV$Depth == "DCM"] <- "10"
counts_metadata_for_ODV$Layer[counts_metadata_for_ODV$Depth == "EPI"] <- "100"
counts_metadata_for_ODV$Layer[counts_metadata_for_ODV$Depth == "MESO"] <- "1000"
counts_metadata_for_ODV$Layer[counts_metadata_for_ODV$Depth == "BATHY"] <- "2500"
counts_metadata_for_ODV$Layer[counts_metadata_for_ODV$Depth == "ABYSS"] <- "5500"
counts_metadata_for_ODV$conc.mn <- log10(counts_metadata_for_ODV$conc.mn)

counts_metadata_for_ODV <- counts_metadata_for_ODV[,c("StationID", "StationName", "sample_ID", "Latitude..degrees_north.", "Longitude..degrees_east.", "Depth", "Layer", "Depth.m.", "Bot_Depth..m.", "Temperature", "Salinity", "Chla_fluor", "NO3", "NO3.NO2", "NH4", "NO2", "conc.mn")]

write.csv(counts_metadata_for_ODV, "~/CARD-FISH/CARD-FISH_water_project/Concentration/counts_metadata_for_ODV.csv")

#############################################################
# Analyses with relative abundances (as a proportion of DAPI)
#############################################################

## Get proportions as a % of DAPI

dapi.count.rel <- counts.aggregated_DAPI
colnames(dapi.count.rel)[9] <- 'dapi.mn'
colnames(dapi.count.rel)[11] <- 'dapi.sd' 

fish.counts.rel <- counts.aggregated_FISH
colnames(fish.counts.rel)[9] <- 'fish.mn'
colnames(fish.counts.rel)[11] <- 'fish.sd' 

counts.rel.m <- dapi.count.rel %>% left_join(fish.counts.rel, by=c("StationName", "Domain", "Depth"))

counts.rel.m <- counts.rel.m[,c("StationName", "Depth", "Domain", "dapi.mn", "fish.mn")]

counts.rel.m$proportion <- counts.rel.m$fish.mn*100/counts.rel.m$dapi.mn

## PCA DCM  for Fish counts

counts.rel.pca <- subset(counts.rel.m[counts.rel.m$Depth %in% c("DCM"),])
counts.rel.pca <- counts.rel.pca[,c("Depth", "Domain", "StationName", "proportion")]
counts.rel.pca %>% spread(Domain, proportion) -> counts.rel.pca.agg

pca.rel <- counts.rel.pca.agg
rownames(pca.rel) <- pca.rel$StationName

pca.rel <- pca.rel[,c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")]

pca.re.data.prcomp <- prcomp(pca.rel, center = TRUE,scale. = TRUE)
summary(pca.re.data.prcomp)
str(pca.re.data.prcomp)
library(ggbiplot)

ggbiplot(pca.re.data.prcomp, labels=rownames(pca.rel))

pca.data.rel.Region <- c(rep("EGC", 2), rep("HG", 6), rep("N", 3), rep("SV", 1))

ggbiplot(pca.re.data.prcomp,ellipse=TRUE,  labels=rownames(pca.rel), groups=pca.data.rel.Region)+theme_plot

## RDA surface for Fish counts

counts.rel.ra <- subset(counts.rel.m[counts.rel.m$Depth %in% c("DCM", "EPI"),])
counts.rel.ra <- counts.rel.ra[,c("Depth", "Domain", "StationName", "proportion")]
counts.rel.ra %>% spread(Domain, proportion) -> counts.rel.ra.agg
counts.rel.ra.agg <- subset(counts.rel.ra.agg, StationName!="HG1")

env.surf.DCM.EPI.only <- subset(env.surf[env.surf$Depth %in% c("DCM", "EPI"),])

counts.rel.surf.metadata <- env.surf.DCM.EPI.only %>% left_join(counts.rel.ra.agg, by=c("StationName", "Depth"))
counts.rel.surf.metadata$Depth.x <- NULL

counts.rel.rda.stats <- counts.rel.surf.metadata[,c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")]
rownames(counts.rel.rda.stats) <- counts.rel.surf.metadata$sample_ID
rownames(env.surf.DCM.EPI.only) <- env.surf.DCM.EPI.only$sample_ID
env.surf.DCM.EPI.only$sample_ID <- as.factor(env.surf.DCM.EPI.only$sample_ID)

ord <- rda(counts.rel.rda.stats ~ Temperature + Chla_fluor + SiO3 + NO3, data=env.surf.DCM.EPI.only, scale = TRUE, center = FALSE) # scale FALSE since it is already scaled

#quartz()
plot(ord)

#generate dataframe for the samples
ps_scores <- vegan::scores(ord)
sites <- data.frame(ps_scores$sites)
sites$sample_ID <- as.factor(rownames(sites))

sites <- sites %>%
  left_join(env.surf.DCM.EPI.only)

############arrows - no need to touch, just run it########
# Now add the environmental variables as arrows
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
#######

#plot
#evals_prop <- 100 * ord$CCA$eig[1:2] # NOT CORRECT VALUE!!!
#use summary(ord) to add manually in the axis the "proportion explained" value
#rda_plot <- ggplot() +
#  geom_point(data = sites, aes(x = RDA1, y = RDA2, shape = Depth),  size = 4) +
#  geom_text(data = sites,aes(x = RDA1, y = RDA2,label = StationName), nudge_y= 0.1,size=3)+
#  guides(col = guide_legend(override.aes = list(size = 3))) +
#  labs(x = sprintf("\nAxis1 [%s%% variance]", round(evals_prop[1], 2)), 
#       y = sprintf("Axis2 [%s%% variance]\n", round(evals_prop[2], 2))) +
#  coord_fixed(sqrt(ord$CCA$eig[2] / ord$CCA$eig[1])) +
#  theme_plot+ 
#  geom_segment(mapping = arrow_mapFL, size = .5, data = arrowdfFL, color = "black", arrow = arrowhead) + 
#  geom_text(mapping = label_mapFL, size = 5, data = arrowdfFL, show.legend = FALSE, nudge_y= 0.1)+
#  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
#  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

## NMDS DCM for FISH counts

data.mnds <- pca.rel
#set.seed(201)
counts.mds <- metaMDS(data.mnds, maxit= 100, trace=TRUE)
str(counts.mds) # gives the stress
plot(counts.mds, type = "t")

data.scores <- as.data.frame(scores(counts.mds))  
data.scores$site <- rownames(data.scores) 
data.scores$Region[data.scores$site %in% EGC] <- "EGC"
data.scores$Region[data.scores$site %in% WSC] <- "WSC"
species.scores <- as.data.frame(scores(counts.mds, "species")) 
species.scores$species <- rownames(species.scores)


nmds_rel_ab <- ggplot()+
  #geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.30) # check if this is possible
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=Region,colour=Region),size=2) + 
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=4,vjust=-0.25)+
  #stat_ellipse(data=data.scores,aes(x=NMDS1,y=NMDS2, fill=Region), alpha=.2,type='t',size =1, geom="polygon")+ #too few points to calculate an ellipse
  scale_colour_manual(values=c("EGC" = "blue", "WSC" = "red")) +
  coord_equal() +
  theme_bw()+
    theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid label
        plot.background = element_blank())+
  annotate("text", x =-0.3, y = -0.25, label=paste('Stress =',round(counts.mds$stress,3)))

## ANOSIM water column with FISH counts

counts.anosim <- counts.rel.m[,c("Depth", "Domain", "StationName", "proportion")]
counts.anosim %>% spread(Domain, proportion) -> counts.anosim.sp

counts.anosim.with.meta <- env %>% left_join(counts.anosim.sp, by=c("StationName", "Depth"))
counts.anosim.with.meta$Region[counts.anosim.with.meta$StationName %in% EGC] <- "EGC"
counts.anosim.with.meta$Region[counts.anosim.with.meta$StationName %in% WSC] <- "WSC"

counts.rel.abs.dist <- dist(counts.anosim.with.meta[,c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")],"euclidean")

adonis(counts.rel.abs.dist ~  Water_mass + Region + Depth, counts.anosim.with.meta)

# Call:
#   adonis(formula = counts.rel.abs.dist ~ Water_mass + Region +      Depth, data = counts.anosim.with.meta) 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#    Df   SumsOfSqs MeanSqs       F.Model   R2    Pr(>F)    
# Water_mass  2   11707.4  5853.7 23.0479 0.49390  0.001 ***
# Region      1     179.6   179.6  0.7072 0.00758  0.524    
# Depth       3    1911.6   637.2  2.5088 0.08064  0.023 *  
# Residuals  39    9905.2   254.0         0.41787           
# Total      45   23703.8                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > 

## ANOSIM water column with DAPI counts

counts.anosim.d <- counts.rel.m[,c("Depth", "Domain", "StationName", "dapi.mn")]
counts.anosim.d %>% spread(Domain, dapi.mn) -> counts.anosim.d.sp

counts.anosim.d.with.meta <- env %>% left_join(counts.anosim.d.sp, by=c("StationName", "Depth"))
counts.anosim.d.with.meta$Region[counts.anosim.d.with.meta$StationName %in% EGC] <- "EGC"
counts.anosim.d.with.meta$Region[counts.anosim.d.with.meta$StationName %in% WSC] <- "WSC"

counts.rel.abs.d.dist <- dist(counts.anosim.d.with.meta[,c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")],"euclidean")

adonis(counts.rel.abs.d.dist ~  Water_mass + Region + Depth, counts.anosim.d.with.meta)

# Call:
#   adonis(formula = counts.rel.abs.d.dist ~ Water_mass + Region +      Depth, data = counts.anosim.d.with.meta) 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#      Df      SumsOfSqs    MeanSqs    F.Model  R2     Pr(>F)    
# Water_mass  2 1.6520e+14 8.2598e+13  55.833 0.44167  0.001 ***
# Region      1 1.9987e+13 1.9987e+13  13.511 0.05344  0.003 ** 
# Depth       3 1.3115e+14 4.3715e+13  29.550 0.35064  0.001 ***
# Residuals  39 5.7695e+13 1.4794e+12         0.15426           
# Total      45 3.7402e+14                    1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## TEST!!! Permanova for water column

# Normality test proportions in DCM: p-values should be greater than the significance level 0.05, some are normality distributed, some not

shapiro.test(counts.anosim.with.meta$ALT) # p-value = 7.149e-07 
shapiro.test(counts.anosim.with.meta$ARCH) # p-value = 0.01214
shapiro.test(counts.anosim.with.meta$BACT) # p-value = 1.185e-07
shapiro.test(counts.anosim.with.meta$CFX) # p-value = 0.01663
shapiro.test(counts.anosim.with.meta$CREN) # p-value = 1.272e-07
shapiro.test(counts.anosim.with.meta$DELTA) # p-value = 0.02678
shapiro.test(counts.anosim.with.meta$EUB) # p-value = 0.3208    <-
shapiro.test(counts.anosim.with.meta$GAM) # p-value = 1.81e-07
shapiro.test(counts.anosim.with.meta$OPI) # p-value = 2.495e-05
shapiro.test(counts.anosim.with.meta$POL) # p-value = 8.699e-09
shapiro.test(counts.anosim.with.meta$ROS) # p-value = 4.434e-06
shapiro.test(counts.anosim.with.meta$SAR11) # p-value = 0.005482
shapiro.test(counts.anosim.with.meta$SAR202) # p-value = 2.255e-05
shapiro.test(counts.anosim.with.meta$SAR324) # p-value = 0.0001698
shapiro.test(counts.anosim.with.meta$SAR406) # p-value = 7.316e-08
shapiro.test(counts.anosim.with.meta$VER) # p-value = 6.419e-07

# We assume anyway that the data is normally distributed since we trasnformed to proportions?

## perMANOVA 

adonis2(counts.rel.abs.dist ~ Water_mass + Depth + Region, counts.anosim.with.meta)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = counts.rel.abs.dist ~ Water_mass + Depth + Region, data = counts.anosim.with.meta)
# Df SumOfSqs      R2       F Pr(>F)    
# Water_mass  2  11707.4 0.49390 23.0479  0.001 ***
#   Depth       3   1873.1 0.07902  2.4584  0.018 *  
#   Region      1    218.1 0.00920  0.8586  0.439    
# Residual   39   9905.2 0.41787                   
# Total      45  23703.8 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## perMANOVA for EUB at DCM 

counts.anosim.eub <- subset(counts.rel.m[counts.rel.m$Domain %in% c("EUB"),])
counts.anosim.eub %>% spread(Domain, proportion) -> counts.anosim.eub


eub.anosim.with.meta <- env %>% left_join(counts.anosim.eub, by=c("StationName", "Depth"))
eub.anosim.with.meta$Region[eub.anosim.with.meta$StationName %in% EGC] <- "EGC"
eub.anosim.with.meta$Region[eub.anosim.with.meta$StationName %in% WSC] <- "WSC"
eub.anosim.with.meta <- subset(eub.anosim.with.meta[eub.anosim.with.meta$Depth %in% c("DCM"),])

counts.eub.rel.abs.dist <- dist(eub.anosim.with.meta[,c("EUB")],"euclidean")

adonis2(counts.eub.rel.abs.dist ~  Water_mass + Region, eub.anosim.with.meta)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = counts.eub.rel.abs.dist ~ Water_mass + Region, data = eub.anosim.with.meta)
# Df SumOfSqs      R2      F Pr(>F)
# Water_mass  1   208.59 0.15882 1.7084  0.227
# Region      1     5.92 0.00450 0.0485  0.825
# Residual    9  1098.88 0.83668              
# Total      11  1313.38 1.00000   

## perMANOVA for DAPI: note shapiro test shown not normal distribution

adonis2(counts.rel.abs.d.dist ~  Water_mass + Region + Depth, counts.anosim.d.with.meta)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = counts.rel.abs.d.dist ~ Water_mass + Region + Depth, data = counts.anosim.d.with.meta)
# Df   SumOfSqs      R2      F Pr(>F)    
# Water_mass  2 1.6520e+14 0.44167 55.833  0.001 ***
#   Region      1 1.9987e+13 0.05344 13.511  0.002 ** 
#   Depth       3 1.3115e+14 0.35064 29.550  0.001 ***
#   Residual   39 5.7695e+13 0.15426                  
# Total      45 3.7402e+14 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



