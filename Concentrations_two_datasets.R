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
EGC <- c("EG1","EG4")
HG <- c("HG9","HG7", "HG5", "HG4", "HG2","HG1")
N <- c("N3", "N4", "N5")
SV <- c("SV2")
WSC <- c("HG9","HG7","HG5","HG4","HG2","HG1","N3","N4","N5", "SV2")

counts.aggregated_DAPI$Region[counts.aggregated_DAPI$StationName %in% EGC] <- "EGC"
counts.aggregated_DAPI$Region[counts.aggregated_DAPI$StationName %in% HG] <- "HG"
counts.aggregated_DAPI$Region[counts.aggregated_DAPI$StationName %in% N] <- "N"
counts.aggregated_DAPI$Region[counts.aggregated_DAPI$StationName %in% SV] <- "SV"

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

#Boxplots of cell concentration all groups all depths

counts.aggregated_FISH$Depth <- factor(counts.aggregated_FISH$Depth, levels = c("DCM","EPI","MESO","BATHY","ABYSS"))

Boxplot_FISH_counts <- ggplot(counts.aggregated_FISH, aes(Region, conc.mn))+geom_boxplot()+facet_grid(Depth~.)+ # add counts.aggregated_FISH[counts.aggregated_FISH$Depth %in% surf,] as data to select surface counts
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

data.eub <- subset(counts.aggregated_FISH[counts.aggregated_FISH$Depth %in% surf,], Domain == "ALT") # change EUB for other domains 
data.eub$conc.mn <- log10(data.eub$conc.mn)
data.eub$smooth <- 1

library(splines)
library(ggplot2)
library(gridExtra)
library(viridis)

Plot_east_west_ARCH.p <- ggplot(na.omit(data.eub), aes(x = long, y = conc.mn, fill = Region, group = smooth))+ 
  geom_point(shape=21, size=5, color = "black")+
  xlab("Longitude [°East]")+
  ylab("log10(cells/ml")+
  geom_smooth(method = "gam", formula = y ~ ns(x,2), show.legend = FALSE, se=FALSE, colour = "black")+
  scale_fill_manual(values=c("EGC"="blue","HG"="red"))+
  ggpmisc::stat_poly_eq(formula = y ~ ns(x,2), parse = TRUE)+
  theme_classic(base_size = 12)

##########################
# Environmental 
##########################

##Data standarization. By default, this function will standardize the data (mean zero, unit variance). To indicate that we just want to subtract the mean, we need to turn off the argument scale = FALSE.
env <- metadata

centre_scale2 <- function(x) {
   scale(x, scale = FALSE)
 }

env$NH4 <- centre_scale2(env$NH4) #run column by column by replacing in env


## data correlation Temperature salinity by Pearson correlation test

library("ggpubr")
ggscatter(env, x = "Temperature", y = "Salinity", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Temperature", ylab = "Salinity")


res <- cor.test(env$Temperature, env$Salinity, 
                method = "pearson")
res

# normality test: data is not normally distributed p-values should be greater than the significance level 0.05

shapiro.test(env$Temperature) # p-value = 1.806e-05
shapiro.test(env$Salinity) # p-value = 3.823e-09
shapiro.test(env$Chla_fluor) # p-value = 0.0001371


## non-parametric tests:

# Kendall rank correlation test

res2 <- cor.test(env$Temperature, env$Salinity,  method="kendall")
res2 # tau 0.3038223 

res2.1 <- cor.test(env$Temperature, env$Chla_fluor,  method="kendall")
res2.1 # tau -0.2577137

res2.2 <- cor.test(env$Salinity, env$Chla_fluor,  method="kendall")
res2.2 # tau -0.3012709

#Spearman rank correlation coefficient

res3 <-cor.test(env$Temperature, env$Salinity,  method = "spearman")
res3 # rho 0.4401271

res3.1 <-cor.test(env$Temperature, env$Chla_fluor,  method = "spearman")
res3.1 # rho -0.3400739 

res3.2 <-cor.test(env$Salinity, env$Chla_fluor,  method = "spearman")
res3.2 # rho -0.4157426  

# -1 indicates a strong negative correlation : this means that every time x increases, y decreases (left panel figure)
# 0 means that there is no association between the two variables (x and y) (middle panel figure)
# 1 indicates a strong positive correlation : this means that y increases with x (right panel figure)

#Temperature and salinity are correlated to each other 
#Temperature and Chlorophyll a are negatively correlated to each other 
#Salinity and Chlorophyll a are negatively correlated to each other 

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
  geom_point(data = data.deep.m, aes(y = Depth, x = conc.mn, colour = Domain))+
  geom_line(data = data.deep.m[data.deep.m$Region == "EGC",], linetype=1, aes(y = Depth, x = conc.mn, colour = Domain, group = Domain),size = 0.8)+ # this geom_line is not connecting the depth right 
  geom_line(data = data.deep.m[data.deep.m$Region == "WSC",], linetype=1, aes(y = Depth, x = conc.mn, colour = Domain, group = Domain),size = 0.8)+
  #scale_x_reverse(breaks = c(50,100,1000,2500,4000,5500))+
  #coord_flip()+
  #scale_color_manual(values = c("EGC" = "blue","WSC"= "red"))+ #uncommented to get default color for more groups
  facet_grid(Region~Domain)+
  theme_plot

