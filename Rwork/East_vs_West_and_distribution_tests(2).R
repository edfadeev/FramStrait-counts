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


#Boxplots of cell concentration DCM:
#FISH
Boxplot_FISH_counts <- ggplot(subset(counts_FISH, Depth %in% c("DCM")), aes(Region, conc.mn))+
  geom_boxplot()+facet_grid(Depth~.)+ geom_text(aes(label=Domain), size=3, position = position_jitter(w = 0.3))+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "cell conc. [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

Boxplot_FISH_counts_all <- ggplot(counts_FISH, aes(Region, conc.mn))+
  geom_boxplot()+facet_grid(Depth~.)+ geom_text(aes(label=Domain), size=3, position = position_jitter(w = 0.3))+
  theme(axis.text.x = element_text(angle = 90))+theme_plot+scale_y_log10(name = "cell conc. [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))


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
#Plot of bacterial (FISH) cell concentration all groups East vs West 
data.fish <-counts_FISH_aggr2 
data.fish$conc.mn <- log10(data.fish$conc.mn)
data.fish$smooth <- 1

Plot_fish_east_west.p <- ggplot(na.omit(subset(data.fish, Depth %in% c("DCM", "EPI", "MESO", "BATHY"))), aes(x = long, y = conc.mn, shape = Depth, group = Depth))+  
  geom_point(aes(colour = Region), size =4)+
  xlab("Longitude [°East]")+
  ylab("log10(cells/ml")+
  geom_text(aes(label = StationName), nudge_y= 0.08)+
  geom_smooth(aes(linetype = Depth), method = "gam", formula = y ~ ns(x,2), 
              show.legend = FALSE, se=FALSE, colour = "black")+
  ggpmisc::stat_poly_eq(formula = y ~ ns(x,2), parse = TRUE)
Plot_fish_east_west.p + scale_colour_manual(values=c("EGC"="blue","WSC"="red"))+
  theme_plot

#Plot of total cells (DAPI)
data.dapi <-counts_DAPI_aggr2 
data.dapi$conc.mn <- log10(data.dapi$conc.mn)
data.dapi$smooth <- 1

Plot_dapi_east_west.p <- ggplot(na.omit(subset(data.dapi, Depth %in% c("DCM", "EPI", "MESO", "BATHY"))), aes(x = long, y = conc.mn, shape = Depth, group = Depth))+  
  geom_point(aes(colour = Region), size =4)+
  xlab("Longitude [°East]")+
  ylab("log10(cells/ml")+
  geom_text(aes(label = StationName), nudge_y= 0.08)+
  geom_smooth(aes(linetype = Depth), method = "gam", formula = y ~ ns(x,2), 
              show.legend = FALSE, se=FALSE, colour = "black")+
  ggpmisc::stat_poly_eq(formula = y ~ ns(x,2), parse = TRUE)
Plot_dapi_east_west.p + scale_colour_manual(values=c("EGC"="blue","WSC"="red"))+
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

##################################
# NMDS for surface layers
##################################
#Manipulate data for nmds at DCM
counts.rel.ab.dcm <- subset(counts.rel.ab[counts.rel.ab$Depth %in% c("DCM"),])
counts.rel.ab.sub.dcm <- counts.rel.ab.dcm[,c("Depth", "Domain", "StationName", "proportion")]
counts.rel.ab.sub.dcm %>% spread(Domain, proportion) -> counts.rel.ab.agg
counts_nmds <- counts.rel.ab.agg
rownames(counts_nmds) <- counts_nmds$StationName
counts_nmds <- counts_nmds[,c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")]

set.seed(201)
counts_nmds <- metaMDS(counts_nmds, maxit= 100, trace=TRUE)
str(counts_nmds) 
#plot(counts_nmds, type = "t")

data.scores <- as.data.frame(scores(counts_nmds))  
data.scores$site <- rownames(data.scores) 
data.scores$Region[data.scores$site %in% EGC] <- "EGC"
data.scores$Region[data.scores$site %in% WSC] <- "WSC"
species.scores <- as.data.frame(scores(counts_nmds, "species")) 
species.scores$species <- rownames(species.scores)

nmds_rel_ab_dcm <- ggplot()+
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size = 5) + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=Region,colour=Region),size=4) + 
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=5,vjust=-0.25)+
  scale_colour_manual(values=c("EGC" = "blue", "WSC" = "red")) +
  coord_equal() +
  theme_plot+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())+
  annotate("text", x =-0.3, y = -0.25, label=paste('Stress =',round(counts_nmds$stress,3)))

#Manipulate data for nmds at EPI
counts.rel.ab.epi <- subset(counts.rel.ab[counts.rel.ab$Depth %in% c("EPI"),])
counts.rel.ab.sub.epi <- counts.rel.ab.epi[,c("Depth", "Domain", "StationName", "proportion")]
counts.rel.ab.sub.epi %>% spread(Domain, proportion) -> counts.rel.ab.agg.epi
counts_nmds_epi <- counts.rel.ab.agg.epi
rownames(counts_nmds_epi) <- counts_nmds_epi$StationName
counts_nmds_epi <- counts_nmds_epi[,c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")]

#set.seed(201)
counts_nmds_epi <- metaMDS(counts_nmds_epi, maxit= 100, trace=TRUE)
str(counts_nmds_epi) 
plot(counts_nmds_epi, type = "t")

data.scores.epi <- as.data.frame(scores(counts_nmds_epi))  
data.scores.epi$site <- rownames(data.scores.epi) 
data.scores.epi$Region[data.scores.epi$site %in% EGC] <- "EGC" #data.scores$Water_mass[data.scores$site %in% AW] <- "AW"
data.scores.epi$Region[data.scores.epi$site %in% WSC] <- "WSC" #data.scores$Water_mass[data.scores$site %in% PSWw] <- "PSWw"
species.scores.epi <- as.data.frame(scores(counts_nmds_epi, "species")) 
species.scores.epi$species <- rownames(species.scores.epi)

nmds_rel_ab_epi <- ggplot()+
  geom_text(data=species.scores.epi,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size = 5) + 
  geom_point(data=data.scores.epi,aes(x=NMDS1,y=NMDS2,shape=Region,colour=Region),size=4) + 
  geom_text(data=data.scores.epi,aes(x=NMDS1,y=NMDS2,label=site),size=4,vjust=-0.25)+
  scale_colour_manual(values=c("EGC" = "blue", "WSC" = "red")) + #("PSWw" = "blue", "AW" = "red")
  coord_equal() +
  theme_plot+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())+
  annotate("text", x =-0.1, y = -0.1, label=paste('Stress =',round(counts_nmds_epi$stress,3)))

##################################
# PCA for surface layers
##################################
pca.rel <- counts.rel.ab.agg
rownames(pca.rel) <- pca.rel$StationName
pca.rel <- pca.rel[,c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")]

pca.re.data.prcomp <- prcomp(pca.rel, center = TRUE,scale. = TRUE)
summary(pca.re.data.prcomp)
str(pca.re.data.prcomp)

ggbiplot(pca.re.data.prcomp, labels=rownames(pca.rel))
pca.data.rel.Region <- c(rep("EGC", 2), rep("WSC", 10))

PCA_plot_rel <- ggbiplot(pca.re.data.prcomp,ellipse=TRUE,  labels=rownames(pca.rel), groups=pca.data.rel.Region)+theme_plot
PCA_plot_rel + scale_colour_manual(values=c("EGC"="blue","WSC"="red"))+
  theme_plot

##################################
# ANOVA (cell concentration is exlpained by water masses and depth over region)
##################################
# data preparation
aov.counts <- data.fish
aov.counts <- aov.counts %>% left_join(env, by=c("StationName", "Depth"))

aov.counts$Region <- as.factor(aov.counts$Region)
aov.counts$Depth <- as.factor(aov.counts$Depth)
aov.counts <- aov.counts[,c("Region", "Depth", "Water_mass","conc.mn")]

#balance/unbalance:
!is.list(replications(conc.mn ~ Region + Depth , aov.counts)) # FALSE

#nested anova
anov.counts1 <-aov(conc.mn ~ Region/Depth, aov.counts)
summary(anov.counts1)
plot(resid(anov.counts1) ~ fitted(anov.counts1))

# Df Sum Sq Mean Sq F value Pr(>F)    
# Region        1  0.104  0.1038   2.333  0.135    
# Region:Depth  7 19.066  2.7237  61.205 <2e-16 ***
#   Residuals    37  1.647  0.0445                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#balance/unbalance:
!is.list(replications(conc.mn ~ Water_mass + Depth, aov.counts)) # FALSE

#nested anova
anov.counts2 <-aov(conc.mn ~ Water_mass/Depth, aov.counts)
summary(anov.counts2)
plot(resid(anov.counts2) ~ fitted(anov.counts2))

# Df Sum Sq Mean Sq F value   Pr(>F)    
# Water_mass        2 17.250   8.625  192.28  < 2e-16 ***
#   Water_mass:Depth  3  1.772   0.591   13.17 4.01e-06 ***
#   Residuals        40  1.794   0.045                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##################################
# ADONIS (Permanova)
##################################
#  Non-parametric multivariate analysis of variance
#  Based on any dissimilarity measure

counts_abs.dist <- dist(aov.counts[,c("conc.mn")],"euclidean")

adon.counts <- adonis(counts_abs.dist ~  Water_mass + Region + Depth, aov.counts)

# Call:
#   adonis(formula = counts_abs.dist ~ Water_mass + Region + Depth,      data = aov.counts) 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Water_mass  2   17.2499  8.6250 213.616 0.82868  0.001 ***
#   Region      1    0.2984  0.2984   7.390 0.01433  0.006 ** 
#   Depth       3    1.6931  0.5644  13.977 0.08133  0.001 ***
#   Residuals  39    1.5747  0.0404         0.07565           
# Total      45   20.8160                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



save.image("card_fish_water_column_ps99.Rdata")

