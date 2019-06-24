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

###################################
## NMDS
###################################
FISH.ra<- counts_FISH[,c("StationName","Depth", "Domain", "Region")]
FISH.ra$FISH.ra <-counts_FISH$conc.mn/counts_DAPI$DAPI_conc.mn 

FISH.ra %>% spread(Domain,FISH.ra)-> FISH.ra.wide

#Calculate NMDS and add labels 
#list all taxa (excluding EUB)
taxa <- c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202","SAR324", "SAR406", "VER")
all_metaMDS <- metaMDS(FISH.ra.wide[,taxa], maxit= 100, trace=TRUE)

data.scores <- as.data.frame(scores(all_metaMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- FISH.ra.wide$StationName  # create a column of site names, from the rownames of data.scores
data.scores$Depth <- FISH.ra.wide$Depth  #  add the grp variable created earlier
data.scores$Region <- FISH.ra.wide$Region

species.scores <- as.data.frame(scores(all_metaMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

#Plot NMDS
ggplot() + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, colour = Region, shape= Depth),size=5) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=3,nudge_y =-0.1) +  # add the site labels
  scale_colour_manual(values=c("WSC" = "red", "EGC" = "blue")) +
  coord_equal() +
  theme_bw()


#Calculate NMDS of only DCM and EPI depths
FISH.ra.wide.SRF <- FISH.ra.wide[FISH.ra.wide$Depth %in% c("DCM","EPI"),]
all_metaMDS.SRF <- metaMDS(FISH.ra.wide.SRF[,taxa], maxit= 100, trace=TRUE)

data.scores.SRF <- as.data.frame(scores(all_metaMDS.SRF))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores.SRF$site <- FISH.ra.wide.SRF$StationName  # create a column of site names, from the rownames of data.scores
data.scores.SRF$Depth <- FISH.ra.wide.SRF$Depth  #  add the grp variable created earlier
data.scores.SRF$Region <- FISH.ra.wide.SRF$Region

species.scores.SRF <- as.data.frame(scores(all_metaMDS.SRF, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores.SRF$species <- rownames(species.scores.SRF)  # create a column of species, from the rownames of species.scores

#Plot NMDS
ggplot() + 
  geom_text(data=species.scores.SRF,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores.SRF,aes(x=NMDS1,y=NMDS2, colour = Region, shape= Depth),size=5) + # add the point markers
  geom_text(data=data.scores.SRF,aes(x=NMDS1,y=NMDS2,label=site),size=3,nudge_y =-0.1) +  # add the site labels
  scale_colour_manual(values=c("WSC" = "red", "EGC" = "blue")) +
  coord_equal() +
  theme_bw()


# #modify data for nmds
# counts.wm.rel <- counts.rel.ab[,c("Depth", "Domain", "StationName", "proportion", "Region")]
# counts.wm.rel %>% spread(Domain, proportion) -> counts.wm.all.depths
# counts.wm.wmetad <- counts.wm.all.depths %>% left_join(env, by=c("StationName", "Depth"))
# 
# c_wm_all <- counts.wm.wmetad[,c("sample_ID", "StationName", "Depth", "Region", "ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202","SAR324", "SAR406", "VER")]
# c_wm_all <- subset(c_wm_all, !Depth == "BATHY")
# c_wm_all <- subset(c_wm_all, !Depth == "ABYSS")
# c_wm_all <- subset(c_wm_all, !Depth == "MESO")
# rownames(c_wm_all) <- c_wm_all$sample_ID
# c_wm_nmds_all <- c_wm_all[,c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202","SAR324", "SAR406", "VER")]
# 
# #nmds all 
# c_wm_nmds_all_metaMDS <- metaMDS(c_wm_nmds_all, maxit= 100, trace=TRUE)
# str(c_wm_nmds_all_metaMDS) 
# plot(c_wm_nmds_all_metaMDS, type = "t")
# 
# #nmds data for plotting
# c_evn_all <- c_wm_all[,c("sample_ID", "StationName", "Depth", "Region")] 
# #c_evn_all$Depth<- factor(c_evn_all$Depth, 
# #                               levels = c("DCM", "EPI"))
# c_evn_all$NMDS1<-c_wm_nmds_all_metaMDS$points[ ,1] 
# c_evn_all$NMDS2<-c_wm_nmds_all_metaMDS$points[ ,2] 
# data.scores.all <- as.data.frame(scores(c_wm_nmds_all_metaMDS))  
# data.scores.all$sample_ID <- rownames(data.scores.all) 
# species.scores.all <- as.data.frame(scores(c_wm_nmds_all_metaMDS, "species")) 
# species.scores.all$species <- rownames(species.scores.all)

# function for ellipses
#taken from the excellent stackoverflow Q+A: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the regions:ice-covered and ice-free
wm.df.all <- data.frame() 
for(g in levels(c_evn_all$Region)){
  wm.df.all <- rbind(wm.df.all, cbind(as.data.frame(with(c_evn_all [c_evn_all$Region==g,],
                                                 veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                              ,Region=g))
}

# data for labelling the ellipse
NMDS.mean.cwm.all=aggregate(c_evn_all[ ,c("NMDS1", "NMDS2")], 
                        list(group = c_evn_all$Region), mean)

# data for labelling the ellipse
NMDS.mean.all=aggregate(c_evn_all[,c("NMDS1", "NMDS2")], 
                    list(group = c_evn_all$Region), mean)
#plot
nmds_all_wm <- ggplot()+
  #geom_path(data = wm.df.all, aes(x = NMDS1, y = NMDS2, color = Region))+ 
  geom_text(data=species.scores.all,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size = 5) + 
  geom_point(data=c_evn_all,aes(x=NMDS1,y=NMDS2,shape=Depth,colour=Region),size=4) +
  stat_ellipse(type = "t")+
  scale_shape_manual(values=c(19, 17), labels= c("Surface", "Epipelagic"))+ 
  geom_text(data=c_evn_all,aes(x=NMDS1,y=NMDS2,label=StationName),size=4,vjust=-0.8, hjust=0.8)+
  scale_colour_manual(values=c("EGC" = "blue", "WSC" = "red"), labels =c("Ice covered", "Ice free")) +
  theme_plot+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())+
  annotate("text", x =-0.4, y = -0.29, label=paste('Stress =',round(c_wm_nmds_all_metaMDS$stress,3)))
nmds_all_wm

###################################
## PCA DCM and EPI
###################################
pca.all.data.prcomp <- prcomp(c_wm_nmds_all, center = TRUE,scale. = TRUE)
summary(pca.all.data.prcomp)
str(pca.all.data.prcomp)
ggbiplot(pca.all.data.prcomp, labels=rownames(c_wm_nmds_all))

pca.data.rel.Region <- c(rep("Ice-covered", 5), rep("Ice-free", 6))

PCA_plot_rel_all <- ggbiplot(pca.all.data.prcomp,ellipse=TRUE,  labels=rownames(c_wm_nmds_all), groups=pca.data.rel.Region)+theme_plot
PCA_plot_rel_all + scale_colour_manual(values=c("Ice-covered"="blue","Ice-free"="red"))+
  theme_plot

###################################
## PCoA
###################################
#Evaluate the distances
EuclDistAT <- dist(c_wm_nmds_all, method = "euclidean")
PcoA2D_AT <- cmdscale(EuclDistAT, k =2, eig=TRUE)
PcoA2D_AT$eig
PcoA2D_AT$eig[1:2]/sum(PcoA2D_AT$eig)*100

explainvar1 <- round(PcoA2D_AT$eig[1] / sum(PcoA2D_AT$eig), 2) * 100
explainvar2 <- round(PcoA2D_AT$eig[2] / sum(PcoA2D_AT$eig), 2) * 100

pcoa.scores.all <- as.data.frame(scores(PcoA2D_AT))  
pcoa.species.scores.all <- as.data.frame(wascores(PcoA2D_AT$points[,1:2], c_wm_nmds_all))

names(pcoa.species.scores.all)[1:2] <- c('PC1', 'PC2')
names(pcoa.scores.all)[1:2] <- c('PC1', 'PC2')

env_pcoa <- env
env_pcoa <- subset(env_pcoa, !Depth == "BATHY")
env_pcoa <- subset(env_pcoa, !Depth == "ABYSS")
env_pcoa <- subset(env_pcoa, !Depth == "MESO")
env_pcoa$Region[env_pcoa$StationName %in% EGC] <- "EGC"
env_pcoa$Region[env_pcoa$StationName %in% WSC] <- "WSC"

PcoA2D_AT$Region <- env_pcoa$Region
PcoA2D_AT$Depth <- env_pcoa$Depth
PcoA2D_AT$StationName <- env_pcoa$StationName
pcoa.species.scores.all$species <- rownames(pcoa.species.scores.all)
pcoa.scores.all$sample_ID <- rownames(pcoa.scores.all)  
pcoa.scores.all$Region <- env_pcoa$Region
pcoa.scores.all$Depth <- env_pcoa$Depth
pcoa.scores.all$Depth <- factor(pcoa.scores.all$Depth, levels = c("DCM","EPI"))
pcoa.scores.all$StationName <- env_pcoa$StationName

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

#plot
Tr_PcoA2 <- ggplot()+
  geom_text(data=pcoa.species.scores.all,aes(x=PC1,y=PC2,label=species), alpha=0.5,size = 3)+ 
  geom_point(data=pcoa.scores.all,aes(x=PC1,y=PC2,shape=Depth, color = Region),size=4) +
  scale_shape_manual(values=c(19, 17), labels= c("Surface", "Epipelagic"))+
  geom_text(data=pcoa.scores.all,aes(x=PC1,y=PC2,label=StationName),size=4,vjust=-0.8, hjust=0.8)+
  scale_colour_manual(values=c("EGC" = "blue", "WSC" = "red"), labels =c("Ice covered", "Ice free")) +
  #coord_equal() +
  xlab("PCoA1 [60 %]")+
  ylab("PCoA2 [25 %]")+
  stat_ellipse(type = "t")+
  theme_plot+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())
Tr_PcoA2 

save.image("card_fish_water_column_ps99.Rdata")
