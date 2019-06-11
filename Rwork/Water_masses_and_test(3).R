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

#modify data for nmds
counts.wm.rel <- counts.rel.ab[,c("Depth", "Domain", "StationName", "proportion")]
counts.wm.rel %>% spread(Domain, proportion) -> counts.wm.all.depths
counts.wm.wmetad <- counts.wm.all.depths %>% left_join(env, by=c("StationName", "Depth"))

c_wm_all <- counts.wm.wmetad[,c("sample_ID", "StationName", "Depth", "Water_mass", "ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202","SAR324", "SAR406", "VER")]
rownames(c_wm_all) <- c_wm_all$sample_ID
c_wm_nmds_all <- c_wm_all[,c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202","SAR324", "SAR406", "VER")]

#nmds all 
c_wm_nmds_all_metaMDS <- metaMDS(c_wm_nmds_all, maxit= 100, trace=TRUE)
str(c_wm_nmds_all_metaMDS) 
plot(c_wm_nmds_all_metaMDS, type = "t")

#nmds data for plotting
c_evn_all <- c_wm_all[,c("sample_ID", "StationName", "Depth", "Water_mass")] 
c_evn_all$Depth<- factor(c_evn_all$Depth, 
                               levels = c("DCM", "EPI", "MESO", "BATHY", "ABYSS"))
c_evn_all$NMDS1<-c_wm_nmds_all_metaMDS$points[ ,1] 
c_evn_all$NMDS2<-c_wm_nmds_all_metaMDS$points[ ,2] 
data.scores.all <- as.data.frame(scores(c_wm_nmds_all_metaMDS))  
data.scores.all$sample_ID <- rownames(data.scores.all) 
species.scores.all <- as.data.frame(scores(c_wm_nmds_all_metaMDS, "species")) 
species.scores.all$species <- rownames(species.scores.all)


# function for ellipses
#taken from the excellent stackoverflow Q+A: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the water mass factor
wm.df.all <- data.frame() #sets up a data frame before running the function.
for(g in levels(c_evn_all$Water_mass)){
  wm.df.all <- rbind(wm.df.all, cbind(as.data.frame(with(c_evn_all [c_evn_all$Water_mass==g,],
                                                 veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                              ,Water_mass=g))
}

# data for labelling the ellipse
NMDS.mean.cwm.all=aggregate(c_evn_all[ ,c("NMDS1", "NMDS2")], 
                        list(group = c_evn_all$Water_mass), mean)

# data for labelling the ellipse
NMDS.mean.all=aggregate(c_evn_all[,c("NMDS1", "NMDS2")], 
                    list(group = c_evn_all$Water_mass), mean)
#plot
nmds_all_wm <- ggplot()+
  geom_path(data = wm.df.all, aes(x = NMDS1, y = NMDS2, group = Water_mass, color = Water_mass))+ 
  geom_text(data=species.scores.all,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size = 5) + 
  geom_point(data=c_evn_all,aes(x=NMDS1,y=NMDS2,shape=Depth,colour=Water_mass),size=4) + 
  geom_text(data=c_evn_all,aes(x=NMDS1,y=NMDS2,label=StationName),size=4,vjust=-0.8, hjust=0.8)+
  scale_colour_manual(values=c("PSWw" = "blue", "AW" = "red", "EBDW" = "purple")) +
  #coord_equal() +
  theme_plot+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())+
  annotate("text", x =-0.45, y = -0.3, label=paste('Stress =',round(c_wm_nmds_all_metaMDS$stress,3)))
nmds_all_wm

#PCA
pca.all.data.prcomp <- prcomp(c_wm_nmds_all, center = TRUE,scale. = TRUE)
summary(pca.all.data.prcomp)
str(pca.all.data.prcomp)
ggbiplot(pca.all.data.prcomp, labels=rownames(c_wm_nmds_all))

#PCoA

#Evaluate the distances
EuclDistAT <- dist(c_wm_nmds_all, method = "euclidean")
PcoA2D_AT <- cmdscale(EuclDistAT, k =2, eig=TRUE)
PcoA2D_AT$eig
PcoA2D_AT$eig[1:2]/sum(PcoA2D_AT$eig)*100

pcoa.scores.all <- as.data.frame(scores(PcoA2D_AT))  
pcoa.species.scores.all <- as.data.frame(wascores(PcoA2D_AT$points[,1:2], c_wm_nmds_all))

names(pcoa.species.scores.all)[1:2] <- c('PC1', 'PC2')
names(pcoa.scores.all)[1:2] <- c('PC1', 'PC2')

env_pcoa <- env

PcoA2D_AT$Water_mass <- env_pcoa$Water_mass
PcoA2D_AT$Depth <- env_pcoa$Depth
PcoA2D_AT$StationName <- env_pcoa$StationName
pcoa.species.scores.all$species <- rownames(pcoa.species.scores.all)
pcoa.scores.all$sample_ID <- rownames(pcoa.scores.all)  
pcoa.scores.all$Water_mass <- env_pcoa$Water_mass
pcoa.scores.all$Depth <- env_pcoa$Depth
pcoa.scores.all$Depth <- factor(pcoa.scores.all$Depth, levels = c("DCM","EPI","MESO","BATHY","ABYSS"))
pcoa.scores.all$StationName <- env_pcoa$StationName

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

#plot
Tr_PcoA2 <- ggplot()+
  geom_text(data=pcoa.species.scores.all,aes(x=PC1,y=PC2,label=species), alpha=0.5,size = 3)+ 
  geom_point(data=pcoa.scores.all,aes(x=PC1,y=PC2,shape=Depth,colour=Water_mass),size=4) + 
  geom_text(data=pcoa.scores.all,aes(x=PC1,y=PC2,label=StationName),size=4,vjust=-0.8, hjust=0.8)+
  scale_colour_manual(values=c("PSWw" = "blue", "AW" = "red", "EBDW" = "purple")) +
  #coord_equal() +
  theme_plot+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())

Tr_PcoA2 + scale_x_discrete(name ="PCoA1 [62.3%]")+ scale_y_discrete(name ="PCoA2 [15.8%]")+stat_ellipse(type = "t")

save.image("card_fish_water_column_ps99.Rdata")
