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
# Env. parameter coorelation test 
##################################

#Normality
shapiro.test(env$Temperature) # p-value = 1.806e-05
shapiro.test(env$Salinity) # p-value = 3.823e-09
shapiro.test(env$Chla_fluor) # p-value = 0.0001371
shapiro.test(env$NO3.NO2) # p-value = 0.002956
shapiro.test(env$NO2) # p-value = 7.592e-05
shapiro.test(env$NO3) # p-value = 0.003924
shapiro.test(env$SiO3) # p-value = 0.001508
shapiro.test(env$PO4) # p-value = 0.02344
shapiro.test(env$NH4) # p-value = 2.357e-07

#Pearson’s product moment correlation

#omitNA measurements
cor.test <- subset(env, StationName!="HG1")
cor.test <-na.omit(cor.test[,c("Temperature", "Salinity", "Chla_fluor", "NO3.NO2", "NO2", "NO3", "SiO3", "PO4", "NH4")])

#Calculate the Pearson’s correlation coefficient in a matrix
cor(cor.test)

pairs.panels(cor.test)

#Correlogram
x <- cor(cor.test)
corrplot(x, type="upper", order="hclust", col=brewer.pal(n=8, name="PuOr"), method="pie", tl.srt=45, addCoef.col = "black")

#Matrix of the p-value of the correlation
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(cor.test)
p.cor.mat <- as.data.frame(p.mat)

#Correlation matrix
chart.Correlation(cor.test, histogram=TRUE, pch=19)

##################################
# RDA for surface waters
##################################
#Redundancy test for constrained variables

## RDA surface for FISH counts
counts.rel.ab.surf <- subset(counts.rel.ab[counts.rel.ab$Depth %in% c("DCM","EPI"),])
counts.rda.surf <- counts.rel.ab.surf[,c("Depth", "Domain", "StationName", "proportion")]
counts.rda.surf %>% spread(Domain, proportion) -> counts.sp.rda
counts.rda <- counts.sp.rda[,c("Depth", "StationName", "ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER")]
env.DCM <- subset(env[env$Depth %in% c("DCM", "EPI"),])
counts.rda.meta <- env.DCM %>% left_join(counts.rda, by=c("StationName", "Depth"))
rownames(counts.rda.meta) <- env.DCM$sample_ID
rownames(counts.rda) <- env.DCM$sample_ID # here check this 
rownames(env.DCM) <- env.DCM$sample_ID
counts.rda.meta$sample_ID <- as.factor(counts.rda.meta$sample_ID)
counts.rda.meta <- counts.rda.meta[colnames(counts.rda.meta) != "X"]
counts.rda.o <- counts.rda[colnames(counts.rda) != "StationName"]
counts.rda.ord <-counts.rda.o[colnames(counts.rda.o) != "Depth"]
counts.rda.meta$Region[counts.rda.meta$StationName %in% EGC] <- "EGC"
counts.rda.meta$Region[counts.rda.meta$StationName %in% WSC] <- "WSC"
counts.rda.meta[is.na(counts.rda.meta)] <- 0

#redundancy test
ord <- rda(counts.rda.ord ~ Temperature + Chla_fluor + NH4 + SiO3, data=counts.rda.meta, scale = TRUE, center = FALSE)

#generate dataframe for the samples
ps_scores <- vegan::scores(ord)
sites <- data.frame(ps_scores$sites)
sites$sample_ID <- as.factor(rownames(sites))

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
  geom_point(data = sites, aes(x = RDA1, y = RDA2, colour = Water_mass, shape = Depth),  size = 4) +
  geom_text(data = sites,aes(x = RDA1, y = RDA2,label = StationName), nudge_y= 0.2,size=4)+
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("\nAxis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]\n", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ord$CCA$eig[2] / ord$CCA$eig[1])) +
  theme_plot+
  geom_segment(mapping = arrow_mapFL, size = .5, data = arrowdfFL, color = "black", arrow = arrowhead) +
  geom_text(mapping = label_mapFL, size = 4, data = arrowdfFL, show.legend = FALSE, nudge_y= 0.2)+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))
rda_plot + scale_colour_manual(values=c("PSWw"="blue","AW"="red", "EBDW"="violet"))+
  theme_plot


save.image("card_fish_water_column_ps99.Rdata")


