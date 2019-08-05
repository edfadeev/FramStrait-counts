#####################################
## set working directory and load preprocessed dataset
#####################################
#set working directory in iOS
setwd("~/CARD-FISH/CARD-FISH_water_project/")

###################################
## Load required libraries
###################################
library(ggplot2)#; packageVersion("ggplot2")
library(dplyr)#; packageVersion("dplyr")
library(ggsignif)#; packageVersion("ggsignif")
library(cowplot)#; packageVersion("cowplot")
library("PerformanceAnalytics")#; packageVersion("PerformanceAnalytics")
library(tidyr)
library(vegan)
library(corrplot)
library(RColorBrewer)


###################################
## defined functions
###################################
#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

#scale parameters
scale_par <- function(x) scale(x, center = FALSE, scale = TRUE)[,1]

# fit text in facets:
swr = function(string, nwrap=15) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
###################################
## Import counts and calculate cell concentrations
###################################
#raw counts
raw.counts.SH <- read.csv("FOV_all_groups_SH.csv", sep = ",", dec = ".", header = TRUE)

#split sample name
raw.counts.SH <- raw.counts.SH %>% 
  separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-")

#remove SV2 station and ABYSS depth
raw.counts.SH <- subset(raw.counts.SH, !StationName == "SV2")
raw.counts.SH <- subset(raw.counts.SH, !Depth == "ABYSS")

###calculate cell concentration for each field of view
#microscope calculation factor
calc.factor <- 99515.5458411807
#DAPI
raw.counts.SH$DAPI_conc <- (raw.counts.SH$DAPI_Nr_Set*calc.factor)/raw.counts.SH$Volume
#FISH
raw.counts.SH$FISH_conc <- (raw.counts.SH$DAPI_Nr_SubSet*calc.factor)/raw.counts.SH$Volume

#calculate DAPI concetration per sample
raw.counts.SH %>% 
  group_by(StationName, Depth, Domain) %>% 
  summarise (DAPI.conc.mn = mean(DAPI_conc),
             DAPI.conc.md =  median(DAPI_conc), 
             DAPI.conc.sd = sd(DAPI_conc),
             FISH.conc.mn = mean(FISH_conc),
             FISH.conc.md =  median(FISH_conc), 
             FISH.conc.sd = sd(FISH_conc),
             n = length(DAPI_conc)) -> counts_all


#Define regions by stations
EGC <- c("EG1","EG4")
N <- c("N3","N4","N5")
WSC <- c("HG9","HG7","HG5","HG4","HG2","HG1")

counts_all$Region <- "EGC"
counts_all$Region[counts_all$StationName %in% N] <- "N"
counts_all$Region[counts_all$StationName %in% WSC] <- "WSC"

###################################
##Check the coverage of EUB and ARCH probes in comparison to DAPI
###################################
#Coverage of EUB and ARCH per depth and station
counts_all%>% 
  filter(Domain %in% c("EUB","ARCH"))%>%
  group_by(Region, StationName, Depth, Domain) %>% 
  summarise(proportion = FISH.conc.mn/DAPI.conc.mn)-> EUB_ARCH_prop

#Total coverage (of EUB and ARCH sum) per depth and station
EUB_ARCH_prop%>%
  group_by(Region, StationName, Depth) %>% 
  summarise(total.coverage = sum(proportion)) -> EUB_ARCH_coverage
#!!!need to figure out how to calculate the SD here !#

#Coverage by regions and depth
EUB_ARCH_coverage %>%
  group_by(Region, Depth) %>% 
  summarise(mean.coverage = mean(total.coverage),
            sd.coverage = se(total.coverage)) -> EUB_ARCH_coverage_by_region

#Significance coverage test between the regions
wilcox.test(EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "WSC",]$total.coverage,
            EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "N",]$total.coverage)
# Wilcoxon rank sum test
# 
# data:  EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "WSC", ]$total.coverage and EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "N", ]$total.coverage
# W = 119, p-value = 0.4164
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "WSC",]$total.coverage,
            EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "EGC",]$total.coverage)

# Wilcoxon rank sum test
# 
# data:  EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "WSC", ]$total.coverage and EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "EGC", ]$total.coverage
# W = 49, p-value = 0.1041
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "N",]$total.coverage,
            EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "EGC",]$total.coverage)

# Wilcoxon rank sum test
# 
# data:  EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "N", ]$total.coverage and EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "EGC", ]$total.coverage
# W = 24, p-value = 0.1422

#Coverage by depth only, takes the mean of the regions
EUB_ARCH_coverage %>%
  group_by(Depth) %>% 
  summarise(mean.coverage = mean(total.coverage),
            sd.coverage = se(total.coverage)) -> EUB_ARCH_coverage_by_depth


###################################
##significance of DAPI  with depth
###################################
counts_all%>% 
  group_by(Region, StationName, Depth) %>% 
  summarise(mean.DAPI = mean(DAPI.conc.mn),
            sd.dapi =se(DAPI.conc.mn))-> DAPI_sum

wilcox.test(DAPI_sum[DAPI_sum$Depth == "DCM",]$mean.DAPI,
            DAPI_sum[DAPI_sum$Depth == "EPI",]$mean.DAPI)

# Wilcoxon rank sum test
# 
# data:  DAPI_sum[DAPI_sum$Depth == "DCM", ]$mean.DAPI and DAPI_sum[DAPI_sum$Depth == "EPI", ]$mean.DAPI
# W = 93, p-value = 0.0336
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(DAPI_sum[DAPI_sum$Depth == "EPI",]$mean.DAPI,
            DAPI_sum[DAPI_sum$Depth == "MESO",]$mean.DAPI)

# Wilcoxon rank sum test
# 
# data:  DAPI_sum[DAPI_sum$Depth == "EPI", ]$mean.DAPI and DAPI_sum[DAPI_sum$Depth == "MESO", ]$mean.DAPI
# W = 121, p-value = 2.835e-06
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(DAPI_sum[DAPI_sum$Depth == "MESO",]$mean.DAPI,
            DAPI_sum[DAPI_sum$Depth == "BATHY",]$mean.DAPI)

# Wilcoxon rank sum test
# 
# data:  DAPI_sum[DAPI_sum$Depth == "MESO", ]$mean.DAPI and DAPI_sum[DAPI_sum$Depth == "BATHY", ]$mean.DAPI
# W = 90, p-value = 0.01272
# alternative hypothesis: true location shift is not equal to 0

DAPI_sum %>% 
  group_by(Region, Depth) %>% 
  summarise(mean = mean(mean.DAPI),
            sd.dapi =se(mean.DAPI))-> DAPI_sum_region
###################################
##significance of Bacteria with depth
###################################
counts_all%>% 
  filter(Domain %in% c("EUB"))%>%
  group_by(Region, StationName, Depth) %>% 
  summarise(mean.FISH = mean(FISH.conc.mn),
            sd.FISH =se(FISH.conc.mn))-> EUB_sum

wilcox.test(EUB_sum[EUB_sum$Depth == "DCM",]$mean.FISH,
            EUB_sum[EUB_sum$Depth == "EPI",]$mean.FISH)

# Wilcoxon rank sum test
# 
# data:  EUB_sum[EUB_sum$Depth == "DCM", ]$mean.FISH and EUB_sum[EUB_sum$Depth == "EPI", ]$mean.FISH
# W = 105, p-value = 0.002447
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(EUB_sum[EUB_sum$Depth == "EPI",]$mean.FISH,
            EUB_sum[EUB_sum$Depth == "MESO",]$mean.FISH)

# Wilcoxon rank sum test
# 
# data:  EUB_sum[EUB_sum$Depth == "EPI", ]$mean.FISH and EUB_sum[EUB_sum$Depth == "MESO", ]$mean.FISH
# W = 121, p-value = 2.835e-06
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(EUB_sum[EUB_sum$Depth == "MESO",]$mean.FISH,
            EUB_sum[EUB_sum$Depth == "BATHY",]$mean.FISH)

# Wilcoxon rank sum test
# 
# data:  EUB_sum[EUB_sum$Depth == "MESO", ]$mean.FISH and EUB_sum[EUB_sum$Depth == "BATHY", ]$mean.FISH
# W = 90, p-value = 0.01272
# alternative hypothesis: true location shift is not equal to 0

###################################
##significance of Groups between the regions 
###################################
counts_all%>% 
  filter(Domain %in% c("POL"))%>% # here change domain to see differences
  group_by(Region, StationName, Depth) %>% 
  summarise(mean.FISH = mean(FISH.conc.mn),
            sd.FISH =se(FISH.conc.mn))-> GAM_sum

wilcox.test(GAM_sum[GAM_sum$Region == "EGC",]$mean.FISH,
            GAM_sum[GAM_sum$Region == "N",]$mean.FISH)

# Wilcoxon rank sum test
# 
# data:  SAR11_sum[SAR11_sum$Region == "EGC", ]$mean.FISH and SAR11_sum[SAR11_sum$Region == "N", ]$mean.FISH
# W = 37, p-value = 0.7108
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(GAM_sum[GAM_sum$Region == "EGC",]$mean.FISH,
            GAM_sum[GAM_sum$Region == "WSC",]$mean.FISH)

# Wilcoxon rank sum test
# 
# data:  SAR11_sum[SAR11_sum$Region == "EGC", ]$mean.FISH and SAR11_sum[SAR11_sum$Region == "WSC", ]$mean.FISH
# W = 75, p-value = 0.6945
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(GAM_sum[GAM_sum$Region == "N",]$mean.FISH,
            GAM_sum[GAM_sum$Region == "WSC",]$mean.FISH)

# Wilcoxon rank sum test
# 
# data:  SAR11_sum[SAR11_sum$Region == "N", ]$mean.FISH and SAR11_sum[SAR11_sum$Region == "WSC", ]$mean.FISH
# W = 154, p-value = 0.7533
# alternative hypothesis: true location shift is not equal to 0

###################################
##significance of deep DAPI between regions
###################################
counts_all%>% 
  filter(Depth %in% c("BATHY"))%>%
  group_by(Region, StationName, Depth) %>% 
  summarise(mean.DAPI = mean(DAPI.conc.mn),
            sd.dapi =se(DAPI.conc.mn))-> DAPI_sum_meso

wilcox.test(DAPI_sum_meso[DAPI_sum_meso$Region == "N",]$mean.DAPI,
            DAPI_sum_meso[DAPI_sum_meso$Region == "WSC",]$mean.DAPI)

##################################
##Plot EUB and ARCH profiles by region
###################################
counts_all%>% 
  filter(Domain %in% c("EUB","ARCH"))%>%
  group_by(Region, StationName, Depth, Domain) -> EUB_ARCH_absolute
EUB_ARCH_absolute$Depth <- factor(EUB_ARCH_absolute$Depth, levels = c("DCM","EPI","MESO","BATHY"))
counts_all$Depth <- factor(EUB_ARCH_absolute$Depth, levels = c("DCM","EPI","MESO","BATHY"))
depths.labs <- c("Surface", "Epipelagic", "Mesopelagic", "Bathypelagic")
names(depths.labs) <- c("DCM", "EPI", "MESO", "BATHY")
domain.labs <- c("Archaea", "Bacteria")
names(domain.labs) <- c("ARCH", "EUB")
EUB_ARCH_absolute$FISH.conc.mn <- log10(EUB_ARCH_absolute$FISH.conc.mn)
EUB_ARCH_vertical_boxplot <- ggplot(EUB_ARCH_absolute, aes(x= Region, y = FISH.conc.mn))+
  geom_boxplot(aes(fill = Region))+
  facet_grid(Domain~Depth, labeller = labeller(Depth = depths.labs, Domain = domain.labs))+
  geom_jitter(aes(x= Region, y = FISH.conc.mn),width = 0.2, alpha = 0.3)+
  #geom_signif(comparisons = list(c("EGC", "WSC"),c("EGC","N"),c("N","WSC")), 
  #            map_signif_level=TRUE, test = "wilcox.test")+
  #scale_y_continuous(trans = 'log10')+
  labs(y = "Cell abundance [log10(Cells/mL)]")+#coord_cartesian(ylim = c(1e+03, 1e+06))+
  labs(x = "")+
  scale_fill_manual(values=c("WSC" = "red", "EGC" = "blue", "N" = "gray")) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "bottom")
EUB_ARCH_vertical_boxplot

DAPI_absolute <- counts_all[, c("DAPI.conc.mn","Region", "Depth")]
DAPI_absolute$DAPI.conc.mn <- log10(DAPI_absolute$DAPI.conc.mn)
DAPI_vertical_boxplot <- ggplot(DAPI_absolute, aes(Region, DAPI.conc.mn))+
  scale_fill_manual(values=c("WSC" = "red", "EGC" = "blue", "N" = "gray"))+
  geom_boxplot(aes(fill = Region))+facet_grid(~Depth, labeller = labeller(Depth = depths.labs))+
  geom_jitter(aes(x= Region, y = DAPI.conc.mn),width = 0.2, alpha = 0.3)+
  labs(y = "Cell abundance [log10(Cells/mL)]")+coord_cartesian(ylim = c(3, 7))+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "bottom")
DAPI_vertical_boxplot 

Figure2 <- plot_grid(EUB_ARCH_vertical_boxplot, DAPI_vertical_boxplot, labels = c("A","B"), nrow = 2)


#save the figure
ggsave("./Figure-2AB.png", 
       plot = Figure2,
       scale = 1,
       units = "cm",
       width = 17.8,
       height = 25,
       dpi = 300)

###################################
# Overview of phylogenetic group abundances by regions
###################################
#calculate abundance in surface by region
counts_all%>% 
  filter(Depth == "DCM")%>%
  group_by(Region, Domain) %>%
  summarise(mean.abund = mean(FISH.conc.mn),
            se.abund = se(FISH.conc.mn),
            n= length(FISH.conc.mn))-> FISH_abundance_surface_by_region

#make wide table
FISH_abundance_surface_by_region %>%
  select(Region,Domain,mean.abund)%>%
  spread(Region,mean.abund) -> FISH_abundance_surface_by_region_wide

#calculate abundances by region
counts_all%>% 
  group_by(Region, Domain, Depth) %>%
  summarise(mean.abund = mean(FISH.conc.mn),
            se.abund = se(FISH.conc.mn),
            n= length(FISH.conc.mn))-> FISH_abundances_by_region
FISH_abundances_by_region$Depth <- factor(FISH_abundances_by_region$Depth, levels = c("DCM","EPI","MESO","BATHY"))
FISH_abundances_by_region %>%
  select(Region,Domain,Depth,mean.abund)%>%
  spread(Domain,mean.abund) -> FISH_abundance_by_region_wide
#write.csv(FISH_abundance_by_region_wide, file="~/CARD-FISH/CARD-FISH_water_project/cell_densities.csv")

#calculate proportional abundance of the different groups at the surface
counts_all%>% 
  filter(Depth %in% c("DCM"))%>%
  #left_join(.,total_abundance_surface_by_station[,c("StationName","total.abund")], by = "StationName") %>%
  group_by(Region, Depth, StationName, Domain)%>%
  mutate(proportion = (FISH.conc.mn / DAPI.conc.mn)*100) %>%
  select(Region, StationName, Domain, FISH.conc.mn, proportion)  -> surface_FISH_proportion

#summarize by regions
surface_FISH_proportion%>%
  group_by(Region, Domain) %>% 
  summarise (mean.abund = mean(FISH.conc.mn),
             se.abund = se(FISH.conc.mn), 
             mean.prop = mean(proportion),
             se.prop = se(proportion)) -> surface_FISH_proportion_by_regions

#calculate proportional abundance of the different groups
counts_all%>% 
  group_by(Region, StationName, Domain)%>%
  mutate(proportion = (FISH.conc.mn / DAPI.conc.mn)*100) %>%
  group_by(Region, Domain, Depth) %>%
  summarise(mean.prop = mean(proportion),
            se.prop = se(proportion))%>%
  select(Region, Domain, mean.prop, Depth) %>%
  spread(Domain,mean.prop) ->  FISH_proportion

##################################
# NMDS plot
##################################
#list all taxa (excluding EUB and ARCH)
taxa.nmds <- c("ALT", "BACT", "CFX", "CREN", "DELTA", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202","SAR324", "SAR406", "VER")

#generate wide abundance table 
surface_FISH_proportion %>% 
  select(Region, StationName, Domain, FISH.conc.mn) %>%
  group_by(Region, StationName) %>%
  filter(Domain %in% taxa.nmds) %>% 
  spread(Domain, FISH.conc.mn) -> surface_FISH_wide

#Calculate distances and generate NMDS
all_metaMDS <- metaMDS(surface_FISH_wide[,taxa.nmds], maxit= 100, trace=TRUE)

data.scores <- as.data.frame(scores(all_metaMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- surface_FISH_wide$StationName  # create a column of site names, from the rownames of data.scores
data.scores$Region <- surface_FISH_wide$Region

species.scores <- as.data.frame(scores(all_metaMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

#Plot NMDS
NMDS_surface <- ggplot() + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2, colour = Region),size=5) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=3,nudge_y =-0.03) +  # add the site labels
  scale_colour_manual(values=c("WSC" = "red", "N" = "black", "EGC" = "blue")) +
  annotate(geom="text", x=0.25, y=-0.2, label= paste("Stress =", round(all_metaMDS$stress, 3), sep = " "),color="black")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "bottom")
NMDS_surface

#save the figure
ggsave("./Figure-NMDS.png", 
       plot = NMDS_surface,
       scale = 1,
       units = "cm",
       width = 19.5,
       height = 17.4,
       dpi = 300)

##################################
# Env. parameter coorelation test 
##################################
## import environmental metadata 
metadata <- read.csv("PS99_samples_meta_EF_MC_nutrient_corrected.csv", sep = ",", dec = ".", header = TRUE)

#remove station SV2 from the dataset
metadata <- subset(metadata, !StationName == "SV2")

#list environmental parameters that need to be scaled, NH4 negative values were transformed to 0, NO2 not included
env.par <- c("Temperature","Salinity", 
             "Chla_fluor", "d.NO3.NO2.", "d.NO3.", "d.NH4.", "d.PO4.", "d.SiO3.")

#drop rows with NA and scale the env. parameters 
metadata %>% 
  drop_na() %>% 
  mutate_at(env.par, scale_par) -> env

#Calculate the Pearson’s correlation coefficient in a matrix
chart.Correlation(env[,env.par], histogram=TRUE, pch=19)
cor.env <- cor(env[,env.par])
corrplot(cor.env, type="upper", order="hclust", col=brewer.pal(n=8, name="PuOr"), method="pie", addCoef.col = "black",srt = 45, number.cex= 0.75, tl.col = "white")
text(1:8, 9.5, expression("Temperature [°C]","Chlorophyll a [µg/L]", "Salinity [PSU]", "ΔSiO"[3], "ΔNH"[4],"ΔPO"[4], "ΔTotal N", "ΔNO"[3]), srt = 45)
text(-0.5, 1:8,  expression("                                                                                                                                                                      ΔNO"[3],
                            "                                                                                                                                           ΔTotal N",
                            "                                                                                                                  ΔPO"[4],
                            "                                                                                               ΔNH"[4], 
                            "                                                                          ΔSiO"[3], 
                            "                                 Salinity [PSU]",
                            "                          Chlorophyll a [µg/L]", "Temperature [°C]" ))   

##################################
# correlation between env. par. and absolute counts
##################################

taxa.sp <- c("CREN","ARCH", "DELTA", "POL","SAR202","SAR324","SAR406","EUB","BACT","CFX","GAM","OPI","ROS","SAR11", "VER","ALT")
#select only surface samples
env %>%
  filter (Depth == "DCM") -> env.SRF

#env.SRF <- env.SRF[,c("StationName", "Depth", "Temperature", "Salinity", "Chla_fluor", "d.NO3.NO2.", "d.NO3.", "d.NO2.","d.SiO3.", "d.PO4.", "d.NH4.")]

#generate wide abundance table 
counts_all%>% 
  filter(Depth == "DCM")%>%
  group_by(StationName, Domain) %>%
  summarise(mean.abund = mean(FISH.conc.mn),
            se.abund = se(FISH.conc.mn),
            n= length(FISH.conc.mn)) -> counts_for_cor

counts_for_cor$Domain <- factor(counts_for_cor$Domain, levels = c("CREN","ARCH", "DELTA", "POL","SAR202","SAR324","SAR406","EUB","BACT","CFX","GAM","OPI","ROS","SAR11", "VER","ALT"))


select(counts_for_cor, c(Domain, StationName, mean.abund))%>%
  spread(Domain, mean.abund) -> surface_FISH_abundances_wide  

#drop samples with no env. data
surface_FISH_abundances_wide %>%
  filter(StationName %in% as.vector(env.SRF$StationName)) -> surface_FISH_counts

#calculate pearson correlation between each taxa and env. parameters.
cor.env.counts <- cor(env.SRF[,env.par], surface_FISH_counts[,taxa.sp], method = "spearman")
table.spearman <- t(as.data.frame(cor.env.counts))

write.csv(table.spearman, file="~/CARD-FISH/CARD-FISH_water_project/spearman_table.abs.csv")


#multiple correlation tests with significance p-values
#install the packages if necessary
if(!require("tidyverse")) install.packages("tidyverse")
if(!require("broom")) install.packages("broom")
if(!require("fs")) install.packages("fs")
if(!require("lubridate")) install.packages("lubridate")

#load packages
library(tidyverse)
library(broom)
library(fs)
library(lubridate)


#merge counts and environmental data
data_all <- left_join(surface_FISH_abundances_wide, env.SRF[,c("StationName",env.par)] , by = "StationName")

#generate long table
data <- gather(data_all, Domain, Abund, taxa.sp)%>%
  gather(variable, value, env.par)

#nest the table according to taxa and env. variable
data_nest <- group_by(data, Domain, variable) %>% nest()
data_nest

#define function for correlation 
cor_fun <- function(df) cor.test(df$Abund, df$value, method = "spearman") %>% tidy()

#nested correlations tests
data_nest <- mutate(data_nest, model = map(data, cor_fun))
data_nest

#summary table
corr_pr <- select(data_nest, -data) %>% unnest()
corr_pr

##################################
# Depth profiles by taxa with log10 of absolute counts
##################################

swr = Vectorize(swr)
domain.depth.labs = swr(domain.depth.labs)

FISH_abundances_by_region$Depth<- factor(FISH_abundances_by_region$Depth, 
                                         levels = rev(c("DCM", "EPI", "MESO", "BATHY")))
FISH_abundances_by_region$log <- log10(FISH_abundances_by_region$mean.abund)
FISH_abundances_by_region <- subset(FISH_abundances_by_region, Domain != "EUB") 
FISH_abundances_by_region <- subset(FISH_abundances_by_region, Domain != "ARCH") 
FISH_abundances_by_region$Domain <- factor(FISH_abundances_by_region$Domain, levels =c("GAM", "ALT", "BACT", "POL", "ROS", "SAR11", "DELTA", "SAR324", "CFX", "SAR202", "VER", "OPI", "SAR406", "CREN"))
domain.depth.labs <- c( "Gammaproteobacteria", "Alteromonadaceae Colwelliaceae Pseudoalteromonadaceae", "Bacteroidia", "Polaribacter", "Rhodobacteraceae", "SAR11", "Deltaproteobacteria", "SAR324", "Chloroflexi", "SAR202", "Verrucomicrobia", "Opitutale", "SAR406", "Thaumarcheota")
names(domain.depth.labs) <- c("GAM", "ALT", "BACT", "POL", "ROS", "SAR11", "DELTA", "SAR324", "CFX", "SAR202", "VER", "OPI", "SAR406", "CREN")


Depth_profiles_abs_abundance.p <- ggplot(FISH_abundances_by_region, aes(y = log, x = Depth))+
  geom_line(aes(color =Region, group=Region),linetype=1)+
  coord_flip()+
  facet_wrap(.~Domain, labeller = labeller(Domain = domain.depth.labs))+
  labs(y = "Cell abundance [log10(Cells/mL)]", x = "Depth zone")+
  scale_x_discrete(labels=c("DCM" = "Surface", "EPI" = "Epipelagic","MESO" = "Mesopelagic", "BATHY" = "Bathypelagic", "ABYSS" = "Abyssopelagic"))+
  scale_color_manual(breaks = c("EGC", "N", "WSC"), values = c("blue", "grey", "red"))+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(face = "italic"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

Depth_profiles_abs_abundance.p

#save the figure
ggsave("./Figure-3.png", 
       plot = Depth_profiles_abs_abundance.p,
       scale = 1,
       units = "cm",
       width = 25,
       height = 25,
       dpi = 300)
