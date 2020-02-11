###################################
## Load required libraries
###################################
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(ggsignif); packageVersion("ggsignif")
library(cowplot); packageVersion("cowplot")
library(tidyr); packageVersion("tidyr")
library(PerformanceAnalytics); packageVersion("PerformanceAnalytics")

library(tidyverse)
library(broom)
library(fs)
library(lubridate)

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

###################################
## Import counts and calculate cell concentrations
###################################
#raw counts
raw.counts.SH <- read.csv("./Rwork/FOV_all_groups_SH.csv", sep = ",", dec = ".", header = TRUE)

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
#coverage by probes
counts_all%>% 
  filter(Domain %in% c("EUB","ARCH"))%>%
  group_by(Region, StationName, Depth, Domain) %>% 
  summarise(proportion = FISH.conc.mn/DAPI.conc.mn)-> EUB_ARCH_prop

#total coverage
EUB_ARCH_prop%>%
  group_by(Region, StationName, Depth) %>% 
  summarise(total.coverage = sum(proportion)) -> EUB_ARCH_coverage
#!!!need to figure out how to calculate the SD here !#

#proportions by regions
EUB_ARCH_coverage %>%
  group_by(Region, Depth) %>% 
  summarise(mean.coverage = mean(total.coverage),
            sd.coverage = se(total.coverage)) -> EUB_ARCH_coverage_by_region

#siginificance tests
wilcox.test(EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "WSC",]$total.coverage,
            EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "N",]$total.coverage)

wilcox.test(EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "WSC",]$total.coverage,
            EUB_ARCH_coverage[EUB_ARCH_coverage$Region == "EGC",]$total.coverage)


#proportions by depth
EUB_ARCH_coverage %>%
  group_by(Depth) %>% 
  summarise(mean.coverage = mean(total.coverage),
            sd.coverage = se(total.coverage)) -> EUB_ARCH_coverage_by_depth

###################################
##Plot EUB and ARCH profiles by region
###################################
counts_all$Depth <- factor(counts_all$Depth, levels = c("DCM","EPI","MESO","BATHY"))

DAPI_vertical_boxplot <- ggplot(counts_all, aes(x= Region, y = DAPI.conc.mn))+
  geom_boxplot(aes(fill = Region))+
  facet_grid(Depth~.)+
  geom_jitter(aes(x= Region, y = DAPI.conc.mn),width = 0.2, alpha = 0.3)+
  #geom_signif(comparisons = list(c("EGC", "WSC"),c("EGC","N"),c("N","WSC")), 
  #            map_signif_level=TRUE, test = "wilcox.test")+
  scale_y_log10(name = "cell concentration (Cells/mL)")+
  scale_fill_manual(values=c("WSC" = "red", "EGC" = "blue", "N" = "gray")) +
  theme_bw()+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")

counts_all%>% 
  filter(Domain %in% c("EUB","ARCH"))%>%
  group_by(Region, StationName, Depth, Domain) -> EUB_ARCH_absolute

EUB_ARCH_absolute$Depth <- factor(EUB_ARCH_absolute$Depth, levels = c("DCM","EPI","MESO","BATHY"))

p <- list()
for (i in c("EUB","ARCH")){
p[[i]] <-  ggplot(EUB_ARCH_absolute[EUB_ARCH_absolute$Domain== i,], aes(x= Region, y = FISH.conc.mn))+
  geom_boxplot(aes(fill = Region))+
  facet_grid(Depth~.)+
  geom_jitter(aes(x= Region, y = FISH.conc.mn),width = 0.2, alpha = 0.3)+
  #geom_signif(comparisons = list(c("EGC", "WSC"),c("EGC","N"),c("N","WSC")), 
  #            map_signif_level=TRUE, test = "wilcox.test")+
  scale_y_log10(name = "cell concentration (cells/mL)")+
  scale_fill_manual(values=c("WSC" = "red", "EGC" = "blue", "N" = "gray")) +
  theme_bw()+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")
}

plot_grid(DAPI_vertical_boxplot, p[["EUB"]],p[["ARCH"]], ncol =3)

#save the figure
ggsave("./figures/Figure-2.pdf", 
       plot = last_plot(),
       scale = 1,
       units = "cm",
       #width = 17.8,
       #height = 17.4,
       dpi = 300)


#calculate mean cell densities in surface by station
counts_all%>% 
  filter(Domain %in% c("EUB","ARCH") & Depth == "DCM")%>%
  group_by(Region, StationName) %>%
  summarise(total.abund = sum(FISH.conc.mn))-> total_abundance_surface_by_station

total_abundance_surface_by_station %>%
  group_by(Region) %>%
  summarise(mean.abund = mean(total.abund),
            se.abund = se(total.abund))-> mean_abundance_surface_by_region

###################################
# Overview of taxonomic group abundances by regions
###################################
#calculate abundnace in surface by region
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

#calculate proportional abundance of the different groups
counts_all%>% 
  filter(Depth %in% c("DCM"))%>%
  left_join(.,total_abundance_surface_by_station[,c("StationName","total.abund")], by = "StationName") %>%
  group_by(Region, StationName, Domain)%>%
  mutate(proportion = (FISH.conc.mn / total.abund)) %>%
  select(Region, StationName, Domain, FISH.conc.mn, proportion)  -> surface_FISH_proportion


#summarize by regions
surface_FISH_proportion%>%
  group_by(Region, Domain) %>% 
  summarise (mean.abund = mean(FISH.conc.mn),
             se.abund = se(FISH.conc.mn), 
             mean.prop = mean(proportion),
             se.prop = se(proportion)) -> surface_FISH_by_regions


##################################
# NMDS plot
##################################
#list all taxa (excluding EUB and ARCH)
taxa <- c("ALT", "BACT", "CFX", "CREN", "DELTA", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202","SAR324", "SAR406", "VER")

#generate wide abundance table 
surface_FISH_proportion %>% 
  select(Region, StationName, Domain, FISH.conc.mn) %>%
  group_by(Region, StationName) %>%
  filter(Domain %in% taxa) %>% 
  spread(Domain, FISH.conc.mn) -> surface_FISH_wide

#Calculate distances and generate NMDS
all_metaMDS <- metaMDS(surface_FISH_wide[,taxa], maxit= 100, trace=TRUE)

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
  annotate(geom="text", x=-0.20, y=0.25, label= paste("Stress =", round(all_metaMDS$stress, 3), sep = " "),color="black")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "bottom")

#save the figure
ggsave("./figures/Figure-NMDS.png", 
       plot = NMDS_surface,
       scale = 1,
       units = "cm",
       #width = 17.8,
       #height = 17.4,
       dpi = 300)


##################################
# Env. parameter coorelation test 
##################################
## import environmental metadata 
metadata <- read.csv("./Rwork/PS99_samples_meta_EF_MC_nutrient_corrected.csv", sep = ",", dec = ".", header = TRUE)

#remove station SV2 from the dataset
metadata <- subset(metadata, !StationName == "SV2")

#list environmental parameters that need to be scaled
env.par <- c("Temperature","Salinity", 
             "Chla_fluor", "d.NO3.NO2.", "d.NO3.", "d.NH4.", "d.PO4.", "d.SiO3.")

#drop rows with NA and scale the env. parameters 
metadata %>% 
  drop_na() %>% 
  mutate_at(env.par, scale_par) -> env

#Calculate the Pearson’s correlation coefficient in a matrix
chart.Correlation(env[,env.par], histogram=TRUE, pch=19)

##################################
# correlation between env. par. and counts
##################################
#list all taxa (excluding EUB and ARCH)
taxa <- c("ALT", "BACT", "CFX", "CREN", "DELTA", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202","SAR324", "SAR406", "VER")

#generate wide abundance table 
counts_all%>% 
  filter(Depth %in% c("DCM"))%>%
  group_by(Region, StationName, Domain)%>%
  select(Region, StationName, Domain, FISH.conc.mn)  %>% 
  spread(Domain, FISH.conc.mn) -> surface_FISH_abundance_wide

#select only surface samples
env %>%
  filter (Depth == "DCM") -> env.SRF

#drop samples with no env. data
surface_FISH_abundance_wide %>%
  filter(StationName %in% env.SRF$StationName) -> surface_FISH_abundance_wide

#merge counts and environmental data
data_all <- left_join(surface_FISH_abundance_wide, env.SRF[,c("StationName",env.par)] , by = "StationName")

#generate long table
data <- gather(data_all, Domain, Abund, taxa)%>%
  gather(variable, value, env.par)

#nest the table according to taxa and env. variable
data_nest <- group_by(data, Domain, variable) %>% nest()
data_nest

#define function for correlation 
cor_fun <- function(df) cor.test(df$Abund, df$value, method = "pearson",conf.level = 0.95) %>% tidy()


#nested correlations tests
data_nest <- mutate(data_nest, model = map(data, cor_fun))
data_nest

#summary table
corr_pr <- select(data_nest, -data) %>% unnest()
corr_pr

corr_pr.sign <- filter(corr_pr, p.value < 0.05)


##################################
# Depth profiles by taxa 
##################################
counts_all%>%ungroup() %>% 
  group_by(Region, Depth, Domain) %>% 
  summarise (mean.abund = mean(FISH.conc.mn),
             se.abund = se(FISH.conc.mn))-> counts_FISH_by_regions


counts_FISH_by_regions %>% ungroup() %>% mutate_if(is.factor, as.character) %>% 
                          mutate(Depth=replace(Depth, Depth=="DCM", 20),
                                 Depth=replace(Depth, Depth=="EPI", 100),
                                  Depth=replace(Depth, Depth=="MESO", 1000),
                                  Depth=replace(Depth, Depth=="BATHY", 2000)) %>% 
                                mutate(Depth = as.numeric(Depth))  -> counts_FISH_by_regions

counts_FISH_by_regions$mean.abund <- as.numeric(counts_FISH_by_regions$mean.abund)


## Get points
ps <- data.frame(xspline(counts_FISH_by_regions[,c()], shape=-0.2, lwd=2, draw=F))


Depth_profiles_abs_abundance.p <- ggplot(counts_FISH_by_regions, aes(y = mean.abund, x = Depth))+
  geom_point(aes(color =Region, group=Region),shape=1)+
  geom_smooth(aes(colour = Region), method = "auto")+
  scale_y_log10()+
  scale_x_reverse()+
  coord_flip()+
  facet_wrap(.~Domain)+
  labs(y = "Cell abundance [log10(Cells/mL)]", x = "Depth zone")+
  #scale_x_discrete(labels=c("DCM" = "Surface", "EPI" = "Epipelagic","MESO" = "Mesopelagic", "BATHY" = "Bathypelagic", "ABYSS" = "Abyssopelagic"))+
  scale_color_manual(breaks = c("EGC", "N", "WSC"), values = c("blue", "gray", "red"))+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(face = "italic"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())






















##################################
# RDA for surface waters (not sure that this is the correct method to use, becuase it requires compositional data)
##################################





#change character to factors and merge with metadata
FISH.ra.SRF %>% spread(Domain,FISH.ra) %>% mutate_at(vars(StationName,Depth), funs(factor)) %>%
  dplyr::left_join(env, by=c("StationName", "Depth")) -> test

#list all taxa (excluding EUB)
taxa <- c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202","SAR324", "SAR406", "VER")
#generate an RDA plot 
ord <- rda(na.omit(test)[,taxa] ~ Temperature + Salinity + Chla_fluor + NH4 + SiO3 + PO4 + NO3.NO2, data=na.omit(test), scale = TRUE, center = FALSE)

#extract values for ggplot
rda.scores <- vegan::scores(ord,display=c("sp","wa","lc","bp","cn"))
rda.sites <- data.frame(rda.scores$sites)
rda.sites$StationName <- as.character(na.omit(test)$StationName)
rda.sites$Region <- as.character(na.omit(test)$Region)
rda.sites$Depth <- as.character(na.omit(test)$Depth)

#taxa
rda.species <- data.frame(rda.scores$species)
rda.species$species <- rownames(rda.species)


#Draw biplots
rda.arrows<- rda.scores$biplot
colnames(rda.arrows)<-c("x","y")
rda.arrows <- as.data.frame(rda.arrows)
rda.evals <- 100*(ord$CCA$eig / sum(ord$CCA$eig))

#Plot RDA
rda.plot <- ggplot() +
  geom_point(data = rda.sites, aes(x = RDA1, y = RDA2, color = Region, shape = Depth), 
             size = 4) +
  scale_color_manual(values=c("WSC" = "red", "EGC" = "blue"))+
  geom_text(data = rda.sites,aes(x = RDA1, y = RDA2,label = StationName), 
            nudge_y= -0.2,size=3)+
  geom_text(data=rda.species,aes(x=RDA1,y=RDA2,label=species),alpha=0.5) +  # add the species labels
  labs(x = sprintf("RDA1 [%s%%]", round(rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(rda.evals[2], 2))) +
  geom_segment(data=rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(rda.arrows*1.2),
            aes(x, y, label = rownames(rda.arrows)),color="black",alpha=0.5)+
  theme_bw()



#################################
## Plot cell density box plot
###################################
#order the depth layers
counts_all$Depth <- factor(counts_all$Depth, levels = c("DCM","EPI","MESO","BATHY"))

#plot DAPI in the different regions along the water column with siginificance test
boxplot_DAPI <- ggplot(counts_all, aes(x= Region, y = DAPI.conc.mn, label = StationName))+
  geom_boxplot(aes(fill = Region))+
  geom_jitter(width = 0.2, alpha = 0.3)+
  facet_grid(~Depth)+
  geom_signif(comparisons = list(c("EGC", "WSC"),c("EGC","N"),c("N","WSC")), 
              map_signif_level=TRUE, test = "wilcox.test")+
  scale_y_log10(name = "cell conc. [log10(Cells/mL)]")+
  scale_fill_manual(values=c("WSC" = "red", "EGC" = "blue", "N" = "gray")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "bottom")


#calculate mean cell densities 
counts_DAPI_overview <-as.data.frame(as.list(aggregate(DAPI.conc.mn~Region+Depth, 
                                                         counts_all, 
                                                         FUN = function(x) c(mean = mean(x), se = se(x), count=length(x)))))




#calculate mean cell densities by depth
counts_all %>% 
  filter(Domain %in% c("EUB","ARCH"))%>%
  group_by(Depth, Domain) %>% 
  summarise (cell.density.mean = mean(FISH.conc.mn),
             cell.density.se = se(FISH.conc.mn), 
             num.samples = length(FISH.conc.mn)) -> EUB_ARCH_absolute_by_depth

#calculate mean cell densities by region
counts_all %>% 
  filter(Domain %in% c("EUB","ARCH"))%>%
  group_by(Region, Depth, Domain) %>% 
  summarise (cell.density.mean = mean(FISH.conc.mn),
             cell.density.se = se(FISH.conc.mn), 
             num.samples = length(FISH.conc.mn)) -> EUB_ARCH_absolute_by_region



###################################
#Calculate and plot total cell abundance 
###################################





#calculate proportions of EUB and ARCH based on their total
#calculate the proportions
counts_FISH%>% 
  filter(Depth %in% c("DCM") & Domain %in% c("EUB","ARCH"))%>%
  left_join(.,epipelagic_total_EUB_ARCH[,c("StationName","total")], by = "StationName") %>%
  mutate(proportion = (conc.mn / total)) -> epipelagic_EUB_ARCH_proportion
#summarize by regions
epipelagic_EUB_ARCH_proportion%>%
  group_by(Region,Domain) %>% 
  summarise (mean.proportion = mean(proportion),
             se.proportion = se(proportion), 
             num.samples = length(proportion)) -> epipelagic_EUB_ARCH_proportion_by_regions






#plot the different regions with siginificance test
Boxplot_epipelagic_DAPI <- ggplot(subset(counts_DAPI, Depth %in% c("DCM")), aes(x= Region, y = conc.mn, label = StationName))+
  geom_boxplot(aes(fill = Region))+
  geom_jitter(width = 0.2, alpha = 0.3)+
  geom_signif(comparisons = list(c("EGC", "WSC"),c("EGC","N"),c("N","WSC")), 
              map_signif_level=TRUE, test = "wilcox.test")+
  scale_y_log10(name = "cell conc. [log10(Cells/mL)]")+
  scale_fill_manual(values=c("WSC" = "red", "EGC" = "blue", "N" = "gray")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "bottom")


#calculate mean cell densities 
counts_DAPI %>% 
  filter(Depth %in% c("DCM"))%>%
  group_by(Region) %>% 
  summarise (cell.density.mean = mean(conc.mn),
             cell.density.se = se(conc.mn), 
             num.samples = length(conc.mn)) -> epipelagic_cell_densities_DAPI












#plot DAPi conc. according to stations (perhaps for SI?)
counts_DAPI$StationName <- factor(counts_DAPI$StationName, levels = c("EG1","EG4","N5","N4","N3","HG9","HG7","HG5","HG4","HG2","HG1"))

Boxplot_epipelagic_stations_DAPI <- ggplot(subset(counts_DAPI, Depth %in% c("DCM","EPI")), aes(x= StationName, y = conc.mn, label = StationName))+
  geom_boxplot()+
  geom_jitter(aes(x= StationName, y = conc.mn, label = StationName),width = 0.2)+
  #geom_text()+
  facet_grid(Depth~.)+
  theme(axis.text.x = element_text(angle = 90))+
  theme_plot+
  scale_y_log10(name = "cell conc. [log10(Cells/mL)]")+
  theme(panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5))+
  theme(panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))




###################################
#Bar plot
###################################
#calculate mean abundace by taxa
counts_FISH %>% 
  filter(Depth %in% c("DCM") & Domain %in% taxa)%>%
  group_by(Region,Domain) %>% 
  summarise (cell.density.mean = mean(conc.mn),
             cell.density.se = se(conc.mn), 
             num.samples = length(conc.mn)) -> epipelagic_cell_densities_FISH


#plot the different regions with siginificance test
Boxplot_epipelagic_FISH <- ggplot(epipelagic_cell_densities_FISH, aes(x= Region, y = cell.density.mean, label = Domain))+
  geom_col(aes(fill = Region))+
  geom_errorbar(data=epipelagic_cell_densities_FISH,aes(ymin = cell.density.mean-cell.density.se, ymax = cell.density.mean+cell.density.se), width = 0.2) +   
  facet_grid(.~Domain)+
  #scale_y_log10(name = "cell conc. [log10(Cells/mL)]")+
  scale_fill_manual(values=c("WSC" = "red", "EGC" = "blue", "N" = "gray")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "bottom")







###################################
## NMDS
###################################
counts_FISH[,c("StationName","Depth", "Domain", "Region","conc.mn")] %>% spread(Domain,conc.mn) %>% filter (Depth =="DCM")-> FISH.wide.SRF

#Calculate NMDS of only DCM depths
#list all taxa (excluding EUB and ARCH)
taxa <- c("ALT", "BACT", "CFX", "CREN", "DELTA", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202","SAR324", "SAR406", "VER")

FISH.wide.SRF$total <- FISH.wide.SRF$EUB+FISH.wide.SRF$ARCH

FISH.wide.SRF%>% mutate_at(vars(-c(StationName, Depth, Region, total)), funs(. / total)) -> FISH.wide.SRF.ra

all_metaMDS.SRF <- metaMDS(FISH.wide.SRF.ra[,taxa], maxit= 100, trace=TRUE)

data.scores.SRF <- as.data.frame(scores(all_metaMDS.SRF))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores.SRF$site <- FISH.wide.SRF.ra$StationName  # create a column of site names, from the rownames of data.scores
data.scores.SRF$Depth <- FISH.wide.SRF.ra$Depth  #  add the grp variable created earlier
data.scores.SRF$Region <- FISH.wide.SRF.ra$Region

species.scores.SRF <- as.data.frame(scores(all_metaMDS.SRF, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores.SRF$species <- rownames(species.scores.SRF)  # create a column of species, from the rownames of species.scores

#Plot NMDS
NMDS_DCM <- ggplot() + 
  geom_text(data=species.scores.SRF,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores.SRF,aes(x=NMDS1,y=NMDS2, colour = Region, shape= Depth),size=5) + # add the point markers
  geom_text(data=data.scores.SRF,aes(x=NMDS1,y=NMDS2,label=site),size=3,nudge_y =-0.01) +  # add the site labels
  scale_colour_manual(values=c("WSC" = "red", "EGC" = "blue", "N" = "gray")) +
  coord_equal() +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "bottom")


#combined plot
plot_grid(Boxplot_epipelagic_DAPI, NMDS_DCM, labels = c("A","B"), ncol = 2, align = "hv")

ggsave("./figures/Figure-2.pdf", 
       plot = ggplot2::last_plot(),
       scale = 1,
       units = "cm",
       #width = 17.8,
       #height = 17.4,
       dpi = 300)

     