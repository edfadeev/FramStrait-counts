###################################
## Load required libraries
###################################
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")



library(ggsignif); packageVersion("ggsignif")


library(PerformanceAnalytics); packageVersion("PerformanceAnalytics")

library(tidyverse); packageVersion("tidyverse")
library(broom); packageVersion("broom")
library(fs); packageVersion("fs")
library(lubridate); packageVersion("lubridate")

###################################
## defined functions
###################################
#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

#microscope calculation factor
calc.factor <- 99515.5458411807

###################################
## Import counts and calculate cell concentrations
###################################
#raw counts
raw.counts <- read.csv("./Rwork/FOV_all_groups_SH.csv", sep = ",", dec = ".", header = TRUE)

#split sample name and remove SV2 station and ABYSS depth and calculate concentrations per FOV
raw.counts.SH <- raw.counts %>% 
  separate(SAMPLE_NAME, c("StationName", "Depth", "Domain"),"-") %>% 
  filter(StationName != "SV2" & Depth != "ABYSS") %>% 
  mutate(DAPI_conc= (DAPI_Nr_Set*calc.factor)/Volume,
         FISH_conc = (DAPI_Nr_SubSet*calc.factor)/Volume)


#calculate concetrations per sample and add regions
counts_all <-  raw.counts.SH %>% 
  mutate(Region = ifelse(StationName %in% c("EG1","EG4"), "EGC", 
                         ifelse(StationName %in% c("N3","N4","N5"), "N", "WSC")),
         Depth = factor(Depth, levels = c("SRF","EPI","MESO","BATHY"))) %>% 
          group_by(Region, StationName, Depth, Domain) %>% 
          summarise(DAPI.conc.mn = mean(DAPI_conc),
             DAPI.conc.md =  median(DAPI_conc), 
             DAPI.conc.sd = sd(DAPI_conc),
             FISH.conc.mn = mean(FISH_conc),
             FISH.conc.md =  median(FISH_conc), 
             FISH.conc.sd = sd(FISH_conc),
             n = length(DAPI_conc))
  

###################################
##Plot EUB and ARCH profiles by region
###################################
DAPI_vertical_boxplot <- ggplot(counts_all, aes(x= Region, y = DAPI.conc.mn))+
  geom_boxplot(aes(fill = Region))+
  facet_grid(Depth~.)+
  geom_jitter(aes(x= Region, y = DAPI.conc.mn),width = 0.2, alpha = 0.3)+
  #geom_signif(comparisons = list(c("EGC", "WSC"),c("EGC","N"),c("N","WSC")), 
  #            map_signif_level=TRUE, test = "wilcox.test")+
  scale_y_log10(name = "cell concentration (cells/mL)")+
  scale_fill_manual(values=c("WSC" = "red", "EGC" = "blue", "N" = "gray")) +
  theme_bw()+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")

counts_all%>% 
  filter(Domain %in% c("EUB","ARCH"))%>%
  group_by(Region, StationName, Depth, Domain) -> EUB_ARCH_absolute

p <- list()
for (i in c("EUB","ARCH")){
  p[[i]] <-  ggplot(counts_all[counts_all$Domain== i,], aes(x= Region, y = FISH.conc.mn))+
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


###################################
##Check the coverage of EUB and ARCH probes in comparison to DAPI
###################################
#coverage by probes
counts_all%>% 
  filter(Domain %in% c("EUB","ARCH"))%>%
  group_by(Region, StationName, Depth, Domain) %>% 
  summarise(proportion = FISH.conc.mn/DAPI.conc.mn)-> EUB_ARCH_prop

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

