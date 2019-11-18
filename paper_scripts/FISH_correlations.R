###################################
## Load required libraries
###################################
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(PerformanceAnalytics); packageVersion("PerformanceAnalytics")
library(tidyverse); packageVersion("tidyverse")
library(broom); packageVersion("broom")

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

#define function for correlation 
cor_fun <- function(df) cor.test(df$Abund, df$value, method = "pearson",conf.level = 0.95) %>% tidy()
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
  filter(FISH_conc > 0) %>% 
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

##################################
# Env. parameter coorelation test 
##################################
## import environmental metadata 
metadata.raw <- read.csv("./Rwork/PS99_metadata.csv", sep = ",", dec = ".", header = TRUE)

#list environmental parameters that need to be scaled
env.par <- c("Temperature","Salinity", 
             "Chla_fluor", "d.NO3.", "d.NH4.", "d.PO4.", "d.SiO3.")

#remove station SV2 from the dataset, drop rows with NA and scale the env. parameters 
metadata.scaled <- drop_na(metadata.raw) %>% 
  mutate_at(env.par, scale_par)

#Calculate the Pearson’s correlation coefficient in a matrix
chart.Correlation(metadata.scaled[,env.par], method = "pearson", histogram=TRUE, pch=19)

##################################
# correlation between env. par. and counts
##################################
#list all taxa (excluding EUB and ARCH)
taxa <- c("ALT", "BACT", "CFX", "CREN", "DELTA", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202","SAR324", "SAR406", "VER")

#select only surface samples
metadata.scaled %>%
  filter (Depth == "SRF") -> metadata.scaled.SRF

#generate wide abundance table and env. par.
data_all <- counts_all%>% 
  filter(Depth %in% c("SRF"))%>%
  group_by(Region, StationName, Domain)%>%
  select(Region, StationName, Domain, FISH.conc.mn)  %>% 
  spread(Domain, FISH.conc.mn)  %>%
#drop samples with no env. data
  filter(StationName %in% metadata.scaled.SRF$StationName) %>% 
#merge counts and environmental data
left_join(., metadata.scaled.SRF[,c("StationName",env.par)] , by = "StationName")


#generate long table and nest the table according to taxa and env. variable
data_nest <- gather(data_all, Domain, Abund, taxa)%>%
  gather(variable, value, env.par) %>% 
  group_by(Domain, variable) %>% nest()

#nested correlations tests
data_nest <- mutate(data_nest, model = map(data, cor_fun))

#summary table
corr_pr.sign <- select(data_nest, -data) %>% unnest() %>% 
  #extract only significant results 
              filter(p.value < 0.05)

write.csv(corr_pr.sign, "./Rwork/corr_table.csv")
