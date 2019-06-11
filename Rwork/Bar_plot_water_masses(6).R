#####################################
#working directory
#####################################

setwd("~/CARD-FISH/CARD-FISH_water_project/")

load("card_fish_water_column_ps99.Rdata")

###################################
## Load required libraries
###################################

library(ggplot2)
library(dplyr)


###################################
## Pie chart
###################################

source("~/PhD/Projects/Color_palettes.R") 

counts.rel.ab.wenv <- counts_FISH %>% left_join(env, by=c("StationName", "Depth"))

counts_pie <- counts.rel.ab.wenv[,c("Water_mass", "Domain", "conc.mn")]
counts_pie <- counts_pie[!counts_pie$Domain %in% c("EUB"), ]
counts_pie <- counts_pie[!counts_pie$Domain %in% c("SAR11"), ]

counts_pie <- transform(counts_pie, Tot.conc.mn = ave(conc.mn, Water_mass, FUN = sum))
counts_pie <- transform(counts_pie, percentage = 100 * conc.mn/Tot.conc.mn)

bar_AW <- ggplot(subset(counts_pie, Water_mass %in% c("AW")), aes(x = Domain, y = percentage, fill = Domain)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = tol21rainbow)+
  theme_plot+
  theme(axis.title.x = element_blank())+theme(axis.title.y = element_blank())+theme(axis.ticks.x = element_blank())+theme(axis.ticks.y = element_blank())
bar_AW 

save.image("card_fish_water_column_ps99.Rdata")
                     