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
library(cowplot)

###################################
## Bar plot with rel.abundance
###################################

counts.rel.ab.dcm <- merge(fish.rel, dapi.rel ,by =c("StationName", "Domain", "Depth"))
counts.rel.ab.dcm <- subset(counts.rel.ab.dcm, Depth %in% c("DCM"))
counts.rel.ab.dcm.egc <- subset(counts.rel.ab.dcm, Region.x %in% c("EGC"))
counts.rel.ab.dcm.wsc <- subset(counts.rel.ab.dcm, Region.x %in% c("WSC"))

counts.rel.ab.dcm.egc <- counts.rel.ab.dcm.egc[,c("Domain", "fish.mn", "dapi.mn")]
counts.rel.ab.dcm.egc <- aggregate(.~Domain, data=counts.rel.ab.dcm.egc, FUN = mean)

counts.rel.ab.dcm.wsc <- counts.rel.ab.dcm.wsc[,c("Domain", "fish.mn", "dapi.mn")]
counts.rel.ab.dcm.wsc <- aggregate(.~Domain, data=counts.rel.ab.dcm.wsc, FUN = mean)

counts.rel.ab.dcm.egc$proportion.egc <- counts.rel.ab.dcm.egc$fish.mn*100/counts.rel.ab.dcm.egc$dapi.mn
counts.rel.ab.dcm.wsc$proportion.wsc <- counts.rel.ab.dcm.wsc$fish.mn*100/counts.rel.ab.dcm.wsc$dapi.mn
counts.rel.ab.dcm.egc$Region <- c("EGC")
counts.rel.ab.dcm.wsc$Region <- c("WSC")

bar_abs_plot <- merge(counts.rel.ab.dcm.egc, counts.rel.ab.dcm.wsc ,by =c("Domain"))
bar_abs_plot<- bar_abs_plot[,c("Domain", "proportion.egc", "proportion.wsc")]

bar_abs_plot <- bar_abs_plot %>%
  gather(Total, Value, -Domain)

bar_abs_p <- ggplot(bar_abs_plot, aes(x = Domain, y = Value, fill = Total)) +
  geom_col(position = "dodge")+
  scale_x_discrete(breaks=c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER"),
                   labels=c("Alteromonadales", "Archaea", "Bacteroidetes", "Chloroflexi", "Thaumarcheota", "Deltaproteobacteria", "Bacteria", "Gammaproteobacteria", "Opitutae", "Polaribacter", "Rhodobacteraceae", "SAR11", "SAR202 clade", "SAR324 clade", "Marinimicrobia (SAR406 clade)", "Verrucomicrobia"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =14), axis.text.y = element_text(size=14), axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65))+
  scale_fill_manual(name ="Region", labels=c("Ice-covered", "Ice-free"), values = c("blue","red"))+
  theme(legend.text = element_text(size=13),
        legend.position = "none",
        legend.title = element_text(size=13),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1),
        panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5),
        panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))

###################################
## Bar plot with cell concentration
###################################
counts.abs.dcm.egc <- counts.rel.ab.dcm.egc 
colnames(counts.abs.dcm.egc)[2] <- 'fish.egc.mn'
counts.abs.dcm.wsc <- counts.rel.ab.dcm.wsc
colnames(counts.abs.dcm.wsc)[2] <- 'fish.wsc.mn'

bar_cellcon_plot <- merge(counts.abs.dcm.egc, counts.abs.dcm.wsc ,by =c("Domain"))
bar_cellcon_plot<- bar_cellcon_plot[,c("Domain", "fish.egc.mn", "fish.wsc.mn")]

bar_cellcon_plot <- bar_cellcon_plot %>%
  gather(Total, Value, -Domain)

bar_con_plot <- ggplot(bar_cellcon_plot, aes(x = Domain, y = Value, fill = Total)) +
  geom_col(position = "dodge")+
  scale_x_discrete(breaks=c("ALT", "ARCH", "BACT", "CFX", "CREN", "DELTA", "EUB", "GAM", "OPI", "POL", "ROS", "SAR11", "SAR202", "SAR324", "SAR406", "VER"),
                   labels=c("Alteromonadales", "Archaea", "Bacteroidetes", "Chloroflexi", "Thaumarcheota", "Deltaproteobacteria", "Bacteria", "Gammaproteobacteria", "Opitutae", "Polaribacter", "Rhodobacteraceae", "SAR11", "SAR202 clade", "SAR324 clade", "Marinimicrobia (SAR406 clade)", "Verrucomicrobia"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =14), axis.text.y = element_text(size=14), axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_y_continuous(expand = c(0, 0),breaks=c(seq(0,9e+05,1e+05)), labels = function(x) format(x, scientific = TRUE))+
  scale_fill_manual(values = c("blue","red"))+
  theme(legend.text = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", size = 1),
        panel.grid.minor = element_line(linetype = 'dotted', color = 'white', size= 0.5),
        panel.grid.major = element_line(linetype = 'dotted', color = 'white', size= 0.5))
bar_con_plot
                  
both_plots <- cowplot::plot_grid(bar_con_plot, bar_abs_p, 
          labels = c("A", "B"), hjust = -45,
          ncol = 2, nrow = 1)
both_plots

save.image("card_fish_water_column_ps99.Rdata")
