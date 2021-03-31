

#import counts and calculate proportions for surface waters
source("./paper_scripts/Scripts/import_cell_counts.R")

surface_FISH_proportion<- counts_all%>% 
  group_by(Region, Depth, StationName, Domain)%>%
  mutate(proportion = (FISH.conc.mn / DAPI.conc.mn)*100,
         Domain = case_when(Domain== "SAR11" ~ "SAR11_clade",
         Domain== "ROS" ~ "Rhodobacteraceae",
         Domain== "GAM" ~ "Gammaproteobacteria",
         Domain== "ALT" ~ "Alteromonadales",
         Domain== "BACT" ~ "Bacteroidetes",
         Domain== "POL" ~ "Polaribacter",
         Domain== "VER" ~ "Verrucomicrobiae",
         Domain== "OPI" ~ "Opitutales")) %>%
  select(Region, Depth, StationName, Domain, proportion) %>% 
  filter(Domain %in% c("SAR11_clade","Rhodobacteraceae","Gammaproteobacteria",
                       "Alteromonadales","Bacteroidia","Polaribacter",
                       "Verrucomicrobiae","Opitutales"),
         Depth =="SRF")

prop16S <- read.csv("./paper_scripts/Data/PS99_16s_taxa.csv", sep =",", header=T, row.names = 1)

prop16S_long <- prop16S %>%
  gather(Domain, proportion_16S, Alteromonadales:Verrucomicrobiae) %>%
  separate(Sample, c("StationName", "Depth")," ") %>%
  mutate(Region = ifelse(StationName %in% c("EG1","EG4"), "EGC", 
                         ifelse(StationName %in% c("N3","N4","N5"), "N", 
                                ifelse(StationName %in% c("N3","N4","N5"), "N", "WSC"))),
         Depth = factor(Depth, levels = c("SRF","EPI","MESO","BATHY"))) %>% 
  filter(Domain %in% c("SAR11_clade","Rhodobacteraceae","Gammaproteobacteria",
                       "Alteromonadales","Bacteroidia","Polaribacter",
                       "Verrucomicrobiae","Opitutales"),
         Depth =="SRF")



test <- merge(surface_FISH_proportion, prop16S_long, by = c("Region","StationName","Depth","Domain")) %>% 
      melt(measure.vars = c("proportion","proportion_16S")) %>% 
      mutate(StationName =factor(StationName, 
                          levels = c("EG1","EG4","N5","N4","N3","HG9","HG7","HG5", "HG4","HG2","HG1")))

ggplot(test)+
  geom_col(aes(x=StationName, y = value, fill = Region))+
  facet_grid(variable~Domain)+
  scale_fill_manual(values=c("WSC" = "red", "EGC" = "blue", "N" = "gray")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

ggsave("./proportion_comparison.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)
