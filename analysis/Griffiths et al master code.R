## Griffiths et al code Mary O'Connor
## June 29 2017


### load packages

library(tidyverse)
library(ggplot2)
library(dplyr)
library(reshape2)
library(plyr)
library(stringr)


### load data
smith.bio <- read_csv("./data/Smith.csv")
compare <- read.csv("./data/choked_seagrass_smithora.csv")
ang <- read.csv("./data/macrophyte_biomass_zeros_20170709.csv", stringsAsFactors=F)
sites <- read.csv("./data/choked_sites.csv")

### compare
## from meadow surveys, identify edge and inner samples
edge <- filter(compare, site == "choked_edge_wolf")
inner <- filter(compare, site == "choked_inner_wolf")
compare_edge <- bind_rows(inner,edge)

###Making data nicer for plotting
head(compare)
View(compare)
compare_edge2 <- compare_edge
compare_edge2$site <- as.character(compare_edge2$site)
compare_edge2$site[compare_edge$site == "choked_edge_wolf"] <- "Edge"
compare_edge2$site <- as.character(compare_edge2$site)
compare_edge2$site[compare_edge$site == "choked_inner_wolf"] <- "Interior"

### ANG  - These are ang's survey data. We still have a few questions about how to align epiphyte and shoot biomass. This current code, which I think is right, produces highlights that there are a few (~ 10) samples with epiphyte data but no Zostera data. So we are trying to understand this, and in the mean time will proceed with monitoring data (below). [7/5/2017]

ang2 <- merge(ang, sites, by.x = "site", by.y = "name2")
## this next bit selects the sites we want, and removes a couple of rows that appeared to be duplicates, zeros, or unclearly assigned to a shoot. 
ang <- ang2 %>%
  filter(group.x != "kelp") %>%
  select(-c(notes, oven_date, QAQC_date,QAQC_Initials, Lat, Long, foil_wgt_g, wet_wgt_g, dry_wgt_g)) %>%
  mutate(uniqueID = paste(date, site, blade_sample_id, sep = ".")) %>%
  filter(final_dry_g != "0") %>%
  filter(final_dry_g > "0") %>%
  filter(uniqueID != "29-May.inner sandspit.75") %>%
  mutate(uniqueID2 = paste(uniqueID, final_dry_g, sep = ".")) %>%
  filter(uniqueID2 != "29-May.alcove.38.0.1") %>%
  select(-uniqueID2)

ang.smith <- subset(ang, ang$macrophyte %in% "smithora")
ang.zostera<- subset(ang, ang$macrophyte %in% "zostera")

ang4 <- merge(ang.smith, ang.zostera, by = c("blade_sample_id", "year","date", "site", "distance_m"))
ang5 <- ang4[,c("site", "group.y", "date", "blade_sample_id", "distance_m", "macrophyte.x","final_dry_g.x", "macrophyte.y", "final_dry_g.y")]
ang5$ratio <- ang5$final_dry_g.x / ang5$final_dry_g.y
ang5$date <- as.Date(ang5$date, format="%d-%b")

### Figure 1C: Smithora load in edge v interior sites  
smith.load <- ggplot(ang5, aes(x = group.y, y = ratio)) + 
  geom_point(size = 6, colour = "gray") +
  geom_boxplot(size = 1, fill = "transparent") + 
  scale_y_continuous(name = "Smithora (g dwt) / Zostera (g dwt)", limits = c(0, 4)) +
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  xlab(expression("Location")) + 
  geom_text(label = "C", x = 2.4, y = 4) +
  geom_text(label = "C", x = 2.4, y = 4) +
  theme(axis.title.x = element_text(vjust = -1, size = 12)) + 
  theme(axis.title.y = element_text(vjust = -1, size = 12))

ggsave("smith.load.ang.jpg", plot = shoot.dwt, width = 3, height = 3)

#### Monitoring data
## started with our data, but there isn't an inner site here, really...

shoots2015 <- read.csv("./data/biodiversitysurvey.2015.eelgrass.single.shoots.monitoring.csv")
shoots2015$site_id <- revalue(shoots2015$site, c("goose_east"="GsE", "goose_west"="GsW", "mcmullins_north" = "MMn", "mcmullins_south" = "MMs", "triquet_south" = "TqS", "triquet_north"= "TqN", "choked_south_pigu" = "Cspg", "choked_sandspit" = "Csan"))
View(shoots2015)

shoots <- shoots2015 %>%
  filter(site_id %in% c("Csan", "Cspg")) %>%
  replace(is.na(.), 0) %>%
  mutate(sm.zm = smithora/total.dry.weight)

sm.zm <- ggplot(shoots, aes(x = site_id, y = sm.zm)) + 
  geom_point(size = 6, colour = "gray") +
  geom_boxplot(size = 1, fill = "transparent") + 
  scale_y_continuous(name = "Smithora / Zostera shoot weight (g/g dry wt)", limits = c(0, 2)) +
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  xlab(expression("Location")) + 
  geom_text(label = "A", x = 2.4, y = 15) +
  theme(axis.title.x = element_text(vjust = -1, size = 12)) + 
  theme(axis.title.y = element_text(vjust = -1, size = 12))
sm.zm

## FIGURE 2: ZOSTERA AND SMITHORA ABUNDANCE at transplant sites before experiment

## A) Shoot density
shoot.density <- ggplot(compare_edge2, aes(x = site, y = number.shoots.quadrat)) + 
  geom_point(size = 6, colour = "gray") +
  geom_boxplot(size = 1, fill = "transparent") + 
  scale_y_continuous(name = "Zostera Shoots (N/0.0625 m^2)", labels=c("0","5","10","15"), limits = c(0, 15)) +
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  xlab(expression("Location")) + 
  geom_text(label = "A", x = 2.4, y = 15) +
  theme(axis.title.x = element_text(vjust = -1, size = 12)) + 
  theme(axis.title.y = element_text(vjust = -1, size = 12))

ggsave("shoot.density.jpg", plot = shoot.density, width = 3, height = 3)

## B) Zostera above ground biomass (dry weight)
shoot.dwt <- ggplot(compare_edge2, aes(x = site, y = total.dry.macro.weight)) + 
  geom_point(size = 6, colour = "gray") +
  geom_boxplot(size = 1, fill = "transparent") + 
  scale_y_continuous(name = "Zostera shoot dry weight (g)", labels=c("0","10","20", "30"), limits = c(0, 30)) +
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  xlab(expression("Location")) + 
  geom_text(label = "B", x = 2.4, y = 30) +
  theme(axis.title.x = element_text(vjust = -1, size = 12)) + 
  theme(axis.title.y = element_text(vjust = -1, size = 12))

ggsave("shoot.dwt.jpg", plot = shoot.dwt, width = 3, height = 3)


## C) Zostera above ground biomass (dry weight)
smithora <- ggplot(compare_edge2, aes(x = site, y = macroepiphyte.weight/total.dry.macro.weight)) + 
  geom_point(size = 6, colour = "gray") +
  geom_boxplot(size = 1, fill = "transparent") + 
  scale_y_continuous(name = "Smithora weight (g / dry g SG)", limits = c(0, 0.5)) +
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  xlab(expression("Location")) + 
  geom_text(label = "C", x = 2.4, y = 0.5) +
  theme(axis.title.x = element_text(vjust = -1, size = 12)) + 
  theme(axis.title.y = element_text(vjust = -1, size = 12))

ggsave("smithora.jpg", plot = smithora, width = 3, height = 3)

## FIGURE 3: Smithora results post experiment - remake this figure.
head(smith.bio)
smith.bio$Source <- revalue(smith.bio$Source, c("edge"="Edge (high Smithora)",  "inner"="Interior (low Smithora)"))
smith.bio$Type <- revalue(smith.bio$Type, c("ambient"="Unmanip.", "experiment"="Transplanted", "control"="Control"))     

experiment <- ggplot(smith.bio, aes(x = Type, y = biomass)) + 
  geom_point(size = 6, colour = "gray") +
  geom_boxplot(size = 1, fill = "transparent") + 
  facet_wrap(~Source) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  theme(strip.background = element_rect(fill="white")) + 
  scale_y_continuous(name = expression(paste(italic("Smithora"), " (g dry wt / g Zostera dry wt)")), limits = c(0, 1)) +
  xlab("Shoot Treatment") +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(title = "Site of Experimental Planting") + 

experiment
ggsave("experiment.jpg", plot = experiment, width = 6, height = 3)

## Experimental analysis: 
hist(smith.bio$biomass)
hist(log(smith.bio$biomass))

mod1 <- lm(smith.bio$biomass ~ smith.bio$Type, family = poisson)

### Figure 5: Experimental results


