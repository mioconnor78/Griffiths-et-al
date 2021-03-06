## Griffiths et al code Mary O'Connor
## June 29 2017


### load packages

library(tidyverse)
library(ggplot2)
library(dplyr)
library(reshape2)
library(plyr)
library(stringr)
library(nlme)
library(MuMIn)

### load data
smith.bio <- read_csv("./data/Smith.csv")
compare <- read.csv("./data/choked_seagrass_smithora.csv")
ang <- read.csv("./data/macrophyte_biomass_zeros_20170709.csv", stringsAsFactors=F)
sites <- read.csv("./data/choked_sites.csv")
smith.expt <- read_csv("./data/Transplant_raw_smithora_wt.csv")
grazertraits <- read_csv("./data/grazertraitsmaster.csv")
grazers <- read_csv("./data/Grazers_EdgeVSInterior_V1_Aug09_2017.csv")

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

## function to estimate sds for boxplots
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

### FIGURE 1: SPATIAL TRENDS IN SMITHORA LOAD

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

ang4 <- merge(ang.smith, ang.zostera, by = c("blade_sample_id", "year","date", "site", "distance_m"), all = TRUE)
ang5 <- ang4[,c("site", "group.y", "date", "blade_sample_id", "distance_m", "macrophyte.x","final_dry_g.x", "macrophyte.y", "final_dry_g.y")]
ang5 <- ang4 %>%
  select(c(site, group.y, date, blade_sample_id, distance_m, macrophyte.x,final_dry_g.x, macrophyte.y, final_dry_g.y)) %>%
  filter(macrophyte.y == "zostera") %>%
  replace_na(list(macrophyte.x = "smithora", final_dry_g.x = "0")) %>%
  mutate(ratio = as.numeric(final_dry_g.x) / final_dry_g.y)

ang5$ratio <- ang5$final_dry_g.x / ang5$final_dry_g.y
ang5$date <- as.Date(ang5$date, format="%d-%b")

summary<-ddply(ang5, .(site, date), summarise, length((ratio)))
hist(log(ang5$ratio))

### Figure 1C: Smithora load in edge v interior sites  
smith.load <- ggplot(ang5, aes(x = group.y, y = ratio)) + 
  geom_jitter(position = position_jitter(width=.2), size=3, colour = "gray60") +
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", fill = "transparent") + 
  scale_y_continuous(name = "Smithora (g dwt) / Zostera (g dwt)", limits = c(-0.15, 4), breaks = c(0,1,2,3,4), labels = c(0,1,2,3,4)) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  xlab(expression("Location")) + 
  geom_text(label = "C", x = 2.4, y = 4) +
  geom_text(label = "C", x = 2.4, y = 4) +
  theme(axis.title.x = element_text(vjust = -1, size = 12)) + 
  theme(axis.title.y = element_text(vjust = -1, size = 12))
smith.load
ggsave("smith.load.ang.jpg", plot = smith.load, width = 3, height = 3)

hist(log(ang5$ratio+0.01))
mod1 <- glm(ang5$ratio ~ ang5$group.y, family = quasipoisson) ## should probably be done with site as a ranef...

ang5$site <- as.factor(ang5$site)
mod3 <- lme(I(log(ratio+0.01)) ~ group.y, random = ~1|site, data = ang5)
mod3.1 <- lm(I(log(ratio+0.01)) ~ group.y, data = ang5)

mod3.2 <- lm(I(log(ratio+0.01)) ~ group.y*site, data = ang5)
model.sel(mod3.1, mod3.2)

mean(ang5[(ang5$site == "wolf"),]$ratio)
mean(ang5[(ang5$site == "inner ang"),]$ratio)


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

## FIGURE 2: ZOSTERA, SMITHORA and GRAZER ABUNDANCE at transplant sites before experiment

## A) Shoot density (moved to appendix)
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

## analysis
hist(compare_edge2$number.shoots.quadrat)
mod1 <- lm(number.shoots.quadrat ~ site, data = compare_edge2)
anova(mod1)

## A) Zostera above ground biomass (dry weight)
shoot.dwt <- ggplot(compare_edge2, aes(x = site, y = total.dry.macro.weight)) + 
  geom_jitter(position = position_jitter(width=.2), size=3, colour = "gray60") +
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", fill = "transparent") + 
  scale_y_continuous(name = "Zostera shoot dry weight (g)", labels=c("0","10","20", "30"), limits = c(0, 30)) +
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  xlab(expression("Location")) + 
  geom_text(label = "A", x = 2.4, y = 30) +
  theme(axis.title.x = element_text(vjust = -1, size = 12)) + 
  theme(axis.title.y = element_text(vjust = -1, size = 12))

ggsave("shoot.dwt.jpg", plot = shoot.dwt, width = 3, height = 3)

hist(compare_edge2$total.dry.macro.weight)
mod1 <- lm(total.dry.macro.weight ~ site, data = compare_edge2)
anova(mod1)

## B) Smithora load (dry weight)
smithora <- ggplot(compare_edge2, aes(x = site, y = macroepiphyte.weight/total.dry.macro.weight)) + 
  geom_jitter(position = position_jitter(width=.2), size=3, colour = "gray60") +
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", fill = "transparent") + 
  scale_y_continuous(name = expression(paste(italic("Smithora"), " (g dry wt / g Zostera dry wt)")), limits = c(0, 0.5)) +
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  xlab(expression("Location")) + 
  geom_text(label = "B", x = 2.4, y = 0.5) +
  theme(axis.title.x = element_text(vjust = -1, size = 12)) + 
  theme(axis.title.y = element_text(vjust = -1, size = 12))

ggsave("smithora.jpg", plot = smithora, width = 3, height = 3)

## D) GRAZER ABUNDANCE -- placeholder, waiting for final plot labels from Gwen
head(grazers)
grazers[is.na(grazers)] <- 0
grazers1 <- melt(grazers, id = c(1,2,3,4))
names(grazers1) <- c("Site", "Source", "Plot", "Size", "Species", "Abundance")
grazers2 <- grazers1 %>%
  filter(Site != "0")
grazers3 <- ddply(grazers2, .(Site, Source, Plot, Species), summarise, sum(Abundance))
grazers4 <- grazers3 %>%
  filter(Species != "byozoan 1cm") %>%
  filter(Species != "NA") %>%
  filter(Species != "stalked jelly") %>%
  filter(Species != "flat worm") %>%
  filter(Species != "anemone") %>%
  filter(Species != "fruit loop") %>%
  filter(Species != "Dungeness") %>%
  filter(Plot != "ASU") 
names(grazers4) <- c("Site", "Source", "Plot", "Species", "Abundance")
grazers5 <- ddply(grazers4, .(Site, Source, Plot), summarise, sum(Abundance))
names(grazers5) <- c("Site", "Source", "Plot","Abundance")

grazers <- ggplot(grazers5, aes(x = Source, y = Abundance)) + 
  geom_jitter(position = position_jitter(width=.2), size=3, colour = "gray60") +
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", fill = "transparent") + 
  scale_y_continuous(name = "Epifaunal Grazer Density (N / XX)") + #limits = c(0, 0.5)
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  xlab(expression("Location")) + 
  geom_text(label = "C", x = 2.4, y = 650) +
  theme(axis.title.x = element_text(vjust = -1, size = 12)) + 
  theme(axis.title.y = element_text(vjust = -1, size = 12))

ggsave("grazers.jpg", plot = grazers, width = 3, height = 3)

## now once we know the treatments for each plot, we can plot total abundance


## FIGURE 3: Smithora results post experiment - remake this figure.

head(smith.expt)
smith.expt$Source <- revalue(smith.expt$Source, c("Edge"="Edge (high Smithora)",  "Interior"="Interior (low Smithora)"))
smith.expt$Source <- as.factor(smith.expt$Source)
smith.expt$Type <- revalue(smith.expt$`Category in transplant`, c("Ambient"="Unmanip.", "Experiment"="Transplanted", "Control"="Control"))
smith.expt$Type <- as.factor(smith.expt$Type)
smith.expt$ratio <- smith.expt$`Ratio smithora/leaves`

experiment <- ggplot(smith.expt[(smith.expt$Type != "Unmanip."),], aes(x = Type, y = ratio)) + 
  geom_jitter(position = position_jitter(width=.2), size=3, colour = "gray60") +
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", fill = "transparent") + 
  facet_wrap(~Source) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  theme(strip.background = element_rect(fill="white")) + 
  scale_y_continuous(name = expression(paste(italic("Smithora"), " (g dry wt / g Zostera dry wt)")), limits = c(-0.1, 1), breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  xlab("") +
  labs(title = "Site of Experimental Planting") 

experiment
ggsave("experiment.jpg", plot = experiment, width = 6, height = 3)

## Experimental analysis: 

## first analyze control vs unmanipulated
hist((smith.expt[(smith.expt$Type != "Transplanted"),]$ratio))
mod1 <- lm((smith.expt[(smith.expt$Type != "Transplanted"),]$ratio) ~ smith.expt[(smith.expt$Type != "Transplanted"),]$Type + smith.expt[(smith.expt$Type != "Transplanted"),]$Source)
anova(mod1)

## then analyze control v transplanted
hist(log(smith.expt[(smith.expt$Type != "Unmanip."),]$ratio + 1))
mod2 <- lm(log(smith.expt[(smith.expt$Type != "Unmanip."),]$ratio + 1) ~ smith.expt[(smith.expt$Type != "Unmanip."),]$Type * smith.expt[(smith.expt$Type != "Unmanip."),]$Source)
anova(mod2)

