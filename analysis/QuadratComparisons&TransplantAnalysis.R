###Comparisons of zostera and smithora biomass at edge and interior sites
getwd()
setwd("C:/Users/Gwen/Desktop/Rpractice")

compare <- read.csv("choked_seagrass_smithora.csv")

library(ggplot2)
library(dplyr)

edge <- filter(compare, site == "choked_edge_wolf")
inner <- filter(compare, site == "choked_inner_wolf")
compare_edge <- bind_rows(inner,edge)

head(inneredge)

###Making data nicer for plotting
head(compare)
View(compare)
compare_edge2 <- compare_edge
compare_edge2$site <- as.character(compare_edge2$site)
compare_edge2$site[compare_edge$site == "choked_edge_wolf"] <- "Edge"
compare_edge2$site <- as.character(compare_edge2$site)
compare_edge2$site[compare_edge$site == "choked_inner_wolf"] <- "Interior"

##Shoot density
comparisonsshoots <- ggplot(compare_edge2, aes(x = site, y = number.shoots.quadrat)) + geom_boxplot(size = 1) + theme_bw() + theme(panel.grid = element_blank()) + ylab(expression("Number of shoots")) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(legend.position = "none") + theme(strip.text.y = element_text(face = 4)) + xlab(expression("Location")) + theme(axis.title.x = element_text(vjust = -1, size =14, color = "black")) + theme(axis.title.y = element_text(vjust = -1, size = 14, color = "black")) + theme(axis.text.x = element_text(color = "black", vjust = -1, size = 12)) + theme(axis.text.y = element_text(color = "black", size = 8))

print(comparisonsshoots)
dev.off()
##shoot biomass
comparisonsbiomass <- ggplot(compare_edge, aes(x = site, y = total.dry.macro.weight)) + geom_boxplot(size = 1)
print(comparisonsbiomass)
##Smithora Biomass
comparisonssmithora <- ggplot(compare_edge2, aes(x = site, y = macroepiphyte.weight)) + geom_boxplot(size = 1) + theme_bw() + theme(panel.grid = element_blank()) + ylab(expression(paste(italic("Smithora"), "biomass (g)"))) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(legend.position = "none") + theme(strip.text.y = element_text(face = 4)) + xlab(expression("Wolf Location")) + theme(axis.title.x = element_text(vjust = -1, size =14, color = "black")) + theme(axis.title.y = element_text(vjust = -1, size = 14, color = "black")) + theme(axis.text.x = element_text(color = "black", vjust = -1, size = 12)) + theme(axis.text.y = element_text(color = "black", size = 8))
print(comparisonssmithora)

hist(compare_edge$number.shoots.quadrat)
hist(edge$total.dry.macro.weight)
shapiro.test(edge$total.dry.macro.weight)

hist(inner$total.dry.macro.weight)

shapiro.test(inner$total.dry.macro.weight)
##Proceed with linear model
model_smith <- lm(compare_edge$macroepiphyte.weight ~ compare_edge$site)

anova(model_smith)
summary(model_smith)

model_shoots <- lm(compare_edge$number.shoots.quadrat ~ compare_edge$site)
anova(model_shoots)

###establish with figure that there is significantly higher smithora biomass at the edge of the meadow
##Significantly higher density of shoots at the edge vs. interior

model_biomass <- lm(compare_edge$total.dry.macro.weight ~ compare_edge$site)
anova(model_biomass)
##edge has denser shoots, more biomass, more shoots biomass

compare_edge$shoot_smithora_ratio <- c(compare_edge$macroepiphyte.weight/compare_edge$total.dry.macro.weight)

comparisonsratio <- ggplot(compare_edge, aes(x = site, y = shoot_smithora_ratio )) + geom_boxplot(size = 1)
print(comparisonsratio)

model_ratio <- lm(compare_edge$shoot_smithora_ratio ~ compare_edge$site)
anova(model_ratio)

##but the ratio of smithora to shoot biomass is still higher for the edge. Something is happening here which is why we did a reciprocal transplant.

smithora <- read.csv('Smithora_Plot_Dataset.csv')
head(smithora)
View(smithora)

hist(smithora$Smithora)
shapiro.test(smithora$Smithora)


Biomass <- read.csv("Smith.csv")
head(Biomass)
View(Biomass)
Biomasstrans <- mutate(Biomassexp, transformsmith = log((biomass)+2))
Biomassexp <- slice(Biomass, 5:19)
hist(Biomasstrans$transformsmith)
shapiro.test(Biomasstrans$transformsmith)
str(Biomasstrans)
experiments <- filter(Biomassexp, Type == "experiment")
##analysis with controls
model_smith <- lm(Biomasstrans$transformsmith ~ Biomasstrans$Source + Biomasstrans$Type + Biomasstrans$Source*Biomasstrans$Type)
##analysis with just experiments
model_smith2 <- lm(experiments$biomass~ experiments$Source)
summary(model_smith2)
boxplot(experiments$biomass,experiments$Source)
summary(model_smith)
anova(model_smith)

##significant effect of source and experiment (whether it was moved out of its environment or kept in place on determining smithora biomass)
library(ggplot2)

comparisons <- ggplot(smithora2)
View(smithora2)

comparisons <- ggplot(smithora2, aes(x = Type, y = Smithora)) + geom_boxplot(size = 1) + facet_wrap(~Location)+ theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression(paste(italic("Smithora"), "biomass (g)"))) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(legend.position = "none") + theme(strip.text.y = element_text(face = 4)) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, colour = "black" )) + theme(axis.title.x = element_blank())

?element_text


?xlab

print(comparisons)

##Making data labels nice for plots

smithora2 <- smithora
smithora2$Location <- as.character(smithora2$Location)
smithora2$Location[smithora$Location == "edge"] <- "Meadow Edge"
smithora2$Location <- as.character(smithora2$Location)
smithora2$Location[smithora$Location == "interior"] <- "Meadow Interior"
smithora2$Type <- as.character(smithora2$Type)
smithora2$Type[smithora$Type == "control"] <- "Placed in original location"
smithora2$Type <- as.character(smithora2$Type)
smithora2$Type[smithora$Type == "experiment"] <- "Transplanted"

smithora2$Type <- as.character(smithora2$Type)
smithora2$Type[smithora$Type == "ambient"] <- "Surrounding shoots"
