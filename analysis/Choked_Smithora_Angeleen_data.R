##Analyzing choked smithora 
getwd()
getwd()
setwd("C:/Users/Gwen/Desktop/Rpractice")

##Import Macrophyte Data

Macrophyte_data <- read.csv("macrophyte.biomass.zeros.csv")
View(Macrophyte_data)
hist(Macrophyte_data$shoot_length_cm)
shapiro.test(Macrophyte_data$shoot_length_cm)

length_site_model <- lm(shoot_length_cm ~ site, data = Macrophyte_data)
weight_site_model <- lm(final_dry_g ~ site, data = macro_smithora)
anova(weight_site_model)
summary(weight_site_model)

anova(length_site_model)
summary(length_site_model)


Macrophyte_data_long <- melt(Macrophyte_data, id.vars = "blade_sample_id")

macro_site <- filter(Macrophyte_data_long, variable == "site")
macro_length <- filter(Macrophyte_data_long, variable == "shoot_length_cm")
macro_width <- filter(Macrophyte_data_long, variable == "shoot_width_cm")

join <- inner_join(macro_length,macro_width, by = "blade_sample_id")
head(join)
join2 <- inner_join(join,macro_site, "blade_sample_id")
join3 <- inner_join(join2,macro_smithora, "blade_sample_id")
finaljoin <- inner_join(join3,macro_zostera,"blade_sample_id")

head(finaljoin)
View(finaljoin)

macro_smithora <- filter(Macrophyte_data, macrophyte == "smithora")





macro_zostera <- filter(Macrophyte_data, macrophyte == "zostera")

smith_length <- ggplot(finaljoin)
smith_length <- smith_length + geom_point((aes(value.y,final_dry_g.x))) + facet_wrap(~value)
print(smith_length)


smith_length <- lm(final_dry_g.x ~ value.y + site.y + value.y*site.y, data = finaljoin)
summary(smith_length)

par(mfrow = c(2, 2))
plot(smith_length)
dev.off()





head(Macrophyte_data)

View(Macrophyte_data_long)

library(dplyr)

Smithora_only <- filter(Macrophyte_data, macrophyte == "smithora")

hist(Smithora_only$final_dry_g)

library(ggplot2)

Chokedsmith <- ggplot(Smithora_only)
Chokedsmith <- Chokedsmith + geom_boxplot(aes(site,final_dry_g))
print(Chokedsmith)
Chokedsmith <- Chokedsmith + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression(paste(italic("Smithora"), "biomass (g)"))) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(legend.position = "none") + theme(strip.text.y = element_text(face = 4)) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, colour = "black" )) + theme(axis.title.x = element_blank())

Seagrass_only <- filter(Macrophyte_data, macrophyte == "zostera")
Chokedzos <- ggplot(Seagrass_only)
Chokedzos <- Chokedzos + geom_boxplot(aes(site,shoot_length_cm))
Chokedzos <- Chokedzos + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression(paste(italic("Zostera"), " length (cm)"))) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(legend.position = "none") + theme(strip.text.y = element_text(face = 4)) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, colour = "black" )) + theme(axis.title.x = element_blank())

print(Chokedzos)

library(reshape2)

smith_length <- ggplot(Macrophyte_data)
smith_lenth <- smith_length + geom_point(aes(shoot_length_cm, final_dry_g)) + facet_wrap(~macrophyte)
print(smith_length)

Seagrass_only_long <- melt(Seagrass_only, id.vars = "blade_sample_id", measure.vars = "final_dry_g","shoot_length_cm")

?melt

str(Seagrass_only)

##Look at relationship between shoots and smithora
Zos_Smith_Metrics <- read.csv("Zos_Smith_Metrics.csv")

View(Zos_Smith_Metrics)
hist(Zos_Smith_Metrics$smithora)
##get residuals

trylinear <- lm(smithora ~ LAI + site + LAI*site, data = Zos_Smith_Metrics)
residuals(trylinear)
visreg(trylinear)
trylinear_res <- residuals(trylinear)
plot(Zos_Smith_Metrics$LAI, trylinear_res)

summary(trylinear)

View(Zos_Smith_Metrics)
hist(Zos_Smith_Metrics$LAI)

shapiro.test(Zos_Smith_Metrics$LAI)

Zos_Smith_Metrics$Transform <- ((Zos_Smith_Metrics$smithora)+2)
hist(Zos_Smith_Metrics$smithora)

hist(Zos_Smith_Metrics$Transform)

library(nlme)

Zos_Smith_Metrics$standard <- (Zos_Smith_Metrics$smithora/Zos_Smith_Metrics$LAI)
hist(Zos_Smith_Metrics$standard)

##Viewing standardized smithora biomass. See if Smithora biomass scaled to how much space available to colonize is actually varying by site
library(ggplot2)
ChokedSmithStandard <- ggplot(Zos_Smith_Metrics)
ChokedSmithStandard <- ChokedSmithStandard + geom_boxplot(aes(site,Transform)) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, colour = "black" ))
print(ChokedSmithStandard)
ChokedSmithLAI <- ggplot(Zos_Smith_Metrics)
ChokedSmithLAI <- ChokedSmithLAI + geom_boxplot(aes(site,LAI))
print(ChokedSmithLAI)


modelLAIgaus <- glm(LAI ~ site, family = gaussian(link="identity"), data = Zos_Smith_Metrics)
plot(modelLAIgaus)
summary(modelLAIgaus)


modelLAIgamma <- glm(LAI ~ site, family = Gamma(link="inverse"), data = Zos_Smith_Metrics)
summary(modelLAIgamma)


modelLAIinverse <- glm(LAI ~ site, family = inverse.gaussian(link="1/mu^2"), data = Zos_Smith_Metrics)
summary(modelLAIinverse)

###looks like LAI does not vary significantly by site, barely significant at outer sandspit. This is actually pretty cool. 

modelstandardgaus <- glm(Transform ~ site, family = gaussian(link="identity"), data = Zos_Smith_Metrics)
summary(modelstandardgaus)
modelstandardgamma <- glm(Transform ~ site, family = Gamma(link="log"), data = Zos_Smith_Metrics)
summary(modelstandardgamma)

modelstandardinverse <- glm(Transform ~ site, family = inverse.gaussian(link="1/mu^2"), data = Zos_Smith_Metrics)
summary(modelstandardinverse)
###when you standardize smithora to LAI and use the length, width, and blade number to explain variation only one site remains different from the others and that is NORTH PIGU
library(AICcmodavg)
##Basically outer sandspit is the only site that varies significantly, and North PIGU is the only site that varies significantly based on Smithora standardized to the LAI it comes from
###With all terms
smithmodelstandardgaus <- glm(Transform ~ site + LAI + LAI*site, family = gaussian(link="identity"), data = Zos_Smith_Metrics)
summary(smithmodelstandardgaus)
smithmodelstandardgamma <- glm(Transform ~ site + LAI + LAI*site, family = Gamma(link="inverse"), data = Zos_Smith_Metrics)
summary(smithmodelstandardgamma)
smithmodelstandardinverse <- glm(Transform ~ site + LAI + LAI*site, family = inverse.gaussian(link="1/mu^2"), data = Zos_Smith_Metrics)
summary(smithmodelstandardinverse)
##With jsut site
sitesmithmodelstandardgaus <- glm(Transform ~ site, family = gaussian(link="identity"), data = Zos_Smith_Metrics)
summary(sitesmithmodelstandardgaus)
sitesmithmodelstandardgamma <- glm(Transform ~ site, family = Gamma(link="inverse"), data = Zos_Smith_Metrics)
summary(sitesmithmodelstandardgamma)
sitesmithmodelstandardinverse <- glm(Transform ~ site, family = inverse.gaussian(link="1/mu^2"), data = Zos_Smith_Metrics)
summary(sitesmithmodelstandardinverse)
##With just LAI
LAIsmithmodelstandardgaus <- glm(Transform ~ LAI, family = gaussian(link="identity"), data = Zos_Smith_Metrics)
summary(LAIsmithmodelstandardgaus)
LAIsmithmodelstandardgamma <- glm(Transform ~ LAI, family = Gamma(link="inverse"), data = Zos_Smith_Metrics)
summary(LAIsmithmodelstandardgamma)
LAIsmithmodelstandardinverse <- glm(Transform ~ LAI, family = inverse.gaussian(link="1/mu^2"), data = Zos_Smith_Metrics)
summary(LAIsmithmodelstandardinverse)
##Leaf area index does not effectively explain model
##better fit with site
##With just interaction term
AIC(sitesmithmodelstandardgamma, k = 12)
AIC(sitesmithmodelstandardgaus, k = 12)
AIC(sitesmithmodelstandardinverse, k = 12)
AIC(LAIsmithmodelstandardinverse)
AIC(LAIsmithmodelstandardgamma)
AIC(LAIsmithmodelstandardgaus)
AIC(smithmodelstandardinverse, k = 25)
AIC(smithmodelstandardgamma, k = 25)
AIC(smithmodelstandardgaus, k = 25)

LAIsmith <- ggplot(Zos_Smith_Metrics)
LAIsmith <- LAIsmith + geom_point(aes(LAI,smithora)) + xlim(0,1500) + facet_wrap(~site)
print(LAIsmith) 

##based on this graph I do not understand how thier could no effect of site on the relationship between smithora and LAI
##On some sites smithora cover increases linearly with space for colonization and at other sites it does not
###interior sites seem to not reach a higher level deeper sites?? 
##Ok for my purposes just show that Smithora is varying in choked pass ALOT and we realy don't know why

citation("phyloseq")

###AIC model selection to explain smithora abundance
Cand.mod <- list()

Cand.mod[[1]] <- glm(Transform ~ 1 , family = inverse.gaussian(link="1/mu^2"), data = Zos_Smith_Metrics)
##model seatable and control and urchin
Cand.mod[[2]] <- glm(Transform ~ site , family = inverse.gaussian(link="1/mu^2"), data = Zos_Smith_Metrics)
##model control,urchin and fed
Cand.mod[[3]] <- glm(Transform ~ LAI, family = inverse.gaussian(link="1/mu^2"), data = Zos_Smith_Metrics)
##model control,urchin and temp
Cand.mod[[4]] <- glm(Transform ~ site*LAI, family = inverse.gaussian(link="1/mu^2"), data = Zos_Smith_Metrics)

Modnames <- c("null", "site",
              "LAI", "LAI + site")
library(AICcmodavg)
##model selection table based on AICc

aictab(cand.set = Cand.mod, modnames = Modnames)
evidence(aictab(cand.set = Cand.mod, modnames = Modnames))
###so smihora variation is best explained by location in the meadow not by the area for colonization on the shoot



