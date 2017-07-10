library(dplyr)
library(stats)
library(base)
getwd()
setwd('C:/Users/Gwen/Desktop/Rpractice')

##species.data is with amphipods species.data2 is without because they messed
##up the diversity calculation

species.data <- read.csv('grazerssansamphi.csv')
Grazers <- read.csv('grazersformatted.csv')
head(Grazers)
head(species.data)
Grazers
species.data2 <- Grazers %>% 
  select(-X)%>%
  slice(1:6)

species.data2
library(base)

citation(package = "vegan", lib.loc = NULL)

Grazers2 <- read.csv('grazerstransformed.csv')

Grazers2
 

species.data
names(species.data)
site.data <- read.csv('grazersite.csv')
library(vegan)
library(broom)
getwd()

names(site.data)

head(site.data)

library(dplyr)

site.data
site.data2 <- site.data %>% 
  select(- contains("X")) %>%
  slice(1:6)

head(site.data2)

site.data2$location2 <- as.factor(c("Inner","Edge","Inner","Edge","Edge","Inner"))

  

site.data2[ , -which(names(site.data2) %in% c("X.1","X.2","X.3","X.4"))]

View(site.data)
head(species.data)
View(Grazers)
##diversity calculation 
diversity(species.data, index = "shannon")
diversity(species.data, index = "simpson")
fisher.alpha(species.data)

site.data2$fisher<-fisher.alpha(species.data)

head(site.data2)

locationdiversity<-lm(fisher~Location, data =site.data2) 

anova(locationdiversity)
par(mfrow=c(2,2))
plot(locationdiversity)

dev.off()

shapiro.test(site.data2$fisher)

##eveness
diversity(species.data, index = "shannon")/log(specnumber(species.data))
head(site.data2)
##rarefaction
min(rowSums(species.data))
rarefy(species.data, 62)
site.data2$rarefy<-(rarefy(species.data, 62))
shapiro.test(site.data2$rarefy)
locationrarefy <- lm(rarefy~Location, data =site.data2)
summary(locationrarefy)
par(mfrow=c(2,2))
plot(locationrarefy)
anova(locationrarefy)
View(species.data)

##abundances over all sites
par(mfrow=c(1,1))
plot(fisherfit(colSums(species.data)))

plot(rad.lognormal(colSums(species.data))) #I first summed abundances across all sites with colSums, then fit a log normal, then plotted it

Int.data<-species.data[site.data$Location %in% "Interior",]
Edge.data<-species.data[site.data$Location %in% "Edge",]
#I made mini data frames of edge and int ^ to see abundances within sites
plot(rad.lognormal(colSums(Int.data)))
plot(fisherfit(colSums(Int.data)))
plot(rad.lognormal(colSums(Edge.data)))
plot(fisherfit(colSums(Edge.data)))

##look at a variety of models to see which one best fits my data
##lowest AIC is the best so Preemption
radlattice(radfit(colSums(species.data)))

str(site.data)
library(BiodiversityR)

RankAbun.1 <- rankabundance(species.data)
RankAbun.1 # a dataframe of the rank of each species, needed for next command
rankabunplot(RankAbun.1,scale='abundance', addit=FALSE, specnames=c(1,2,3)) #rank abudnance plot, labelling the most common 3 species
#or you can overlay the rank-abundance lines for each level of a factor as follows:
par(mfrow=c(1,1))
head(species.data)
head(site.data2)
species.data
Grazers
rankabuncomp(species.data, y=site.data2, factor='Location',scale='proportion', legend=TRUE)
class(site.data2)
library(vegan)
Grazers
head(Grazers2)

dissim.mat <- vegdist(species.data2, method="bray")
dissim.mat

?vegdist #opens up the documentation that explains all other dissimilarity metrics you could use in place of "bray"
#these include: "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".

fit <- hclust(dissim.mat, method="average") #uses the average linkage method described in class
plot(fit)
dev.off()
rect.hclust(fit, k=4, border="red") # if you wanted to emphasize 4 clusters, draw red boxes around them

plot(fit);
rect.hclust(fit, h=0.5, border="red") #if you wanted to emphasize clusters <0.5 different

myNMDS<-metaMDS(species.data2,k=2)
myNMDS #most important: is the stress low?
stressplot(myNMDS) #low stress means that the observed dissimilarity between site pairs matches that on the 2-D plot fairly well (points hug the line)

plot(myNMDS)#sites are open circles and species are red +'s ...but it might be nice to label these.

#the following commands create layers on a plot, and should be run sequentially
ordiplot(myNMDS,type="n") #this clears the symbols from the plot
orditorp(myNMDS,display="species",col="red",air=0.01) #this adds red species names
orditorp(myNMDS,display="sites",cex=0.75,air=0.01) #this adds black site labels, cex is the font size

site.mat<-vegdist(site.data2$drift.seaweed.weight, method="euclidean") #A1 is how thick the A1 soil horizon was, and is continuous
site.mat

head(site.data2)

site.mat<-vegdist(site.data2$macroepiphyte.weight, method="euclidean") #A1 is how thick the A1 soil horizon was, and is continuous
site.mat


mantel(dissim.mat, site.mat, method="pearson", permutations=9999)

head(site.data2)

library(ggplot2)


ordihull(myNMDS,groups=site.data2$Location,draw="polygon", col=1:4, lwd=3)
orditorp(myNMDS,display="sites", labels = as.character(site.data2$Site), cex=0.75,air=0.1)

ord.fit <- envfit(myNMDS ~ location2 + macroepiphyte.weight, data=site.data2, perm=999)
ord.fit
text(myNMDS, display="sites", labels = as.character(Location))
orditorp(myNMDS,display="species",col="red",air=0.01)
orditorp(myNMDS,display="sites", labels = as.character(site.data2$Site), cex=0.75,air=0.01)

ordiplot(myNMDS,type="n") 
ordispider(myNMDS,groups=site.data$Location,spiders="centroid",col="black",lab=F)#this is the new line
orditorp(myNMDS,display="species",col="red",air=0.01) 
orditorp(myNMDS,display="sites",cex=0.75,air=0.01)

#your friend here is adonis

adonis( dissim.mat ~ Location, data=site.data2, permutations=9999)



##looking at amphipods my gammarid amphipods differed strongly between sites
amphipod <- read.csv('amphipods.csv')


View(amphipod)
amphipod
shapiro.test(amphipod$Amphipod.Abundance)
qqnorm(amphipod$Amphipod.Abundance)
qqline(amphipod$Amphipod.Abundance)
library(ggplot2)
library(broom)''
amphi <- ggplot(amphipod)
comparisonsamphipod <- ggplot(amphipod, aes(x = Location, y = Amphipod.Abundance)) + geom_boxplot(size = 1) + theme_bw() + theme(panel.grid = element_blank()) + ylab(expression("Abundance")) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(legend.position = "none") + theme(strip.text.y = element_text(face = 4)) + xlab(expression("Location")) + theme(axis.title.x = element_text(vjust = -1, size =14, color = "black")) + theme(axis.title.y = element_text(vjust = -1, size = 14, color = "black")) + theme(axis.text.x = element_text(color = "black", vjust = -1, size = 12)) + theme(axis.text.y = element_text(color = "black", size = 8))

print(comparisonsamphipod)


aggregate(amphipod, list(Location = "Interior"), mean)
amphipodsite <- lm(Amphipod.Abundance~Location,data = amphipod)
summary(amphipodsite)

library(ggplot2)

##Location had a significant effect on amphipod abundance 
anova(amphipodsite)

library(ggplot2)

STDEV <- aggregate(amphipod$Amphipod.Abundance, list(Location = amphipod$Location), sd)
STDEV
std <- function(x) sd(x)/sqrt(length(x))
StandardError <- aggregate(amphipod$Amphipod.Abundance, list(Location = amphipod$Location), std)
StandardError

citation(package = "vegan", lib.loc = NULL)


means <- as.numeric(c(66,511))
Location <- as.factor(c("Interior","Edge"))
amphicompare <- data.frame(Location,means,STDEV,StandardError)
upperlimit<-as.numeric(c((66+24.7),(511+37.12)))
lowerlimit<-as.numeric(c((66-24.7),(511-37.12)))
amphicompare2 <- data.frame(Location,means,upperlimit,lowerlimit)
amphicompare2




ggplot(amphicompare2, aes(y=amphicompare2$means, x=amphicompare2$Location)) + geom_bar(stat= "identity", position= "dodge") + labs(x="Location") + labs(y="Abundance")+ theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white")) + geom_errorbar(aes(ymin=amphicompare2$LowerLimit, ymax=amphicompare2$UpperLimit), width=.2,position=position_dodge(.9))
 
ggplot(amphicompare2, aes(y=amphicompare2$means, x=amphicompare2$Location)) + geom_bar(stat= "identity", position= "dodge") 

ggplot(amphicompare2, aes(y=amphicompare2$means, x=amphicompare2$Location)) + geom_bar(stat= "identity", position= "dodge") + labs(x="Location") + labs(y="Gammaridian amphipod abundance")+ theme_bw() + theme(panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_line(colour = "white")) + geom_errorbar(aes(ymin=amphicompare2$lowerlimit, ymax=amphicompare2$upperlimit),  width=.2,position=position_dodge(.9))


