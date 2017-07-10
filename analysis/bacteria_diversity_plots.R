

##name the new data set by what it is organizing the taxa by
##find out how mny colours I will need
set.seed(3)
Familycolour <- c("#c7a2ab","#2bf04a","#403fe8","#76e70f","#7f44f2","#00f483","#b738eb","#82ff8c","#0017a6","#fff64b","#0165f5","#d2be00","#b861ff","#5aa400","#d600b5","#02df83","#740084","#00a944","#f584ff","#1c7c00","#6e7bff","#83a500","#00268b","#fdff8a","#000e69","#c1ffa6","#3f0063","#65ffc2","#ff3b26","#00f7fd","#9e0002","#02d8c4","#ff4150","#99fff2","#a4001c","#8bebf","#be003d","#bcffcb","#29003d","#fdffd1","#0c0034","#ffc571","#00225f","#b18400","#cb88ff","#3f7200","#e4267c","#00613f","#ff3e7c","#003c2d","#ff546f","#37abff","#cb6900","#0266be","#a66500","#7eafff","#942c00","#a7a4ff","#716a00","#730e5b","#f3fffb","#5a0005","#c0eaff","#b90051","#00747e","#ff876c","#001b2b","#ffa66a","#002f55","#ffc197","#001813","#f7d8ff","#1b2300","#ff7ba3","#2f4100","#ffa7ab","#003f52","#743200","#0084b0","#5a1a00","#8d3f5a","#4e3500","#66002e","#290a00")

Gwensmith_longfamily_experiments2$before_after_transplant <- factor(Gwensmith_longfamily_experiments2$before_after_transplant,levels = c("Before detachment and transplant","After detachment and transplant"))
levels(Gwensmith_longfamily_experiments2$before_after_transplant)

tryingbarexp <- ggplot(Gwensmith_longfamily_experiments2)

tryingbarexp <- ggplot(Gwensmith_longfamily_experiments2, aes(x = blade_ID, y = Abundance, fill = Family)) + geom_bar(stat = "identity") + scale_fill_manual(values = Familycolour) + facet_wrap(~before_after_transplant) + ylab(expression("Relative Abundance")) + ggtitle("Relative abundance of bacterial species by Family") + theme(panel.grid = element_blank()) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Relative Abundance")) + theme(axis.text.x = element_text(angle = 60, size = 12, colour = "black",hjust = 1), axis.title.x = element_blank())

##open up new window
par(mfrow = c(1,2))
windows(width=10, height=6)
print(tryingbarexp)

###try with just controls

Gwensmith_longfamily_controls2$before_after_transplant <- factor(Gwensmith_longfamily_controls2$before_after_transplant,levels = c("Before detachment","After detachment"))

tryingbarcont <- ggplot(Gwensmith_longfamily_controls2)

tryingbarcont <- ggplot(Gwensmith_longfamily_controls2, aes(x = blade_ID, y = Abundance, fill = Family)) + geom_bar(stat = "identity") + scale_fill_manual(values = Familycolour) + facet_grid(~ before_after_transplant) + ylab(expression("Relative Abundance")) + ggtitle("Relative abundance of bacterial species by Family") + theme(panel.grid = element_blank()) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Relative Abundance")) + theme(axis.text.x = element_text(angle = 60, size = 12, colour = "black",hjust = 1), axis.title.x = element_blank())

##open up new window
windows(width=10, height=6)
print(tryingbarcont)

tryingbaramb <- ggplot(Gwensmith_longfamily_ambient2)

tryingbaramb <- ggplot(Gwensmith_longfamily_ambient2, aes(x = blade_ID, y = Abundance, fill = Family)) + geom_bar(stat = "identity") + scale_fill_manual(values = Familycolour) + ylab(expression("Relative Abundance")) + ggtitle("Relative abundance of bacterial species by Family") + theme(panel.grid = element_blank()) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Relative Abundance")) + theme(axis.text.x = element_text(angle = 60, size = 12, colour = "black",hjust = 1), axis.title.x = element_blank())

##open up new window
windows(width=10, height=6)
print(tryingbaramb)

tryingbarsmith <- ggplot(Gwensmith_longfamily_smith2)

tryingbarsmith <- ggplot(Gwensmith_longfamily_smith2, aes(x = blade_ID, y = Abundance, fill = Family)) + geom_bar(stat = "identity") + scale_fill_manual(values = Familycolour) + ylab(expression("Relative Abundance")) + ggtitle("Relative abundance of bacterial species by Family") + theme(panel.grid = element_blank()) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Relative Abundance")) + theme(axis.text.x = element_text(angle = 60, size = 12, colour = "black",hjust = 1), axis.title.x = element_blank())

##open up new window
windows(width=10, height=6)
print(tryingbarsmith)

##Yay plots done!!





##Smithora blades have quite low diversity at the class level which is interesting. Lots of gammaproteobactera, and flavobacteria




###Diversity stuff, not reliable so will not use but if anyone asks you still wrote the code to fdeal with it
##To look at diversity I found a blog that subsamples at the minimum reads 100 times. Then takes the average of these subsamples. It seems like a good technique. Its like rarefying a bunch. Spreading out your error
###plots look interesting lets look at diversity
min_lib <- min(sample_sums(Gwensmith))
# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(Gwensmith)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(Gwensmith)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(Gwensmith)

# It is always important to set a seed when you subsample so your result is replicable 
set.seed(3)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(Gwensmith, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
}

Sample <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(Sample, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
Sample <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(Sample, mean, sd, measure)

##Making a big data frame with all the stuff
alpha <- rbind(rich_stats, even_stats)
s <- data.frame(sample_data(map))
alphadiv <- merge(alpha, s, by = "Sample")
View(alphadiv)

alphadiv_long <- melt(alphadiv, id.vars = "blade_ID")
alphadiv_diversity <- filter(alphadiv, measure == "Inverse Simpson")
View(alphadiv_diversity)
alphadiv_diversity_long <- melt(alphadiv_diversity, id=c("blade_ID","before_after_transplant", "Smithora_status","type_swab"), measure.vars = c("mean","sd"))

alphadiv_justexp <- filter(alphadiv_diversity_long, before_after_transplant != "ambient" & variable == "mean")

alphadiv_after <- filter(alphadiv_justexp, before_after_transplant == "after")
alphadiv_before <- filter(alphadiv_justexp, before_after_transplant == "before")

alphadiv_before$divbefore <- c(alphadiv_before$value)
alphadiv_after$divafter <- c(alphadiv_after$value)

total_div <- full_join(alphadiv_before, alphadiv_after, by = "blade_ID")
View(total_div)
dev.off()
hist(total_div$value.x)
hist(total_div$value.y)


total_div_exp <- filter(total_div, type_swab.x == "experiment")

###A paired t.test showed that transplant had no effect on bacterial diversity. However when comparsing samples before transplant (shoots recently removed from meadow, and ambient shoots smithora seems to be correlated with lower bacterial diversity)

smithora_div_before <- ggplot(total_div)
smithora_div_before <- smithora_div_before + geom_boxplot(aes(x = Smithora_status.x, y = divbefore))
print(smithora_div_before)
smithora_div_after <- ggplot(total_div)
smithora_div_after <- smithora_div_after + geom_boxplot(aes(x = Smithora_status.y, y = divafter))
print(smithora_div_after)

##No effect of transplant on bacterial diversity 
t.test(total_div$divbefore,total_div$divafter, paired = TRUE)

boxplot(total_div$divbefore,total_div$divafter)

smithoramodel<- lm(total_div$divbefore ~ total_div$Smithora_status.x)
boxplot(total_div$divbefore,total_div$Smithora_status.x)
boxplot(total_div$divafter,total_div$Smithora_status.y)
##This is whta phyloseq produces
boxplots <- plot_richness(Gwensmith, x = "before_after_transplant", color = "Smithora_status") + geom_boxplot()
boxplots <- boxplots + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Diversity Value")) + xlab(expression("Type of Sample")) + ggtitle("Diversity Values from Various Measures") + theme(plot.title = element_text(vjust= 2, size = 14)) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3)) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
windows(width=10, height=6) 
print(boxplots)

t.test(alphadiv$mean,y2,paired=TRUE)

trygroup <- %>%
  group_by(alphadiv,blade_ID) %>%
  View()

alphadiv_long <- melt(alphadiv, idvar = blade_ID)
library(nlme)

dev.off()
hist(alphadiv_diversity_seagrass$mean)

alphadiv_diversity_seagrass <- filter(alphadiv_diversity, before_after_transplant != "ambient")
model_div_linear <- lm(mean ~ before_after_transplant*Smithora_status, data = alphadiv_diversity_seagrass)
AIC(model_div_linear)
summary(model_div_linear)
library(AICcmodavg)
model_div_gamma <- glm(mean ~ before_after_transplant*Smithora_status, family = Gamma, data = alphadiv_diversity_seagrass)
AIC(model_div_gamma)
model_div_inverse <- glm(mean ~ before_after_transplant*Smithora_status, family = inverse.gaussian, data = alphadiv_diversity_seagrass)
AIC(model_div_inverse)



Cand.mod <- list()

Cand.mod[[1]] <- glm(mean ~ 1, family = Gamma, data = alphadiv_diversity)

Cand.mod[[2]] <- glm(mean ~ before_after_transplant, family = Gamma, data = alphadiv_diversity_seagrass)

Cand.mod[[3]] <- glm(mean ~ Smithora_status, family = Gamma, data = alphadiv_diversity)

Cand.mod[[4]] <- glm(mean ~ Smithora_status*before_after_transplant, family = Gamma, data = alphadiv_diversity)

View(alphadiv_diversity)

alphadiv_diversity$Community <- c("Ambient Zostera Edge","Ambient Zostera Edge","Ambient Zostera Interior","Ambient Zostera Interior","After Detachment Edge","After Detachment and Transplant to Interior","After Detachment Edge","After Detachment Edge","After Detachment and Transplant to Interior","After Detachment and Transplant to Interior","After Detachment and Transplant to Edge","After Detachment and Transplant to Edge","After Detachment Interior","After Detachment Interior","After Detachment and Transplant to Edge","After Detachment Interior","Before Detachment Edge","Before Detachment and Transplant to Interior","Before Detachment Edge","Before Detachment and Transplant to Interior","Before Detachment and Transplant to Interior","Before Detachment and Transplant to Edge","Before Detachment and Transplant to Edge","Before Detachment Interior","Before Detachment Interior","Before Detachment and Transplant to Interior","Before Detachment Interior","Ambient Smithora Blade","Ambient Smithora Blade","Ambient Smithora Blade","Ambient Smithora Blade")

dev.off()
hist(alphadiv_diversity_long$mean)



Modnames <- c("null", "transplant",
              "smith","smithTrans")
library(AICcmodavg)
##model selection table based on AICc

aictab(cand.set = Cand.mod, modnames = Modnames)
evidence(aictab(cand.set = Cand.mod, modnames = Modnames))

View(alphadiv_diversity_seagrass)
alphadiversity_controls2$before_after_transplant <- factor(alphadiversity_controls2$before_after_transplant,levels = c("Before detachment","After detachment (1 month)"))

str(alphadiv_diversity$Smithora_status)
levels(alphadiv_diversity$Smithora_status)
alphadiv_diversity$Smithora_status <- factor(alphadiv_diversity$Smithora_status,levels = c("Present","Absent","smithora","unknown"))

alphadiversity_justexp <- filter(alphadiv_diversity, before_after_transplant != "ambient")
alphadiversity_justexp <- filter(alphadiv_diversity, blade_ID != "edge_d")
alphadiversity_controls <- filter(alphadiversity_justexp, type_swab == "control")
alphadiversity_experiments <- filter(alphadiv_diversity, type_swab == "experiment")
View(alphadiversity_experiments)
alphadiversity_ambient <- filter(alphadiv_diversity, before_after_transplant== "ambient" & host_species == "Zostera")
alphadiversity_smith <- filter(alphadiv_diversity, host_species == "Smithora")

alphadiversity_smith2 <- alphadiversity_smith
alphadiversity_smith2$blade_ID <- as.character(alphadiversity_smith2$blade_ID)
alphadiversity_smith2$blade_ID[alphadiversity_smith$blade_ID == "random_blade_1"] <- "Smithora Thallus 1"
alphadiversity_smith2$blade_ID[alphadiversity_smith$blade_ID == "random_blade_2"] <- "Smithora Thallus 2"
alphadiversity_smith2$blade_ID[alphadiversity_smith$blade_ID == "random_blade_3"] <- "Smithora Thallus 3"
alphadiversity_smith2$blade_ID[alphadiversity_smith$blade_ID == "random_blade_4"] <- "Smithora Thallus 4"

alphadiversity_ambient2 <- alphadiversity_ambient
alphadiversity_ambient2$blade_ID <- as.character(alphadiversity_ambient2$blade_ID)
alphadiversity_ambient2$blade_ID[alphadiversity_ambient$blade_ID == "ambient_edge_1"] <- "Edge Zostera 1"
alphadiversity_ambient2$blade_ID[alphadiversity_ambient$blade_ID == "ambient_edge_2"] <- "Edge Zostera 2"
alphadiversity_ambient2$blade_ID[alphadiversity_ambient$blade_ID == "ambient_inner_1"] <- "Interior Zostera 1"
alphadiversity_ambient2$blade_ID[alphadiversity_ambient$blade_ID == "ambient_inner_2"] <- "Interior Zostera 2"

alphadiversity_controls2 <- alphadiversity_controls

alphadiversity_controls2$before_after_transplant <- factor(alphadiversity_controls2$before_after_transplant,levels = c("Before detachment","After detachment (1 month)"))

alphadiversity_controls2$before_after_transplant <- as.character(alphadiversity_controls2$before_after_transplant)
alphadiversity_controls2$before_after_transplant[alphadiversity_controls$before_after_transplant == "after"] <- "After detachment (1 month)"
alphadiversity_controls2$before_after_transplant[alphadiversity_controls$before_after_transplant == "before"] <- "Before detachment"

alphadiversity_controls2$blade_ID <- as.character(alphadiversity_controls2$blade_ID)
alphadiversity_controls2$blade_ID[alphadiversity_controls$blade_ID == "edge_b"] <- "Control edge B"
alphadiversity_controls2$blade_ID[alphadiversity_controls$blade_ID == "edge_f"] <- "Control edge F"
alphadiversity_controls2$blade_ID[alphadiversity_controls$blade_ID == "inner_f"] <- "Control inner F"
alphadiversity_controls2$blade_ID[alphadiversity_controls$blade_ID == "inner_g"] <- "Control inner G"
alphadiversity_controls2$blade_ID[alphadiversity_controls$blade_ID == "inner_k"] <- "Control inner K"

View(alphadiversity_controls2)

View(alphadiversity_experiments2)

alphadiversity_experiments$before_after_transplant <- factor(alphadiversity_experiments$before_after_transplant,levels = c("before","after"))

alphadiversity_experiments2 <- alphadiversity_experiments

alphadiversity_experiments2$before_after_transplant <- as.character(alphadiversity_experiments2$before_after_transplant)
alphadiversity_experiments2$before_after_transplant[alphadiversity_experiments$before_after_transplant == "after"] <- "Detached & transplanted 1 month "
alphadiversity_experiments2$before_after_transplant[alphadiversity_experiments$before_after_transplant == "before"] <- "Before detachment and transplant"


alphadiversity_experiments2$blade_ID <- as.character(alphadiversity_experiments2$blade_ID)
alphadiversity_experiments2$blade_ID[alphadiversity_experiments$blade_ID == "edge_c"] <- "Experiment edge C"
alphadiversity_experiments2$blade_ID[alphadiversity_experiments$blade_ID == "edge_h"] <- "Experiment edge H"
alphadiversity_experiments2$blade_ID[alphadiversity_experiments$blade_ID == "edge_k"] <- "Experiment edge K"
alphadiversity_experiments2$blade_ID[alphadiversity_experiments$blade_ID == "inner_a"] <- "Experiment inner A"
alphadiversity_experiments2$blade_ID[alphadiversity_experiments$blade_ID == "inner_d"] <- "Experiment inner D"
alphadiversity_experiments2$blade_ID[alphadiversity_experiments$blade_ID == "inner_h"] <- "Experiment inner H"

alphadiversity_experiments2$before_after_transplant <- factor(alphadiversity_experiments2$before_after_transplant,levels = c("Before detachment and transplant","Detached & transplanted 1 month"))

View(alphadiversity_experiments2)
dev.off()
windows(width=10, height=6)
par(mfrow=c(2,2))
par(mfcol=c(2,2))

diversityboxplot <- ggplot(alphadiversity_experiments2, aes(x = blade_ID, y = mean, color = Smithora_status)) + geom_point(size = 1) + geom_errorbar(aes(ymin = mean + sd, ymax = mean - sd)) + facet_wrap(~before_after_transplant) + theme_bw() + scale_fill_manual(values = Familycolour) + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Diversity Value")) + xlab(expression("Type of Sample")) + ggtitle("Diversity Values from Random Sampling Method") + theme(panel.grid = element_blank()) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + theme(axis.text.x = element_text(angle = 60, size = 12, colour = "black",hjust = 1), axis.title.x = element_blank())

print(diversityboxplot)
#inner conrols had a decrease in divesity after transplant likely cause they were not colonized by smithora everything else increased in diversity following transplant and temporal changes

controlsdiversityboxplot <- ggplot(alphadiversity_controls2, aes(x = blade_ID, y = mean, color = Smithora_status)) + geom_boxplot(size = 1) + facet_wrap(~before_after_transplant) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Diversity Value")) + xlab(expression("Type of Sample")) + ggtitle("Diversity Values from Random Sampling Method") + theme(panel.grid = element_blank()) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + theme(axis.text.x = element_text(angle = 60, size = 12, colour = "black",hjust = 1), axis.title.x = element_blank())
print(controlsdiversityboxplot)

ambientdiversityboxplot <- ggplot(alphadiversity_ambient2, aes(x = blade_ID, y = mean, color = Smithora_status)) + geom_boxplot(size = 1) + facet_wrap(~before_after_transplant) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Diversity Value")) + xlab(expression("Type of Sample")) + ggtitle("Diversity Values from Random Sampling Method") + theme(panel.grid = element_blank()) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + theme(axis.text.x = element_text(angle = 60, size = 12, colour = "black",hjust = 1), axis.title.x = element_blank())
print(ambientdiversityboxplot)

smithoradiversityboxplot <- ggplot(alphadiversity_smith2, aes(x = blade_ID, y = mean, color = Smithora_status)) + geom_boxplot(size = 1) + facet_wrap(~before_after_transplant) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Diversity Value")) + xlab(expression("Type of Sample")) + ggtitle("Diversity Values from Random Sampling Method") + theme(panel.grid = element_blank()) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + theme(axis.text.x = element_text(angle = 60, size = 12, colour = "black",hjust = 1), axis.title.x = element_blank())
print(smithoradiversityboxplot)


