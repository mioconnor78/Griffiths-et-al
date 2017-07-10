source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library("phyloseq")
packageVersion("phyloseq")

library("ggplot2")
packageVersion("ggplot2")

library("scales")
packageVersion("scales")

library("grid")
packageVersion("grid")

theme_set(theme_bw())
getwd()
setwd("C:/Users/Gwen/Desktop/Rpractice")
gwen_griff_biom = system.file("16s_seagrass_OTU_Table.filtered", package = "phyloseq")
gwen_griff = import_biom(gwen_griff_biom)
print(gwen_griff)

vignette("phyloseq_basics")

`?`("phyloseq-package")
`?`(phyloseq)
?system.file
##load packages
library(dplyr)
library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)

##import biom files
getwd()
setwd("C:/Users/Gwen/Desktop/Rpractice/gwen_griffiths")

?system.file
##importing Biom file
##C:/Users/Gwen/Desktop/Rpractice/gwen_griffiths is the file path
Gwen_Biom = file.path("C:/Users/Gwen/Desktop/Rpractice/gwen_griffiths", "16s_seagrass_OTU_Table.filtered.wtaxa.biom")
##you can also add a phylogenetic tree so you can get unifrac distances. Unifrac is a distance measure that incorporates phylogenetic relationships
trees = file.path("C:/Users/Gwen/Desktop/Rpractice/gwen_griffiths","16s_makephylo_fasttree.tre")
##Physeq imports both simultaneously and then names them as one object
physeq <- import_biom(Gwen_Biom,trees)
##will produce list of what your object contains. OTU table, taxa table, and a phylogenetic tree!
physeq

##adding sample map file. This will give you sampling information about your data. Smithora, before or after transplant etc. 
mapfile <- read.csv("map.csv")
##important to designate your map as sample_data
map <- sample_data(mapfile)

# Use merge to add your map to your sample file. Assign rownames to be Samples in your OTU table and taxonomy table
rownames(map) <- map$Sample
##merge your map to your OTU table and phylogenetic tree, make sure the number of samples matches
phylo_merge <- merge_phyloseq(physeq, map)
phylo_merge

##Now I want my taxonomy classes to match my reads for each grouping
colnames(tax_table(phylo_merge))
colnames(tax_table(phylo_merge)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")

head(map)
##Now lets subsample so I am just looking at my samples and lets remove mitochondria and chloroplasts, and any taxa that are 0. 
phylo_sub <- phylo_merge %>%
  subset_samples(survey_type == "Gwen_tansplant") %>%
  prune_taxa(taxa_sums(.) > 0, .)
Gwensmith <- phylo_sub %>%
  subset_taxa(
    Kingdom == "k__Bacteria" &
      Family  != "f__mitochondria" &
      Class   != "c__Chloroplast"
  )

##check subset to make sure it as all my samples (31)
Gwensmith

##This makes a sample summary data frame so you can see your distribution of reads and such
sample_sum_df <- data.frame(sum = sample_sums(Gwensmith))
library(ggplot2)
# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
##I kind of have a weird distribution, is this problematic
##mean, max, and min of read counts
smin <- min(sample_sums(Gwensmith))
smean <- mean(sample_sums(Gwensmith))
smax <- max(sample_sums(Gwensmith))
smin
smean
smax
##I will use the minimum read count later
###I made multiple plots to try and show the composition of the the samples. What bacteria are actually present
##make phylum graph. Grouped by hylum first get rid of taxa at less than 0.02 abundance
Gwensmith_longPhylum <- Gwensmith %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum) 
##name the new data set by what it is organizing the taxa by
##find out how mny colours I will need
length(levels(Gwensmith_longPhylum$Phylum))
##16
##make a vector with discrete colours to map the phylum colour palette 
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

# Plot. Group samples by smithora status to see overall what bacteria are present with different levels of smithora, before/after the transplant and as ambient samples. Please note all labelled as ambient are from august.Compare to inner/edge before to seee differences. Also because I have unequal sample size and unequal reads the bargraphs are not the same height. All abundances are relative abundances compared to what also appeared under the samples in each group. 
tryingbar <- ggplot(Gwensmith_longPhylum, aes(x = Smithora_status, y = Abundance, fill = Phylum)) + geom_bar(stat = "identity") + scale_fill_manual(values = phylum_colors) + facet_wrap(~before_after_transplant) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
tryingbar <- tryingbar + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Relative Abundance")) + xlab(expression("Smithora Status")) + ggtitle("Relative abundance of bacterial species by Phylum ") + theme(plot.title = element_text(vjust= 2, size = 14)) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3))

##open up new window
windows(width=10, height=6)
print(tryingbar)
##At phylum level there doesn't seem to be a large change in time. But you can see some striking differences between the way communties changed. Shoots that had smithora after the treatment looked very diffferent from those that didnt

##make Class graph. See if the similarities hold. Or some more patterns start to emeerge. 
Gwensmith_longClass <- Gwensmith %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at Class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class) 
##Same steps as the phylum graph
##how many colours do we need
length(levels(Gwensmith_longClass$Class))
##make a vector with 30 colours
Classcolour <- c("#00114d",
                 "#dfaa00",
                 "#0157ff",
                 "#b1ad00",
                 "#0149c5",
                 "#009d30",
                 "#a60095",
                 "#5b8800",
                 "#ca7cff",
                 "#ffd869",
                 "#000065",
                 "#d0ffaa",
                 "#ff75f0",
                 "#008754",
                 "#ff438f",
                 "#0197a9",
                 "#d4003d",
                 "#028ad5",
                 "#963100",
                 "#7986ff",
                 "#ff8a54",
                 "#a8a6ff",
                 "#5d2200",
                 "#ade9ff",
                 "#aa002d",
                 "#fdffe3",
                 "#5a0011",
                 "#ffe9fe",
                 "#004d56",
                 "#ffa3eb")
                    
# Plot 
tryingbarclass <- ggplot(Gwensmith_longClass, aes(x = Smithora_status, y = Abundance,fill = Class)) + geom_bar(stat = "identity") + scale_fill_manual(values=Classcolour) + facet_wrap(~before_after_transplant) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
tryingbarclass <- tryingbarclass + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Relative Abundance")) + xlab(expression("Smithora Status")) + ggtitle("Relative abundance of bacterial species by Class ") + theme(plot.title = element_text(vjust= 2, size = 14)) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3))

windows(width=10, height=6)
print(tryingbarclass)
##Smithora blades have quite low diversity at the class level which is interesting. Lots of gammaproteobactera, and flavobacteria

##Make Family Graph
Gwensmith_longfamily <- Gwensmith %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Family) 
##trying family
length(levels(Gwensmith_longfamily$Family))

Familycolour <- c("#c7a2ab","#2bf04a","#403fe8","#76e70f","#7f44f2","#00f483","#b738eb","#82ff8c","#0017a6","#fff64b","#0165f5","#d2be00","#b861ff","#5aa400","#d600b5","#02df83","#740084","#00a944","#f584ff","#1c7c00","#6e7bff","#83a500","#00268b","#fdff8a","#000e69","#c1ffa6","#3f0063","#65ffc2","#ff3b26","#00f7fd","#9e0002","#02d8c4","#ff4150","#99fff2","#a4001c","#8bebf","#be003d","#bcffcb","#29003d","#fdffd1","#0c0034","#ffc571","#00225f","#b18400","#cb88ff","#3f7200","#e4267c","#00613f","#ff3e7c","#003c2d","#ff546f","#37abff","#cb6900","#0266be","#a66500","#7eafff","#942c00","#a7a4ff","#716a00","#730e5b","#f3fffb","#5a0005","#c0eaff","#b90051","#00747e","#ff876c","#001b2b","#ffa66a","#002f55","#ffc197","#001813","#f7d8ff","#1b2300","#ff7ba3","#2f4100","#ffa7ab","#003f52","#743200","#0084b0","#5a1a00","#8d3f5a","#4e3500","#66002e","#290a00")
                  

tryingbarfamily <- ggplot(Gwensmith_longfamily, aes(x = Smithora_status, y = Abundance, fill = Family)) + geom_bar(stat = "identity") + scale_fill_manual(values=Familycolour) + facet_wrap(~before_after_transplant) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
tryingbarfamily <- tryingbarfamily + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Relative Abundance")) + xlab(expression("Smithora Status")) + ggtitle("Relative abundance of bacterial species by Family ") + theme(plot.title = element_text(vjust= 2, size = 14)) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3))
##The before absent seem to look more like the after present. This is because blades that were colonized likely had their unique bactera community that they maintained after being transplanted. This idea could be really interesting, but I don't know how to develop it! There are also some unique families that only smithora blades, and shoots that have smithora share.  
windows(width=10, height=6)
print(tryingbarfamily)

##Now I am going to try two different ordinations using different distance metrics including unifrac. 
###ordinations

library(scales)
library(vegan)
library(grid)

##scale reads assuming I am scaling them to the minimum number of reads "n". I know library size can also mean the length of the sequences. Am I being a dumbass here?
##define a function that
scale_reads <- function(physeq, n) {
  physeq.scale <-
    transform_sample_counts(physeq, function(x) {
      (n * x/sum(x))
    })
  otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}
##scale reads to minimum library size which I think is just the total number of reads
Gwensmith_scale <- scale_reads(Gwensmith,7244)

# Fix smithora levels in sample_data
sample_data(Gwensmith_scale)$Smithora_status <- factor(
  sample_data(Gwensmith_scale)$Smithora_status, 
  levels = c("absent","present","smithora")
)
sample_data(Gwensmith_scale)$before_after_transplant <- factor(
  sample_data(Gwensmith_scale)$before_after_transplant, 
  levels = c("before","after","ambient")
)

##check that my physeq file is still in order
Gwensmith_scale

# Ordinate look at an uncontrained PCA with unifrac
Gwensmith_pcoa <- ordinate(
  physeq = Gwensmith_scale, 
  method = "PCoA", 
  distance = "unifrac"
)
##I get a warning but likely because the tree contained bacteria from all of the samples (bobbys, Rheas etc.). But I am just looking at mine. So the function is alerting me that I am only getting a subset. 
# Plot 
windows(width=10, height=6)
plot_ordination(
  physeq = Gwensmith_scale,
  ordination = Gwensmith_pcoa,
  axes = 1:2,
  color = "Smithora_status",
  shape = "before_after_transplant",
  title = "PCoA of Smithora"
) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                "#4daf4a", "#1919ff", "darkorchid3", "magenta")
  ) +
  geom_point(aes(color = Smithora_status), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) +theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3))

##ambients cluster with samples taken after the transplant.Time has a strong effect. ##Before then transplant, shoots that had or didn't have smithora had more distinct communities. Then after the transplant everything is less distinct. Ambients that werent transplanted still maintain soe distinction which is cool. I think harming the shooot compromises its ability to maintain its own bacterial community. 

##The guy in the blog did this. I think its to get the same 2D image each time. 
set.seed(1)

# Ordinate
Gwensmith_nmds <- ordinate(
  physeq = Gwensmith_scale, 
  method = "NMDS", 
  distance = "unifrac"
)

windows(width=10, height=6)
plot_ordination(
  physeq = Gwensmith_scale,
  ordination = Gwensmith_nmds,
  color = "Smithora_status",
  shape = "before_after_transplant",
  title = "NMDS of Smithora"
) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                "#4daf4a", "#1919ff", "darkorchid3", "magenta")
  ) +
  geom_point(aes(color = Smithora_status), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) + 
  theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3))
###Similar patterns to PCA
##does smithora lower your diversity???? It seems like things without it are way more spread out

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



##This is whta phyloseq produces
boxplots <- plot_richness(Gwensmith, x = "before_after_transplant", color = "Smithora_status") + geom_boxplot()
boxplots <- boxplots + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Diversity Value")) + xlab(expression("Type of Sample")) + ggtitle("Diversity Values from Various Measures") + theme(plot.title = element_text(vjust= 2, size = 14)) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3)) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
windows(width=10, height=6) 
print(boxplots)


diversityboxplot <- ggplot(alphadiv, aes(x = before_after_transplant, y = mean, color = Smithora_status)) + geom_boxplot(size = 1) + facet_wrap(~measure, ncol = 1, scales = "free") + scale_color_manual(values = c("#E96446", "#302F3D", "#87CEFA")) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Diversity Value")) + xlab(expression("Type of Sample")) + ggtitle("Diversity Values from Random Sampling Method") + theme(plot.title = element_text(vjust= 2, size = 14)) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3)) 
windows(width=10, height=6)
print(diversityboxplot)


##Make a constrained Ordination
# Remove data points with missing metadata
Gwensmith_not_na <- Gwensmith_scale %>%
  subset_samples(
    !is.na(Smithora_status) & 
      !is.na(before_after_transplant)
  )

unifrac_gwen <- phyloseq::distance(physeq = Gwensmith_not_na, method = "unifrac")


# CAP ordinate
cap_ord <- ordinate(
  physeq = Gwensmith_not_na,
  method = "CAP",
  distance = unifrac_gwen,
  formula = ~ Smithora_status + before_after_transplant
)

##I guess anova tells you if the axes you chose are significant. I cannot get my caps to map though. 
anova(cap_ord)

# CAP plot
cap_plot <- plot_ordination(
  physeq = Gwensmith_not_na, 
  ordination = cap_ord, 
  color = "Smithora_status", 
  axes = c(1,2)
) + 
  aes(shape = before_after_transplant) + 
  geom_point(aes(colour = Smithora_status), alpha = 0.4, size = 4) + 
  geom_point(colour = "grey90", size = 1.5) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", 
                                "#1919ff", "darkorchid3", "magenta"))

print(cap_plot)
# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display= "cn")
?vegan::scores()
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
arrowdf
rownames(arrowmat)

# Add labels, make a data.frame
arrowdf <- data.frame(arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)
print(arrow_map)
plot(arrow_map)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)
print(label_map)
plot(label_map)

arrowhead = arrow(length = unit(0.02, "npc"))
dev.off()
# Make a new graphic
newcapplot <- cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )

newcapplot <- newcapplot + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ggtitle("Constrained PCA") + theme(plot.title = element_text(vjust= 2, size = 14)) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3)) 

windows(width=10, height=6)
print(newcapplot)


dev.off()

##Try Permanova again

