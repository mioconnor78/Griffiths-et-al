
###Code for reading in OTU table, tree, map 

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
##make phylum graph. Grouped by phylum first get rid of taxa at less than 0.02 abundance
Gwensmith_longPhylum <- Gwensmith %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance 
  arrange(Phylum) 
##Lets try getting the abundances the same
Gwensmithrar <- rarefy_even_depth(Gwensmith, sample.size = 7244, verbose = FALSE, replace = TRUE)
Gwensmithrar_longPhylum <- Gwensmithrar %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance 
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
Gwensmith_longPhylum_justexp <- filter(Gwensmith_longPhylum, before_after_transplant != "ambient")
Gwensmith_longPhylum_justexp <- filter(Gwensmith_longPhylum_justexp, blade_ID != "edge_d")
Gwensmith_longPhylum_controls <- filter(Gwensmith_longPhylum_justexp, type_swab == "control")
Gwensmith_longPhylum_experiments <- filter(Gwensmith_longPhylum_justexp, type_swab == "experiment")
Gwensmith_longPhylum_ambient <- filter(Gwensmith_longPhylum, before_after_transplant== "ambient" & host_species == "Zostera")
Gwensmith_longPhylum_smith <- filter(Gwensmith_longPhylum, host_species == "Smithora")




tryingbarexp <- ggplot(Gwensmith_longPhylum_experiments)

tryingbarexp <- ggplot(Gwensmith_longPhylum_experiments, aes(x = blade_ID, y = Abundance, fill = Phylum)) + geom_bar(stat = "identity") + scale_fill_manual(values = phylum_colors) + facet_wrap(~before_after_transplant) + theme(axis.text.x = element_text(angle = 60, size = 12, colour = "black",hjust = 1))

tryingbarexp <- tryingbar + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Relative Abundance")) + xlab(expression("Blade ID")) + ggtitle("Relative abundance of bacterial species by Phylum ") + theme(plot.title = element_text(vjust= 2, size = 14)) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3))

##open up new window
windows(width=10, height=6)
print(tryingbarexp)

###try with just controls
tryingbarcont <- ggplot(Gwensmith_longPhylum_controls)

tryingbarcont <- ggplot(Gwensmith_longPhylum_controls, aes(x = blade_ID, y = Abundance, fill = Phylum)) + geom_bar(stat = "identity") + scale_fill_manual(values = phylum_colors) + facet_wrap(~before_after_transplant) + theme(axis.text.x = element_text(angle = 60, size = 12, colour = "black",hjust = 1))

tryingbarcont <- tryingbar + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Relative Abundance")) + xlab(expression("Blade ID")) + ggtitle("Relative abundance of bacterial species by Phylum ") + theme(plot.title = element_text(vjust= 2, size = 14)) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3))

##open up new window
windows(width=10, height=6)
print(tryingbarcont)

tryingbaramb <- ggplot(Gwensmith_longPhylum_ambient)

tryingbaramb <- ggplot(Gwensmith_longPhylum_ambient, aes(x = blade_ID, y = Abundance, fill = Phylum)) + geom_bar(stat = "identity") + scale_fill_manual(values = phylum_colors) + facet_wrap(~host_species) + theme(axis.text.x = element_text(angle = 60, size = 12, colour = "black",hjust = 1))

tryingbaramb <- tryingbaramb + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + ylab(expression("Relative Abundance")) + xlab(expression("Blade ID")) + ggtitle("Relative abundance of bacterial species by Phylum ") + theme(plot.title = element_text(vjust= 2, size = 14)) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3))

##open up new window
windows(width=10, height=6)
print(tryingbaramb)










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
Gwensmith_longfamily <- Gwensmithrar %>%
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


###See bacterial diversity plots R. file that uses this family data set to display bacterial community composition 

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
###Gwensmith2 <- Gwensmith %>%
  ###subset_samples(Smithora_status != "smithora")
Gwensmith_scale <- scale_reads(Gwensmith,100)
##Gwensmith_scale2 <- Gwensmith_scale %>%
  ##subset_samples(Smithora_status != "smithora")


# Fix smithora levels in sample_data
sample_data(Gwensmith_scale)$Smithora_status <- factor(
  sample_data(Gwensmith_scale)$Smithora_status, 
  levels = c("absent","present","smithora")
)
sample_data(Gwensmith_scale)$before_after_transplant <- factor(
  sample_data(Gwensmith_scale)$before_after_transplant, 
  levels = c("before","after","ambient")
)

View(sample_data(Gwensmith_scale))

##Try and Add a new row in the data set that will better explain where the bacterial samples are coming from

sample_data(Gwensmith_scale)$Community <- c("Ambient Zostera Edge","Ambient Zostera Edge","Ambient Zostera Interior","Ambient Zostera Interior","After Detachment Edge","After Detachment and Transplant to Interior","After Detachment Edge","After Detachment Edge","After Detachment and Transplant to Interior","After Detachment and Transplant to Interior","After Detachment and Transplant to Edge","After Detachment and Transplant to Edge","After Detachment Interior","After Detachment Interior","After Detachment and Transplant to Edge","After Detachment Interior","Before Detachment Edge","Before Detachment and Transplant to Interior","Before Detachment Edge","Before Detachment and Transplant to Interior","Before Detachment and Transplant to Interior","Before Detachment and Transplant to Edge","Before Detachment and Transplant to Edge","Before Detachment Interior","Before Detachment Interior","Before Detachment and Transplant to Interior","Before Detachment Interior","Ambient Smithora Blade","Ambient Smithora Blade","Ambient Smithora Blade","Ambient Smithora Blade")

sample_data(Gwensmith_scale)$Time <- c("After","After","After","After","After","After","After","After","After","After","After","After","After","After","After","After","Before","Before","Before","Before","Before","Before","Before","Before","Before","Before","Before","Before","Before","Before","Before")
View(sample_data(Gwensmith_scale))

##check that my physeq file is still in order
Gwensmith_scale
Gwensmith_scale2

# Ordinate look at an uncontrained PCA with unifrac
Gwensmith_pcoa <- ordinate(
  physeq = Gwensmith_scale, 
  method = "PCoA", 
  distance = "unifrac"
)
##Do the same with smithora removed to see if it affects the dispersion.

Gwensmith2_pcoa <- ordinate(
  physeq = Gwensmith_scale2, 
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
  title = "PCoA of Smithora"
) + 
  scale_color_manual(values = c("#c7a2ab","#2bf04a","#403fe8","#76e70f","#7f44f2","#00f483","#b738eb","#82ff8c","#0017a6","#fff64b","#0165f5","#4e3500","#66002e","#290a00")) +
  geom_point(aes(color = Community), alpha = 0.7, size = 4) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3))

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

NMDSplot <- plot_ordination(
  physeq = Gwensmith_scale,
  ordination = Gwensmith_nmds,
  title = "Surface bacterial communities",
  color = "Community",
  shape = "Time")
NMDSplot <- NMDSplot + geom_point(size = 5) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill="white")) + theme(axis.title.x = element_text(vjust = -1, size =14)) + theme(axis.title.y = element_text(vjust = -1, size = 14)) + theme(strip.text.y = element_text(face = 3)) 

windows(width=10, height=6)
print(NMDSplot)

windows(width=10, height=6)
p2 = plot_ordination(Gwensmith_scale, Gwensmith_nmds, type="Samples", color="Community_Source", shape = "before_after_transplant") 
p2 + geom_polygon(aes(fill=Community_Source)) + geom_point(size=5) + ggtitle("Surface Bacterial Communities")

windows(width=10, height=6)
print(p2)
##Use NMDS because it can use any similarity metric like unifrac

###Similar patterns to PCA
##does smithora lower your diversity???? It seems like things without it are way more spread out

##Try Permanova again

set.seed(1)

# Calculate bray curtis distance matrix
bacteria_unifrac <- phyloseq::distance(Gwensmith_scale, method = "unifrac")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(Gwensmith_scale))

# Adonis test
adonis(bacteria_unifrac ~ before_after_transplant*Smithora_status, data = sampledf)

##Beta dispersion test for transplant
transplantcomm <- betadisper(bacteria_unifrac, sampledf$Community)
permutest(transplantcomm)
##our betadisper results are not significant, meaning we cannot reject the null hypothesis that our groups have the same dispersions. This means we can be more confident that our adonis result is a real result, and not due to differences in group dispersions
##Beta dispersion test for smithora status
smithorabeta <- betadisper(bacteria_unifrac, sampledf$Smithora_status)
permutest(smithorabeta)


