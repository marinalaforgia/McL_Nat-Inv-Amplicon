### 16S and ITS Analysis April 2020 ###
rm(list=ls())
#### Setup ####

#setwd("/Users/Marina/Documents/UC-Davis/Projects/McL_Nat-Inv-Amplicon") # Marina's Working Directory
#setwd("/Users/Cassie/Documents/Dropbox/Greenhouse_ITS/McL_Nat-Inv-Amplicon") # Cassie's Working Directory

### Load Libraries ###

library(EcolUtils)
library(ggplot2)
#library(plyr)
library(multcomp)
library(FSA)
#library(Rmisc)
library(emmeans)
library(lme4)
library(lmerTest)
library(knitr)
#library(qiimer) not on cran anymore?
library(vegan)
library(ape)
library(Hmisc)
library(phyloseq)
library(RColorBrewer)
library(coin)
library(rmarkdown)
library(reshape)
library(betapart)
library(dada2)
library(extrafont)
library(tidyverse) # for tidy data packages, automatically loads dplyr
library(magrittr) # for piping
library(DESeq2)
library(biobroom)
library(patchwork)
library(VennDiagram)

### Read in Data ###

## Remove background soil for analysis ##

### Not rarefied data ###

## 16S ##
ps.16s.nocontrols <- readRDS("Data/16S/Intermediates/ps-nocontrols-unrare.RDS") # 16s data, Marina's path
ps.16s.nocontrols <- readRDS("Data/ps-nocontrols-unrare.RDS") #16s data, Cassie's path
#132 samples

## ITS ##
ps.ITS.nocontrols <- readRDS("Data/ITS/greenhouse_its_itsx_decontam_controlsremoved.rds") # ITS data
ps.ITS.nocontrols <- readRDS("Data/greenhouse_its_itsx_decontam_controlsremoved.rds") # ITS data, Cassie's path
#132 samples

# Metadata mapping file including biomass and traits (version4)
mapping <- read.csv("/Users/Marina/Documents/UC-Davis/Projects/McL_Nat-Inv-Amplicon/Data/Setup/Grassland-Amplicon-Mapping-File4.csv") 
#mapping <- read.csv("Data/Grassland-Amplicon-Mapping-File4.csv") #Cassie's path

#Adding columns to the mapping file
mapping$FunGroup <- ifelse(mapping$TreatmentName == "Invasive_grass", "Grass", NA)
mapping$FunGroup <- ifelse(mapping$Competion == "SingleSpecies" & mapping$TreatmentName != "Invasive_grass", "Forb", mapping$FunGroup)
mapping[is.na(mapping$FunGroup),]$FunGroup <- "grass_x_forb"
mapping$FunGroup <- as.factor(mapping$FunGroup)

# 16S mapping file does not match up with the mapping file for the ITS data, adjust here 
for(i in sample_names(ps.16s.nocontrols)) {
  sample_names(ps.16s.nocontrols)[which(sample_names(ps.16s.nocontrols) == i)] <- as.character(mapping[mapping$SampleID.16S == i,]$SampleID_Fix)
}

sample_names(ps.16s.nocontrols)
sample_names(ps.ITS.nocontrols) # they match up!

df.16s <- data.frame(sample_data(ps.16s.nocontrols))
df.ITS <- data.frame(sample_data(ps.ITS.nocontrols))

missing <- df.ITS[!(df.ITS$SampleID_Fix %in% df.16s$SampleID_Fix),]
missing2 <- df.16s[!(df.16s$SampleID_Fix %in% df.ITS$SampleID_Fix),]
row.names(mapping) <- mapping$SampleID_Fix
sample_data(ps.ITS.nocontrols) <- mapping[which(mapping$SampleID_Fix %in% ps.ITS.nocontrols@sam_data$SampleID_Fix),]
sample_data(ps.ITS.nocontrols) # fixed!

sample_data(ps.16s.nocontrols) <- mapping

rm(missing, missing2, df.16s, df.ITS, i)


### Fix Taxonomy for DESeq analysis later ###

# Rename NA's & remove prefixes if they exist (e.g. p__)
# This is important if we want to bind taxonomy to deseq results later
# 16S dataset has NA's but no prefixes, ITS has both

## 16S ##

df.16s.tax <- data.frame(tax_table(ps.16s.nocontrols))

df.16s.tax %<>% 
  mutate(Phylum = fct_explicit_na(Phylum, na_level = "Unclassified"), 
         Class = fct_explicit_na(Class, na_level = "Unclassified"), 
         Order = fct_explicit_na(Order, na_level = "Unclassified"), 
         Family = fct_explicit_na(Family, na_level = "Unclassified"), 
         Genus = fct_explicit_na(Genus, na_level = "Unclassified"), 
         Species = fct_explicit_na(Species, na_level = "Unclassified"))

row.names(df.16s.tax) <- row.names(tax_table(ps.16s.nocontrols))
tax.16s <- as.matrix(df.16s.tax)

tax_table(ps.16s.nocontrols) <- tax.16s

## ITS ##

df.ITS.tax <- data.frame(tax_table(ps.ITS.nocontrols))

df.ITS.tax %<>% 
  mutate(Phylum = fct_explicit_na(Phylum, na_level = "p__Unclassified"), 
         Class = fct_explicit_na(Class, na_level = "c__Unclassified"), 
         Order = fct_explicit_na(Order, na_level = "o__Unclassified"), 
         Family = fct_explicit_na(Family, na_level = "f__Unclassified"), 
         Genus = fct_explicit_na(Genus, na_level = "g__Unclassified"), 
         Species = fct_explicit_na(Species, na_level = "s__Unclassified"))

tax.list <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
tax.header <- c(Phylum = "p__", Class = "c__", Order = "o__", 
                Family = "f__", Genus = "g__", Species = "s__")

for (i in tax.list) {
  names <- sapply(strsplit(as.character(df.ITS.tax[[i]]), as.character(tax.header[[i]])), `[`, 2)
  df.ITS.tax[[i]] <- names 
}

row.names(df.ITS.tax) <- row.names(tax_table(ps.ITS.nocontrols))
ITS.tax <- as.matrix(df.ITS.tax)

tax_table(ps.ITS.nocontrols) <- ITS.tax


### Rarefied Data ###

## 16S ##
ps.nocontrols.16s.rare9434 <- rarefy_even_depth(ps.16s.nocontrols, 9434, replace = FALSE, rngseed = 5311)
#saveRDS(ps.nocontrols.16s.rare9434, "Data/ps-nocontrols-rare9434.RDS")
# 131 samples, lose one sample, GLM-0116

# Save for first pass at Halla (https://huttenhower.sph.harvard.edu/halla)
#write.table(t(otu_table(ps.nocontrols.16s.rare9434)), "Data/ps.nocontrols.16s.rare9434.transposed.txt", sep = "\t")

#Save for first pass sourcetracker (https://github.com/danknights/sourcetracker/)
#write.csv(otu_table(ps.nocontrols.16s.rare9434), "Data/ps.nocontrols.16s.rare9434.csv")
#write.csv(sample_data(ps.nocontrols.16s.rare9434), "Data/ps.nocontrols.16s.rare9434.mappingdata.csv")

# Load 
ps.16s.nocontrols.rare <- readRDS("Data/16S/Intermediates/ps-nocontrols-rare9434.RDS") # rarefied 16s data
ps.16s.nocontrols.rare <- readRDS("Data/ps-nocontrols-rare9434.RDS") # rarefied 16s data, Cassie's path


## ITS ##
ps.nocontrols.its.rare7557 <- rarefy_even_depth(ps.ITS.nocontrols, 7557, replace = FALSE, rngseed = 5311)
#saveRDS(ps.nocontrols.its.rare7557, "Data/greenhouse_its_itsx_decontam_controlsremoved_rare7557.rds")
# 132 samples

# Save for first pass at Halla (https://huttenhower.sph.harvard.edu/halla)
#\write.table(t(otu_table(ps.nocontrols.its.rare7557)), "Data/ps.nocontrols.its.rare7557.transposed.txt", sep = "\t")



#Save for first pass sourcetracker (https://github.com/danknights/sourcetracker/)
#write.csv(otu_table(ps.nocontrols.its.rare7557), "Data/ps.nocontrols.its.rare7557.csv")
#write.csv(sample_data(ps.nocontrols.its.rare7557), "Data/ps.nocontrols.its.rare7557.mappingdata.csv")


# Load 
ps.ITS.nocontrols.rare <- readRDS("Data/ITS/greenhouse_its_itsx_decontam_controlsremoved_rare7557.rds") # rarefied ITS data
ps.ITS.nocontrols.rare <- readRDS("Data/greenhouse_its_itsx_decontam_controlsremoved_rare7557.rds") # rarefied ITS data, Cassie's path



#### Prep Datasets ####

## 16s ##

ps.16s.nocontrols.rare.spp <- subset_samples(ps.16s.nocontrols.rare, Competion != "TwoSpecies")

ps.16s.nocontrols.rare.spp.RA <- transform_sample_counts(ps.16s.nocontrols.rare.spp, function(x) x / sum(x)) # relative abundance (used later)

## ITS ##

ps.ITS.nocontrols.rare.spp <- subset_samples(ps.ITS.nocontrols.rare, Competion != "TwoSpecies") 

ps.ITS.nocontrols.rare.spp.RA <- transform_sample_counts(ps.ITS.nocontrols.rare.spp, function(x) x / sum(x)) 

#### Ordination (16S): Soil v Rhiz ####
# Confirm whether rhizosphere samples are different from soil

# Weighted Unifrac (RA, tree)
ps.16s.nc.rare.soil.ord <- ordinate(ps.16s.nocontrols.rare, "PCoA", "wunifrac")

plot_ordination(ps.16s.nocontrols.rare, ps.16s.nc.rare.soil.ord, color = "SampleSubType")  +
  theme_bw(base_size = 20) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = SampleSubType))

DistWU.Soil <- phyloseq::distance(ps.16s.nocontrols.rare, method = "wunifrac", type = "samples")

adonis(DistWU.Soil ~ SampleSubType, as(sample_data(ps.16s.nocontrols.rare), "data.frame"), permutations = 9999) # rhizosphere is sig different from soil

C <- betadisper(DistWU.Soil, as(sample_data(ps.16s.nocontrols.rare), "data.frame")$SampleSubType)
permutest(C, permutations = 9999) # but their spread varies.. so... who knows

rm(C, DistWU.Soil, ps.16s.nc.rare.soil.ord)

#### Ordination (16S): G v SA v ST ####

# Weighted Unifrac (RA, tree)
ps.16s.nc.rare.spp.ord <- ordinate(ps.16s.nocontrols.rare.spp, "PCoA", "wunifrac")

# PCoA Plot
plot_ordination(ps.16s.nocontrols.rare.spp, ps.16s.nc.rare.spp.ord, color = "TreatmentName", shape = "TreatmentName") +
  theme_bw(base_size = 20) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = TreatmentName)) 

DistWU.Spp <- phyloseq::distance(ps.16s.nocontrols.rare.spp, method = "wunifrac", type = "samples")

adonis(DistWU.Spp ~ TreatmentName, as(sample_data(ps.16s.nocontrols.rare.spp), "data.frame"), permutations = 9999) 

C <- betadisper(DistWU.Spp, as(sample_data(ps.16s.nocontrols.rare.spp), "data.frame")$TreatmentName)
permutest(C, permutations = 9999) # spread doesnt vary between competition treatments

pair.am.spp <- adonis.pair(vegdist(DistWU.Spp), as(sample_data(ps.16s.nocontrols.rare.spp), "data.frame")$TreatmentName, corr.method = "BH")
kable(pair.am.spp, caption = "Differences between Functional Groups")

# not sig different...

#### Ordination (16S): G v F ####
ps.16s.nc.rare.spp.ord <- ordinate(ps.16s.nocontrols.rare.spp, "PCoA", "wunifrac")

plot_ordination(ps.16s.nocontrols.rare.spp, ps.16s.nc.rare.spp.ord, shape = "FunGroup", color = "FunGroup") +
  theme_bw(base_size = 20) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = FunGroup)) 

adonis(DistWU.Spp ~ FunGroup, as(sample_data(ps.16s.nocontrols.rare.spp), "data.frame"), permutations = 9999)

C <- betadisper(DistWU.Spp, as(sample_data(ps.16s.nocontrols.rare.spp), "data.frame")$FunGroup)
permutest(C, permutations = 9999) # spread doesnt vary between functional group. grasses and forbs have significantly different communities

#### Ordination (16S): Competition ####
ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion)))

ps.rare.ord.tr <- ordinate(ps.16s.nocontrols.rare, "PCoA", "wunifrac")

plot_ordination(ps.16s.nocontrols.rare, ps.rare.ord.tr, color = "FunGroup", shape = "Competion") +
  theme_bw(base_size = 20) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = FunGroup)) 

DistWU <- phyloseq::distance(ps.16s.nocontrols.rare, method = "wunifrac", type = "samples")

# adonis compares centroids and is sensitive to differences in spread

adonis(DistWU ~ Competion, as(sample_data(ps.16s.nocontrols.rare), "data.frame"), permutations = 9999) # communities shift when going from single species to two species

adonis(DistWU ~ TreatmentName, as(sample_data(ps.16s.nocontrols.rare), "data.frame"), permutations = 9999)

adonis(DistWU ~ FunGroup, as(sample_data(ps.16s.nocontrols.rare), "data.frame"), permutations = 9999) # I think this should be our focus, no on the SA/ST

# Looks like spread doesnt vary significantly so we can be assured that our centroids are telling us something about these are different
C <- betadisper(DistWU, as(sample_data(ps.16s.nocontrols.rare), "data.frame")$Competion)
permutest(C, permutations = 9999) # spread doesnt vary between competition treatments

C <- betadisper(DistWU, as(sample_data(ps.16s.nocontrols.rare), "data.frame")$TreatmentName)
permutest(C, permutations = 9999) # spread doesnt vary between treatments

C <- betadisper(DistWU, as(sample_data(ps.16s.nocontrols.rare), "data.frame")$FunGroup)
permutest(C, permutations = 9999) # spread doesnt vary between treatments

pair.am2 <- adonis.pair(vegdist(DistWU), as(sample_data(ps.16s.nocontrols.rare), "data.frame")$TreatmentName, corr.method = "none")
pair.am2 <- pair.am2[-c(8,7),1:6]
pair.am2 <- cbind(pair.am2, p.val.corr = p.adjust(pair.am2$P.value, method = "BH", n = length(pair.am2$P.value)))
kable(pair.am2)

pair.am3 <- adonis.pair(vegdist(DistWU), as(sample_data(ps.16s.nocontrols.rare), "data.frame")$FunGroup, corr.method = "BH")
kable(pair.am3)

#### Ordination (ITS): Soil v Rhiz ####
GL_pcoa <- ordinate(
  physeq = ps.ITS.nocontrols.rare, 
  method = "PCoA", 
  distance = "bray")

plot_ordination(
  physeq = ps.ITS.nocontrols.rare,
  ordination = GL_pcoa,
  shape = "SampleSubType",
  color = "SampleSubType") + 
  geom_point(size = 6) + 
  theme(text = element_text(size=28)) 

DistBC = phyloseq::distance(ps.ITS.nocontrols.rare, method = "bray", type="samples")

adonis(DistBC ~ SampleSubType, as(sample_data(ps.ITS.nocontrols.rare), "data.frame"), permutations = 9999) # rhizosphere soil is diff from background soil

A <- betadisper(DistBC, as(sample_data(ps.ITS.nocontrols.rare), "data.frame")$SampleSubType)

permutest(A, permutations = 9999) 

#### Ordination (ITS): G v SA v ST ####

GL_pcoa <- ordinate(
  physeq = ps.ITS.nocontrols.rare.spp, 
  method = "PCoA", 
  distance = "bray")

plot_ordination(
  physeq = ps.ITS.nocontrols.rare.spp,
  ordination = GL_pcoa,
  shape = "TreatmentName",
  color = "TreatmentName") + 
  geom_point(size = 6) + 
  theme(text = element_text(size=28)) 

DistBC = phyloseq::distance(ps.ITS.nocontrols.rare.spp, method = "bray", type="samples")

adonis(DistBC ~ TreatmentName, as(sample_data(ps.ITS.nocontrols.rare.spp), "data.frame"), permutations = 9999) # nothing

A <- betadisper(DistBC, as(sample_data(ps.ITS.nocontrols.rare.spp), "data.frame")$TreatmentName)
permutest(A, permutations = 9999)

#### Ordination (ITS): G v F ####
GL_pcoa <- ordinate(
  physeq = ps.ITS.nocontrols.rare.spp, 
  method = "PCoA", 
  distance = "bray")

plot_ordination(
  physeq = ps.ITS.nocontrols.rare.spp,
  ordination = GL_pcoa,
  shape = "FunGroup",
  color = "FunGroup") + 
  geom_point(size = 6) + 
  theme(text = element_text(size=28)) 

DistBC = phyloseq::distance(ps.ITS.nocontrols.rare.spp, method = "bray", type = "samples")

adonis(DistBC ~ FunGroup, as(sample_data(ps.ITS.nocontrols.rare.spp), "data.frame"), permutations = 9999) # no differences


A <- betadisper(DistBC, as(sample_data(ps.ITS.nocontrols.rare.spp), "data.frame")$FunGroup)
permutest(A, permutations = 9999)

#### Ordination (ITS): Competition ####
ps.ITS.nocontrols.rare <- subset_samples(ps.ITS.nocontrols.rare, !(is.na(Competion)))

GL_pcoa <- ordinate(
  physeq = ps.ITS.nocontrols.rare, 
  method = "PCoA", 
  distance = "bray")

plot_ordination(physeq = ps.ITS.nocontrols.rare, ordination = GL_pcoa, color = "FunGroup", shape = "Competion") +
  geom_point(size = 6) + 
  theme(text = element_text(size=28)) +
  stat_ellipse(aes(group = FunGroup)) 

DistBC = phyloseq::distance(ps.ITS.nocontrols.rare, method = "bray", type="samples")

adonis(DistBC ~ TreatmentName, as(sample_data(ps.ITS.nocontrols.rare), "data.frame"), permutations = 9999) # marginal treatment diffs

adonis(DistBC ~ Competion, as(sample_data(ps.ITS.nocontrols.rare), "data.frame"), permutations = 9999) # marginal competition shifts

adonis(DistBC ~ FunGroup, as(sample_data(ps.ITS.nocontrols.rare), "data.frame"), permutations = 9999) # marginal fungroup diffs

# because all these are marginal, we should not be looking at pairwise comparisons


#### ASV (16S): G v F ####
ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion)))
ps.16s.nocontrols.rare.RA <- transform_sample_counts(ps.16s.nocontrols.rare, function(x) x / sum(x))

AvgRA99_ASV <- filter_taxa(ps.16s.nocontrols.rare.RA, function(x) mean(x) > .01, TRUE)
df_ASV <- psmelt(AvgRA99_ASV)
grouped_ASV <- group_by(df_ASV, FunGroup, OTU, Genus, Family, Order, Class, Phylum, Kingdom)
avgs_ASV <- summarise(grouped_ASV, mean=100*mean(Abundance), sd=100*sd(Abundance), se=100*se(Abundance))

ggplot(avgs_ASV, aes(x = OTU, y = mean, fill = Genus)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se)), 
                width = .4, position = position_dodge(.9)) +
  facet_wrap(~FunGroup, nrow = 3) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 20)) + 
  ylab("Mean Relative Abundance")


listofcats_sig = NULL
DataSet = df_ASV
taxa = unique(DataSet$OTU)
df_taxa <- data.frame(taxa = taxa, chisq = rep(NA, length(taxa)), pvals = rep(NA, length(taxa)))

for (cat in df_taxa$taxa) { # for each ASV that has a mean abundance above 1%
  new_df <- subset(DataSet, DataSet$OTU == cat)
  new_df$ST <- as.factor(new_df$FunGroup)
  kw = kruskal.test(Abundance ~ ST, data=new_df) # test whether the abundance of that ASV varies by ST, defined above as functional group
  df_taxa[df_taxa$taxa == cat,]$chisq = kw$statistic # return the chi square stat
  df_taxa[df_taxa$taxa == cat,]$pvals = kw$p.value # and the p-value
  
  if (kw$p.value <= 0.05) { # filter by significant p values
    listofcats_sig = c(cat, listofcats_sig)

  }
}

df_taxa$pvals.bh <- p.adjust(df_taxa$pvals, method="BH")
kable(df_taxa) # this tells us which ASVs vary by functional group

# now let's figure out which functional groups
pvals.dunn = NULL
pvals.dunn.bh = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$OTU == cat)
  new_df$ST <- as.factor(new_df$FunGroup)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bh") # bh adjustment done here, no need to do it below
  
  for (i in 1:length(dT$res$Comparison)) {

    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bh = c(dT$res$P.adj[i], pvals.dunn.bh)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(as.character(dT$res$Comparison)[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)

  }

}

df_taxa <- data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bh)
kable(df_taxa, caption = "Significant Differences between Functional Groups")


#### ASV (16S): G v SA v ST ####
AvgRA99_ASV <- filter_taxa(ps.16s.nocontrols.rare.RA, function(x) mean(x) > .01, TRUE)
df_ASV <- psmelt(AvgRA99_ASV)
grouped_ASV <- group_by(df_ASV, TreatmentName, OTU, Genus, Family, Order, Class, Phylum, Kingdom)
avgs_ASV <- summarise(grouped_ASV, mean=100*mean(Abundance), sd=100*sd(Abundance), se=100*se(Abundance))

ggplot(avgs_ASV, aes(x = OTU, y = mean, fill = Order)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se)), 
                width = .4, position = position_dodge(.9)) +
  facet_wrap(~TreatmentName, nrow = 3) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 20)) + 
  ylab("Mean Relative Abundance")

listofcats_sig = NULL
DataSet = df_ASV
taxa = unique(DataSet$OTU)
df_taxa <- data.frame(taxa = taxa, chisq = rep(NA, length(taxa)), pvals = rep(NA, length(taxa)))

for (cat in df_taxa$taxa) { # for each ASV that has a mean abundance above 1%
  new_df <- subset(DataSet, DataSet$OTU == cat)
  new_df$ST <- as.factor(new_df$TreatmentName)
  kw = kruskal.test(Abundance ~ ST, data=new_df) # test whether the abundance of that ASV varies by ST, defined above as functional group
  df_taxa[df_taxa$taxa == cat,]$chisq = kw$statistic # return the chi square stat
  df_taxa[df_taxa$taxa == cat,]$pvals = kw$p.value # and the p-value
  
  if (kw$p.value <= 0.05) { # filter by significant p values
    listofcats_sig = c(cat, listofcats_sig)

  }
}

df_taxa$pvals.bh <- p.adjust(df_taxa$pvals, method="BH")
kable(df_taxa) # this tells us which ASVs vary by functional group

# now let's figure out which functional groups
pvals.dunn = NULL
pvals.dunn.bh = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$OTU == cat)
  new_df$ST <- as.factor(new_df$TreatmentName)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bh") # bh adjustment done here, no need to do it below
  
  for (i in 1:length(dT$res$Comparison)) {

    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bh = c(dT$res$P.adj[i], pvals.dunn.bh)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(as.character(dT$res$Comparison)[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)

  }

}

df_taxa <- data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bh)
kable(df_taxa, caption = "Significant Differences between Functional Groups")


#### ASV (ITS): G v F ####
ps.ITS.nocontrols.rare <- subset_samples(ps.ITS.nocontrols.rare, !(is.na(Competion)))
ps.ITS.nocontrols.rare.RA <- transform_sample_counts(ps.ITS.nocontrols.rare, function(x) x / sum(x)) 

AvgRA99_ASV <- filter_taxa(ps.ITS.nocontrols.rare.RA, function(x) mean(x) > .015, TRUE)
df_ASV <- psmelt(AvgRA99_ASV)
grouped_ASV <- group_by(df_ASV, FunGroup, OTU, Order, Class, Phylum, Kingdom)
avgs_ASV <- summarise(grouped_ASV, mean=100*mean(Abundance), sd=100*sd(Abundance), se=100*se(Abundance))

ggplot(avgs_ASV, aes(x = OTU, y = mean, fill = Order)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se)), 
                width = .4, position = position_dodge(.9)) +
  facet_wrap(~FunGroup, nrow = 3) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 20)) + 
  ylab("Mean Relative Abundance")


listofcats_sig = NULL
DataSet = df_ASV
taxa = unique(DataSet$OTU)
df_taxa <- data.frame(taxa = taxa, chisq = rep(NA, length(taxa)), pvals = rep(NA, length(taxa)))

for (cat in df_taxa$taxa) { # for each ASV that has a mean abundance above 1%
  new_df <- subset(DataSet, DataSet$OTU == cat)
  new_df$ST <- as.factor(new_df$FunGroup)
  kw = kruskal.test(Abundance ~ ST, data=new_df) # test whether the abundance of that ASV varies by ST, defined above as functional group
  df_taxa[df_taxa$taxa == cat,]$chisq = kw$statistic # return the chi square stat
  df_taxa[df_taxa$taxa == cat,]$pvals = kw$p.value # and the p-value
  
  if (kw$p.value <= 0.05) { # filter by significant p values
    listofcats_sig = c(cat, listofcats_sig)

  }
}

df_taxa$pvals.bh <- p.adjust(df_taxa$pvals, method="BH")
kable(df_taxa) # this tells us which ASVs vary by functional group

# now let's figure out which functional groups
pvals.dunn = NULL
pvals.dunn.bh = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$OTU == cat)
  new_df$ST <- as.factor(new_df$FunGroup)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bh") # bh adjustment done here, no need to do it below
  
  for (i in 1:length(dT$res$Comparison)) {

    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bh = c(dT$res$P.adj[i], pvals.dunn.bh)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(as.character(dT$res$Comparison)[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)

  }

}

df_taxa <- data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bh)
kable(df_taxa, caption = "Significant Differences between Functional Groups")

#### ASV (ITS): G v SA v ST ####
AvgRA99_ASV <- filter_taxa(ps.ITS.nocontrols.rare.RA, function(x) mean(x) > .015, TRUE)
df_ASV <- psmelt(AvgRA99_ASV)
grouped_ASV <- group_by(df_ASV, TreatmentName, OTU, Order, Class, Phylum, Kingdom)
avgs_ASV <- summarise(grouped_ASV, mean=100*mean(Abundance), sd=100*sd(Abundance), se=100*se(Abundance))

ggplot(avgs_ASV, aes(x = OTU, y = mean, fill = Order)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se)), 
                width = .4, position = position_dodge(.9)) +
  facet_wrap(~TreatmentName, nrow = 3) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 20)) + 
  ylab("Mean Relative Abundance")


listofcats_sig = NULL
DataSet = df_ASV
taxa = unique(DataSet$OTU)
df_taxa <- data.frame(taxa = taxa, chisq = rep(NA, length(taxa)), pvals = rep(NA, length(taxa)))

for (cat in df_taxa$taxa) { # for each ASV that has a mean abundance above 1%
  new_df <- subset(DataSet, DataSet$OTU == cat)
  new_df$ST <- as.factor(new_df$TreatmentName)
  kw = kruskal.test(Abundance ~ ST, data=new_df) # test whether the abundance of that ASV varies by ST, defined above as functional group
  df_taxa[df_taxa$taxa == cat,]$chisq = kw$statistic # return the chi square stat
  df_taxa[df_taxa$taxa == cat,]$pvals = kw$p.value # and the p-value
  
  if (kw$p.value <= 0.05) { # filter by significant p values
    listofcats_sig = c(cat, listofcats_sig)

  }
}

df_taxa$pvals.bh <- p.adjust(df_taxa$pvals, method="BH")
kable(df_taxa) # this tells us which ASVs vary by functional group

# now let's figure out which functional groups
pvals.dunn = NULL
pvals.dunn.bh = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$OTU == cat)
  new_df$ST <- as.factor(new_df$TreatmentName)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bh") # bh adjustment done here, no need to do it below
  
  for (i in 1:length(dT$res$Comparison)) {

    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bh = c(dT$res$P.adj[i], pvals.dunn.bh)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(as.character(dT$res$Comparison)[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)

  }

}

df_taxa <- data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bh)
kable(df_taxa, caption = "Significant Differences between Functional Groups")

#### Order (16S): G v F ####
## Is this code correct? double check, ITS uses tax_glom
AvgRA99_Order <- filter_taxa(ps.16s.nocontrols.rare.RA, function(x) mean(x) > .007, TRUE)
df_Order <- psmelt(AvgRA99_Order)
grouped_Order <- group_by(df_Order, FunGroup, Family, Order, Class, Phylum, Kingdom)
avgs_Order <- summarise(grouped_Order, mean=100*mean(Abundance), sd=100*sd(Abundance), se=100*se(Abundance))

ggplot(avgs_Order, aes(x = Family, y = mean, fill = Order)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se)), 
                width = .4, position = position_dodge(.9)) +
  facet_wrap(~FunGroup, nrow = 3) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 20)) + 
  ylab("Mean Relative Abundance")


listofcats_sig = NULL
DataSet = df_Order
taxa = unique(DataSet$Order)
df_taxa <- data.frame(taxa = taxa, chisq = rep(NA, length(taxa)), pvals = rep(NA, length(taxa)))

for (cat in df_taxa$taxa) { # for each Order that has a mean abundance above 1%
  new_df <- subset(DataSet, DataSet$Order == cat)
  new_df$ST <- as.factor(new_df$FunGroup)
  kw = kruskal.test(Abundance ~ ST, data=new_df) # test whether the abundance of that Order varies by ST, defined above as functional group
  df_taxa[df_taxa$taxa == cat,]$chisq = kw$statistic # return the chi square stat
  df_taxa[df_taxa$taxa == cat,]$pvals = kw$p.value # and the p-value
  
  if (kw$p.value <= 0.05) { # filter by significant p values
    listofcats_sig = c(cat, listofcats_sig)

  }
}

df_taxa$pvals.bh <- p.adjust(df_taxa$pvals, method="BH")
kable(df_taxa) # this tells us which Orders vary by functional group

# now let's figure out which functional groups
pvals.dunn = NULL
pvals.dunn.bh = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$Order == cat)
  new_df$ST <- as.factor(new_df$FunGroup)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bh") # bh adjustment done here, no need to do it below
  
  for (i in 1:length(dT$res$Comparison)) {

    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bh = c(dT$res$P.adj[i], pvals.dunn.bh)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(as.character(dT$res$Comparison)[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)

  }

}

df_taxa <- data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bh)
kable(df_taxa, caption = "Significant Differences between Functional Groups")

#### Order (16S): G v SA v ST ####
## Is this code correct? double check, ITS uses tax_glom
AvgRA99_Order <- filter_taxa(ps.16s.nocontrols.rare.RA, function(x) mean(x) > .007, TRUE)
df_Order <- psmelt(AvgRA99_Order)
grouped_Order <- group_by(df_Order, TreatmentName, Family, Order, Class, Phylum, Kingdom)
avgs_Order <- summarise(grouped_Order, mean=100*mean(Abundance), sd=100*sd(Abundance), se=100*se(Abundance))

ggplot(avgs_Order, aes(x = Family, y = mean, fill = Order)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se)), 
                width = .4, position = position_dodge(.9)) +
  facet_wrap(~TreatmentName, nrow = 3) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 20)) + 
  ylab("Mean Relative Abundance")


listofcats_sig = NULL
DataSet = df_Order
taxa = unique(DataSet$Order)
df_taxa <- data.frame(taxa = taxa, chisq = rep(NA, length(taxa)), pvals = rep(NA, length(taxa)))

for (cat in df_taxa$taxa) { # for each Order that has a mean abundance above 1%
  new_df <- subset(DataSet, DataSet$Order == cat)
  new_df$ST <- as.factor(new_df$TreatmentName)
  kw = kruskal.test(Abundance ~ ST, data=new_df) # test whether the abundance of that Order varies by ST, defined above as functional group
  df_taxa[df_taxa$taxa == cat,]$chisq = kw$statistic # return the chi square stat
  df_taxa[df_taxa$taxa == cat,]$pvals = kw$p.value # and the p-value
  
  if (kw$p.value <= 0.05) { # filter by significant p values
    listofcats_sig = c(cat, listofcats_sig)

  }
}

df_taxa$pvals.bh <- p.adjust(df_taxa$pvals, method="BH")
kable(df_taxa) # this tells us which Orders vary by functional group

# now let's figure out which functional groups
pvals.dunn = NULL
pvals.dunn.bh = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$Order == cat)
  new_df$ST <- as.factor(new_df$TreatmentName)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bh") # bh adjustment done here, no need to do it below
  
  for (i in 1:length(dT$res$Comparison)) {

    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bh = c(dT$res$P.adj[i], pvals.dunn.bh)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(as.character(dT$res$Comparison)[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)

  }

}

df_taxa <- data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bh)
kable(df_taxa, caption = "Significant Differences between Functional Groups")

#### Order (ITS):  G v F ####
AvgRA_g <- tax_glom(ps.ITS.nocontrols.rare.RA, taxrank = "Order", NArm = FALSE)
AvgRA99_Order <- filter_taxa(AvgRA_g, function(x) var(x) > .006, TRUE)
df_Order <- psmelt(AvgRA99_Order)
grouped_Order <- group_by(df_Order, FunGroup, Order, Class, Phylum, Kingdom)
avgs_Order <- summarise(grouped_Order, mean=100*mean(Abundance), sd=100*sd(Abundance), se=100*se(Abundance))

ggplot(avgs_Order, aes(x = Order, y = mean, fill = Phylum)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se)), 
                width = .4, position = position_dodge(.9)) +
  facet_wrap(~FunGroup, nrow = 3) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 20)) + 
  ylab("Mean Relative Abundance")


listofcats_sig = NULL
DataSet = df_Order
taxa = unique(DataSet$Order)
df_taxa <- data.frame(taxa = taxa, chisq = rep(NA, length(taxa)), pvals = rep(NA, length(taxa)))

for (cat in df_taxa$taxa) { # for each Order that has a mean abundance above 1%
  new_df <- subset(DataSet, DataSet$Order == cat)
  new_df$ST <- as.factor(new_df$FunGroup)
  kw = kruskal.test(Abundance ~ ST, data=new_df) # test whether the abundance of that Order varies by ST, defined above as functional group
  df_taxa[df_taxa$taxa == cat,]$chisq = kw$statistic # return the chi square stat
  df_taxa[df_taxa$taxa == cat,]$pvals = kw$p.value # and the p-value
  
  if (kw$p.value <= 0.05) { # filter by significant p values
    listofcats_sig = c(cat, listofcats_sig)

  }
}

df_taxa$pvals.bh <- p.adjust(df_taxa$pvals, method="BH")
kable(df_taxa) # this tells us which Orders vary by functional group

# now let's figure out which functional groups (give any vary above which they dont)
# pvals.dunn = NULL
# pvals.dunn.bh = NULL
# Zsc.dunn = NULL
# comparison.dunn = NULL
# cats.dunn = NULL
# 
# for (cat in listofcats_sig){
#   new_df <- subset(DataSet, DataSet$Order == cat)
#   new_df$ST <- as.factor(new_df$FunGroup)
#   dT = dunnTest(Abundance ~ ST, data = new_df, method = "bh") # bh adjustment done here, no need to do it below
#   
#   for (i in 1:length(dT$res$Comparison)) {
# 
#     pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
#     pvals.dunn.bh = c(dT$res$P.adj[i], pvals.dunn.bh)
#     Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
#     comparison.dunn = c(as.character(dT$res$Comparison)[i], comparison.dunn)
#     cats.dunn = c(cat, cats.dunn)
# 
#   }
# 
# }
# 
# df_taxa <- data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bh)
# kable(df_taxa, caption = "Significant Differences between Functional Groups")

#### Order (ITS): G v SA v ST ####
AvgRA_g <- tax_glom(ps.ITS.nocontrols.rare.RA, taxrank = "Order", NArm = FALSE)
AvgRA99_Order <- filter_taxa(AvgRA_g, function(x) var(x) > .006, TRUE)
df_Order <- psmelt(AvgRA99_Order)
grouped_Order <- group_by(df_Order, TreatmentName, Order, Class, Phylum, Kingdom)
avgs_Order <- summarise(grouped_Order, mean=100*mean(Abundance), sd=100*sd(Abundance), se=100*se(Abundance))

ggplot(avgs_Order, aes(x = Order, y = mean, fill = Phylum)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se)), 
                width = .4, position = position_dodge(.9)) +
  facet_wrap(~TreatmentName, nrow = 3) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 20)) + 
  ylab("Mean Relative Abundance")


listofcats_sig = NULL
DataSet = df_Order
taxa = unique(DataSet$Order)
df_taxa <- data.frame(taxa = taxa, chisq = rep(NA, length(taxa)), pvals = rep(NA, length(taxa)))

for (cat in df_taxa$taxa) { # for each Order that has a mean abundance above 1%
  new_df <- subset(DataSet, DataSet$Order == cat)
  new_df$ST <- as.factor(new_df$TreatmentName)
  kw = kruskal.test(Abundance ~ ST, data=new_df) # test whether the abundance of that Order varies by ST, defined above as functional group
  df_taxa[df_taxa$taxa == cat,]$chisq = kw$statistic # return the chi square stat
  df_taxa[df_taxa$taxa == cat,]$pvals = kw$p.value # and the p-value
  
  if (kw$p.value <= 0.05) { # filter by significant p values
    listofcats_sig = c(cat, listofcats_sig)

  }
}

df_taxa$pvals.bh <- p.adjust(df_taxa$pvals, method="BH")
kable(df_taxa) # this tells us which Orders vary by functional group

# now let's figure out which functional groups (give any vary above which they dont)
pvals.dunn = NULL
pvals.dunn.bh = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$Order == cat)
  new_df$ST <- as.factor(new_df$TreatmentName)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bh") # bh adjustment done here, no need to do it below

  for (i in 1:length(dT$res$Comparison)) {

    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bh = c(dT$res$P.adj[i], pvals.dunn.bh)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(as.character(dT$res$Comparison)[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)

  }

}

df_taxa <- data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bh)
kable(df_taxa, caption = "Significant Differences between Functional Groups")

#### Alpha diveristy: 16S ####
plot_richness(ps.16s.nocontrols.rare, measures = c("Observed", "Shannon"), x = "FunGroup", color = "FunGroup") + 
  theme(text = element_text(size=24)) + geom_boxplot() + 
  geom_jitter() +
  theme(axis.text.x = element_text(angle=-70, hjust=0, vjust=.5)) 

## Stats ##

GL_Alpha <- estimate_richness(ps.16s.nocontrols.rare, measures = c("Observed","Shannon", "InvSimpson"))
GL_Alpha2 <- cbind(GL_Alpha, sample_data(ps.16s.nocontrols.rare))

#kruskal_test(Observed ~ FunGroup, distribution = approximate(nresample = 9999), data = GL_Alpha2)

kruskal_test(Shannon ~ FunGroup, distribution = approximate(nresample = 9999), data = GL_Alpha2)
# breaking news: shannon diversity does vary

dunnTest(Shannon ~ FunGroup, data = GL_Alpha2, method = "bh") # forbs have higher alpha diversity than grasses; grasses with forbs have marginally higher diversity than grasses, but not than forbs

#### Alpha diversity: ITS ####
plot_richness(ps.ITS.nocontrols.rare, measures = c("Observed", "Shannon"), x = "FunGroup", color = "FunGroup") + 
  theme(text = element_text(size=24)) + geom_boxplot() + 
  geom_jitter() +
  theme(axis.text.x = element_text(angle=-70, hjust=0, vjust=.5)) 

## Stats ##

GL_Alpha <- estimate_richness(ps.ITS.nocontrols.rare, measures = c("Observed","Shannon", "InvSimpson"))
GL_Alpha2 <- cbind(GL_Alpha, sample_data(ps.ITS.nocontrols.rare))

kruskal_test(Observed ~ FunGroup, distribution = approximate(nresample = 9999), data = GL_Alpha2)

kruskal_test(Shannon ~ FunGroup, distribution = approximate(nresample = 9999), data = GL_Alpha2)
# breaking news: shannon diversity doesn't vary

#### Mantel test: PC v 16S ####
ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion)))

man.pc.16s <- subset_samples(ps.16s.nocontrols.rare, !(is.na(PC.F))) # this test includes forbs in competition

dist_mic <- vegdist(otu_table(man.pc.16s), method = "bray") # establishes Bray-Curtis dissimilarity matrix of 16S community

head(sample_data(man.pc.16s)) # Figure out which columns of your sample data contain the relevant plant data. You don't want zeros in the data if possible. 

dist_PC <- vegdist(sample_data(man.pc.16s)[,33], method = "euclidean") 

mantel(dist_mic, dist_PC, method = "spearman", permutations = 9999) #significant if p < 0.05 (or your p value of choice) # marginally significant


#### Mantel test: is comp comm more correlated with grass or with forb? ####
# ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion)))
# 
# man.g.16s <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "Grass")
# man.f.16s <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "Forb")
# man.c.16s <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "grass_x_forb")
# 
# dist_micg <- vegdist(otu_table(man.g.16s), method = "bray") 
# dist_micf <- vegdist(otu_table(man.f.16s), method = "bray") 
# dist_micc <- vegdist(otu_table(man.c.16s), method = "bray") 
# 
# mantel(dist_micg, dist_micc, method = "spearman") # doesnt work i assume because they are of different sizes... maybe we could sample from the objects so that they are even?

#### Mantel test: sing sp weight v 16S ####
ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion)))

man.b.16s <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Weight.g)))

dist_mic <- vegdist(otu_table(man.b.16s), method = "bray") # establishes Bray-Curtis dissimilarity matrix of 16S community

head(sample_data(man.b.16s))[,32] # Figure out which columns of your sample data contain the relevant plant data. You don't want zeros in the data if possible. 

dist_b <- vegdist(sample_data(man.b.16s)[,32], method = "euclidean") 

mantel(dist_mic, dist_b, method = "spearman", permutations = 9999) #significant if p < 0.05 (or your p value of choice) # weight in single species not correlated with community 

#### Mantel test: forb weight v 16S ####
ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion)))

man.f.16s <- subset_samples(ps.16s.nocontrols.rare, !(is.na(forb.wt)))

dist_micf <- vegdist(otu_table(man.f.16s), method = "bray") # establishes Bray-Curtis dissimilarity matrix of 16S community

head(sample_data(man.f.16s))[,30] # Figure out which columns of your sample data contain the relevant plant data. You don't want zeros in the data if possible. 

dist_b <- vegdist(sample_data(man.f.16s)[,30], method = "euclidean") 

mantel(dist_micf, dist_b, method = "spearman", permutations = 9999) #significant if p < 0.05 (or your p value of choice) # forb weight (alone or in comp) not correlated with community 

# honestly cant tell if these really make sense, can we say that plants with similar weights will have more similar communities? This would make more sense to do it with carbon:nitrogen ratio

#### Mantel test: grass weight v 16S ####
ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion)))

man.g.16s <- subset_samples(ps.16s.nocontrols.rare, !(is.na(grass.wt)))

dist_micg <- vegdist(otu_table(man.g.16s), method = "bray") # establishes Bray-Curtis dissimilarity matrix of 16S community

head(sample_data(man.g.16s))[,31] # Figure out which columns of your sample data contain the relevant plant data. You don't want zeros in the data if possible. 

dist_b <- vegdist(sample_data(man.g.16s)[,31], method = "euclidean") 

mantel(dist_micg, dist_b, method = "spearman", permutations = 9999) #significant if p < 0.05 (or your p value of choice) # grass weight (alone or in comp) IS correlated with community

#### Mantel test: Biomass v ITS ####
# decided not to do this because I'm not sure it really makes sense

#### Mantel test: 16S v ITS ####
man.ITS <- subset_samples(ps.ITS.nocontrols.rare, SampleID_Fix != "GLM-0116")
man.16s <- ps.16s.nocontrols.rare

dist_ITS <- vegdist(otu_table(man.ITS), method = "bray") 
dist_16s <- vegdist(otu_table(man.16s), method = "bray") 

mantel(dist_ITS, dist_16s, method = "pearson", permutations = 9999) # non parametric, they are highly correlated

# do a (generalized) linear model linear model
m.man <- lm(dist_ITS ~ dist_16s)
summary(m.man)

#### Core Functions Needed ####
fo_difference <- function(pos){
  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
  return(left - right)
}

sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
  require(minpack.lm)
  require(Hmisc)
  require(stats4)
  
  options(warn=-1)
  
  #Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))
  
  #Calculate the average relative abundance of each taxa across communities
  if(is.null(pool)){
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  }
  
  #Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  #Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  #Calculate the limit of detection
  d = 1/N
  
  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
  sncm.LL <- function(m, sigma){
    R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
    R = dnorm(R, 0, sigma)
    -sum(log(R))
  }
  m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
  
  ##Calculate Akaike's Information Criterion (AIC)
  aic.fit <- AIC(m.mle, k=2)
  bic.fit <- BIC(m.mle)
  
  ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for binomial model
  bino.LL <- function(mu, sigma){
    R = freq - pbinom(d, N, p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.bino <- AIC(bino.mle, k=2)
  bic.bino <- BIC(bino.mle)
  
  ##Goodness of fit for binomial model
  bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
  Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
  
  bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)
  
  ##Goodness of fit for Poisson model
  pois.pred <- ppois(d, N*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
  
  pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Results
  if(stats==TRUE){
    fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
    fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
    return(fitstats)
  } else {
    A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
    A <- as.data.frame(A)
    colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
    if(is.null(taxon)){
      B <- A[order(A[,1]),]
    } else {
      B <- merge(A, taxon, by=0, all=TRUE)
      row.names(B) <- B[,1]
      B <- B[,-1]
      B <- B[order(B[,1]),]
    }
    return(B)
  }
}

#### Core ASV (16s): All ####

#NOTE: Issues arise when all the packages are loaded, only packages needed are tidyverse, reshape2, and vegan

ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion))) # these were likely removed earlier but just to be sure
test <- psmelt(ps.16s.nocontrols.rare)
otu <- as.data.frame(t(otu_table(ps.16s.nocontrols.rare)))

nReads <- 9434 
map <- mapping

otu_PA <- 1*((otu>0)==1)  # presence-absence data (if OTU is present (greater than 0) assign a 1)
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  # total number of sites that OTU is present in, divided by the number of sites  (occupancy calculation)
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     # relative abundance: For each column divide every entry by the column total (this give relative abundance  of each OTU per sample), then give calculate the mean relative abundance per OTU by calculating the mean of each row 
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance; occupancy of an OTU is average number of samples an OTU occurs in, abundance is the average relative abundance of that OTU within a sample

# for some reason this code replaces the - with a . in our sample names, so creating a new column to make it run (this is used later in PlotDF)
map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)
# 
# PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
#   gather(SampleID.occ, abun, -otu) %>%
#   left_join(map, by = 'SampleID.occ') %>%
#   group_by(otu, FunGroup) %>% # first try it by FunGroup?
#   summarise(plot_freq = sum(abun > 0)/length(abun), # number of samples within a treatment where that OTU was present (greater than 0) divided by total number of samples within that treatment; so this is 1 if the OTU was present in every sample in that treatment and lower otherwise; so basically it's the percentage of subsamples within a treatment that OTU occurs in
#             coreTrt = ifelse(plot_freq == 1, 1, 0), # Core Treatment = 1 if that OTU was present in every sample within a Treatment, 0 otherwise
#             detect = ifelse(plot_freq > 0, 1, 0)) %>%    # 1 if that OTU was detected at all within that Treatment and 0 if not
#   group_by(otu) %>% # group by OTU
#   summarise(sumF = sum(plot_freq), # frequency of an OTU across treatments, so with 3 treatments (grass, forb, and grassxforb), a 3 would indicate that this OTU occurs in every treatment and every subsample of that treatment
#             sumG = sum(coreTrt), # total number of Treatments where that OTU was present in every subsample
#             nS = length(FunGroup)*2, # total number of Treatments (times 2) = 6
#             Index = (sumF + sumG)/nS) # calculating weighting Index based on number of Treatments detected
# 
# otu_ranked2 <- occ_abun %>%
#   left_join(PresenceSum, by = 'otu') %>% # why even need occ_abun here?
#   transmute(otu = otu,
#             rank = Index) %>%
#   arrange(desc(rank))

otu_ranked <- occ_abun %>%
  #left_join(PresenceSum, by='otu') %>%
  transmute(otu=otu,
            rank = (otu_rel) + (otu_occ)) %>%
  arrange(desc(rank))

BCaddition <- NULL

otu_start <- otu_ranked$otu[1] # take the first ranked OTU
start_matrix <- as.matrix(otu[otu_start,]) # extract that OTU's abundance per sample, this should be a one column vector
#start_matrix <- t(start_matrix) # turns it into a one row vector, but mine already is because my original OTU table is a dataframe, theirs is already a matrix

x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) #Error in combn(ncol(start_matrix), 2) : n < m (fixed above); take every combination of samples, and for each combination give the absolute difference in those two OTU abundances, sum them all (sum what?) and divide by 2* number of reads (because 2 samples); i have no idea what the sum function is doing given that inside sum there's only one number, this is relevant later I believe?
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 # i dont fully understand what this column is 
BCaddition <- rbind(BCaddition,df_s)
  
for(i in 2:500){
  otu_add = otu_ranked$otu[i] # for the top 500 ranked OTUs, starting with the second one beacuse we computed the first previously
  add_matrix <- as.matrix(otu[otu_add,])
  #add_matrix <- t(add_matrix) # again dont need this because its a dataframe
  start_matrix <- rbind(start_matrix, add_matrix) # bind together the next ranked OTU with the previously ranked OTU
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads)) # and compute the difference again
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
  df_a <- data.frame(x_names, x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))  
}

# and this code does it for the last one, I think 
x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]] - otu[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse = ' - '))
df_full <- data.frame(x_names, x)
#names(df_full)[2] <- length(rownames(otu)) # this assumes we are looking at all OTUs, not just the top 500 ranked, right?
names(df_full)[2] <- i + 1 #I think
BCfull <- left_join(BCaddition, df_full, by = 'x_names') # each column is an OTU, each row a comparison of two samples, so the first (non-sample) column is the comparison of the first ranked OTU's abundance between the two samples specified, the second column is the difference in the second rank OTU's abundance between teh two samples specified and so on all the way until the 500th ranked OTU; it does this for each combination of samples
 
rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC = mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC = MeanBC/max(MeanBC))

Increase = BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase) + 1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)

BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

(elbow.all.p <- ggplot(BC_ranked[1:250,], aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
  geom_point(aes(y = proportionBC)) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_text(size = 7, angle = 45)) +
  geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
  geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])), 
             lty = 3, col = 'blue', cex = .5) +
  labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
  annotate(geom = "text", 
           x = elbow + 10, 
           y = .15, 
           label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +    
  annotate(geom = "text", 
           x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4, 
           y = .08, 
           label = paste("Last 2% increase (",last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))

occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))]] <- 'core'

# Models for the whole community
obs.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats = FALSE, pool = NULL)
#sta.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats = TRUE, pool = NULL)

(core.all.p <- ggplot() +
  geom_point(data = occ_abun[occ_abun$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun[occ_abun$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: All Samples", x = "log10(mean relative abundance)", y = "Occupancy"))

core.all <- occ_abun$otu[occ_abun$fill == 'core']

otu_relabun <- decostand(otu, method = "total", MARGIN = 2)

plotDF <-  data.frame(otu = as.factor(row.names(otu_relabun)), otu_relabun) %>% # this is changing my sample names for some reason, replacing the "-" with a "."
  gather(SampleID.occ, relabun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  left_join(otu_ranked, by = 'otu') %>%
  filter(otu %in% core.all) %>% 
  group_by(FunGroup, otu) %>%
  summarise(plot_freq = sum(relabun>0)/length(relabun), 
            coreSite = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0))

plotDF$otu <- factor(plotDF$otu, levels = otu_ranked$otu[1:213]) # 1: # OTUs before cutoff of last 2% increase, but there are only 96 in the core so why up to 205? does this just rank them

(all.plotDF <- ggplot(plotDF, aes(x = otu, plot_freq, group = FunGroup, fill = FunGroup)) +    
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(plotDF$otu %in% core.all))) +
  theme(axis.text = element_text(size = 6)) +
  labs(x = 'Ranked OTUs', y = 'Occupancy by site'))

#### Core ASV (16s): Grass ####
grass.core <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "Grass")

otu <- as.data.frame(t(otu_table(grass.core)))

nReads <- 9434 

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 

# for some reason this code replaces the - with a . in our sample names, so creating a new column to make it run (this is used later in PlotDF)
map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

otu_ranked <- occ_abun %>%
  transmute(otu=otu,
            rank = (otu_rel) + (otu_occ)) %>%
  arrange(desc(rank))

BCaddition <- NULL

otu_start <- otu_ranked$otu[1] 
start_matrix <- as.matrix(otu[otu_start,])

x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) 
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
BCaddition <- rbind(BCaddition,df_s)
  
for(i in 2:500){
  otu_add=otu_ranked$otu[i] 
  add_matrix <- as.matrix(otu[otu_add,])
  start_matrix <- rbind(start_matrix, add_matrix) 
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads)) 
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))  
}


x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names, x)
names(df_full)[2] <- i + 1 
BCfull <- left_join(BCaddition,df_full, by = 'x_names') 
 
rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))

Increase <- BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC = c(0, (Increase)), rank = factor(c(1:(length(Increase) + 1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)

BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

(elbow.g.p <- ggplot(BC_ranked[1:350,], 
                       aes(x = factor(BC_ranked$rank[1:350], levels = BC_ranked$rank[1:450]))) +
  geom_point(aes(y = proportionBC)) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_text(size = 7, angle = 45)) +
  geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
  geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])), 
             lty = 3, col = 'blue', cex = .5) +
  labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
  annotate(geom = "text", 
           x = elbow + 10, 
           y = .15, 
           label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +    
  annotate(geom = "text", 
           x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4, 
           y = .08, 
           label = paste("Last 2% increase (", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))

occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)]))]] <- 'core'

obs.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats=FALSE, pool = NULL)
sta.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats = TRUE, pool = NULL)

(core.g.p <- ggplot() +
  geom_point(data = occ_abun[occ_abun$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun[occ_abun$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: All Samples", x = "log10(mean relative abundance)", y = "Occupancy"))

core.g <- occ_abun$otu[occ_abun$fill == 'core']

otu_relabun <- decostand(otu, method = "total", MARGIN = 2)

#### Core ASV (16s): Forb ####
forb.core <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "Forb")

otu <- as.data.frame(t(otu_table(forb.core)))

nReads <- 9434 

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 

# for some reason this code replaces the - with a . in our sample names, so creating a new column to make it run (this is used later in PlotDF)
map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

otu_ranked <- occ_abun %>%
  transmute(otu=otu,
            rank = (otu_rel) + (otu_occ)) %>%
  arrange(desc(rank))

BCaddition <- NULL

otu_start <- otu_ranked$otu[1] 
start_matrix <- as.matrix(otu[otu_start,])

x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) 
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
BCaddition <- rbind(BCaddition,df_s)
  
for(i in 2:500){
  otu_add=otu_ranked$otu[i] 
  add_matrix <- as.matrix(otu[otu_add,])
  start_matrix <- rbind(start_matrix, add_matrix) 
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads)) 
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))  
}


x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names, x)
names(df_full)[2] <- i + 1 
BCfull <- left_join(BCaddition,df_full, by = 'x_names') 
 
rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))

Increase <- BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase)+1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)

BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

(elbow.f.p <- ggplot(BC_ranked[1:250,], 
                       aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
  geom_point(aes(y = proportionBC)) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_text(size = 7, angle = 45)) +
  geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
  geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])), 
             lty = 3, col = 'blue', cex = .5) +
  labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
  annotate(geom = "text", 
           x = elbow + 50, 
           y = .15, 
           label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +    
  annotate(geom = "text", 
           x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 50, 
           y = .08, 
           label = paste("Last 2% increase (", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))

occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)]))]] <- 'core'

obs.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats = FALSE, pool = NULL)
sta.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats = TRUE, pool = NULL)

(core.f.p <- ggplot() +
  geom_point(data = occ_abun[occ_abun$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun[occ_abun$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: All Samples", x = "log10(mean relative abundance)", y = "Occupancy"))

core.f <- occ_abun$otu[occ_abun$fill == 'core']

otu_relabun <- decostand(otu, method = "total", MARGIN = 2)

#### Core ASV (16s): Grass x Forb ####
gf.core <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "grass_x_forb")

otu <- as.data.frame(t(otu_table(gf.core)))

nReads <- 9434 

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 

# for some reason this code replaces the - with a . in our sample names, so creating a new column to make it run (this is used later in PlotDF)
map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

otu_ranked <- occ_abun %>%
  transmute(otu=otu,
            rank = (otu_rel) + (otu_occ)) %>%
  arrange(desc(rank))

BCaddition <- NULL

otu_start <- otu_ranked$otu[1] 
start_matrix <- as.matrix(otu[otu_start,])

x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) 
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
BCaddition <- rbind(BCaddition,df_s)
  
for(i in 2:500){
  otu_add = otu_ranked$otu[i] 
  add_matrix <- as.matrix(otu[otu_add,])
  start_matrix <- rbind(start_matrix, add_matrix) 
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads)) 
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))  
}


x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names, x)
names(df_full)[2] <- i + 1 
BCfull <- left_join(BCaddition, df_full, by = 'x_names') 
 
rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC = mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC = MeanBC/max(MeanBC))

Increase = BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase) + 1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)

BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

(elbow.gf.p <- ggplot(BC_ranked[1:250,], 
                       aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
  geom_point(aes(y = proportionBC)) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_text(size = 7, angle = 45)) +
  geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
  geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])), 
             lty = 3, col = 'blue', cex = .5) +
  labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
  annotate(geom = "text", 
           x = elbow + 10, 
           y = .15, 
           label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +    
  annotate(geom = "text", 
           x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4, 
           y = .08, 
           label = paste("Last 2% increase (", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))

occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)]))]] <- 'core'

obs.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats=FALSE, pool = NULL)
sta.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats = TRUE, pool = NULL)

(core.gf.p <- ggplot() +
  geom_point(data = occ_abun[occ_abun$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun[occ_abun$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: All Samples", x = "log10(mean relative abundance)", y = "Occupancy"))

core.gf <- occ_abun$otu[occ_abun$fill == 'core']

otu_relabun <- decostand(otu, method = "total", MARGIN = 2)

#### Comparison of Core ####
library(VennDiagram)
core.g <- data.frame(FunGroup = "Grass", OTU = core.g)
core.f <- data.frame(FunGroup = "Forb", OTU = core.f)
core.gf <- data.frame(FunGroup = "grass_x_forb", OTU = core.gf)
core.all <- data.frame(FunGroup = "All", OTU = core.all)

# forb = 1
# grass = 2
# grass x forb = 3

nrow(core.f[core.f$OTU %in% core.g$OTU,]) # n12: 206 forb OTUs also in grass
nrow(core.g[core.g$OTU %in% core.gf$OTU,]) # n23: 201 grass OTUs in forbsxgrasses
nrow(core.f[core.f$OTU %in% core.gf$OTU,]) # n13: 178 forb OTUs in forbsxgrasses

test <- core.g[core.g$OTU %in% core.gf$OTU,]
nrow(core.f[core.f$OTU %in% test$OTU,]) # n123

# core <- merge(core.gf, core.g, by = "OTU", all = T)
# core <- merge(core, core.f, by = "OTU", all = T)
# core <- merge(core, core.all, by = "OTU", all = T)

# grid.newpage()  
# draw.triple.venn(area1 = 241,   # Forb core                  
#                  area2 = 309,   # Grass core
#                  area3 = 219,   # g_x_f core
#                  n12 = 206,
#                  n23 = 201,
#                  n13 = 178,
#                  n123 = 168,
#                  fill = c("pink", "green", "orange"),
#                  lty = "blank",
#                  category = c("Forbs", "Grasses", "Forbs & Grasses"))

grid.newpage()  
draw.triple.venn(area1 = 219,   # Forb core                  
                 area2 = 267,   # Grass core
                 area3 = 212,   # g_x_f core
                 n12 = 180,
                 n23 = 188,
                 n13 = 168,
                 n123 = 154,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Forbs", "Grasses", "Forbs & Grasses"))
# grasses have a bigger core, but forbs are more diverse

#### Core Family (16s): All ####

### CANT GET THIS TO WORK AS ON APRIL 29 ###
#NOTE: Issues arise when all the packages are loaded, only packages needed are tidyverse, reshape2, and vegan

ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion))) # these were likely removed earlier but just to be sure
core.all.ord <- tax_glom(ps.16s.nocontrols.rare, taxrank = "Family", NArm = FALSE)

test <- data.frame(otu = rownames(t(otu_table(core.all.ord))), t(otu_table(core.all.ord)))
test2 <- t(data.frame(t(tax_table(core.all.ord))))
test2 <- data.frame(otu = rownames(test2), test2)

# test2[test2$Family == "Unknown_Family" & test2$Phylum == "Acidobacteria",]$Family <- "Unk_p_Acido_c_Sub6"
# test2[test2$Family == "Unknown_Family" & test2$Order == "Oxyphotobacteria_Incertae_Sedis",]$Family <- "Unk_c_Oxyphotobacteria"
# test2[test2$Family == "Unknown_Family" & test2$Order == "Gammaproteobacteria",]$Family <- "Unk_c_Gammaproteobacteria
# "
tmp <- merge(test2[,c(1,6)], test, by = "otu")
tmp <- filter(tmp, Family != "Unknown_Family", Family != "Family_XI")
rownames(tmp) <- tmp$Family
otu <- tmp[ , c(3:128)]

nReads <- 9434 
map <- mapping

otu_PA <- 1*((otu>0)==1)  # presence-absence data (if OTU is present (greater than 0) assign a 1)
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  # total number of sites that OTU is present in, divided by the number of sites  (occupancy calculation)
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     # relative abundance: For each column divide every entry by the column total (this give relative abundance  of each OTU per sample), then give calculate the mean relative abundance per OTU by calculating the mean of each row 
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance; occupancy of an OTU is average number of samples an OTU occurs in, abundance is the average relative abundance of that OTU within a sample

# for some reason this code replaces the - with a . in our sample names, so creating a new column to make it run (this is used later in PlotDF)
map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)
# 
# PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
#   gather(SampleID.occ, abun, -otu) %>%
#   left_join(map, by = 'SampleID.occ') %>%
#   group_by(otu, FunGroup) %>% # first try it by FunGroup?
#   summarise(plot_freq = sum(abun > 0)/length(abun), # number of samples within a treatment where that OTU was present (greater than 0) divided by total number of samples within that treatment; so this is 1 if the OTU was present in every sample in that treatment and lower otherwise; so basically it's the percentage of subsamples within a treatment that OTU occurs in
#             coreTrt = ifelse(plot_freq == 1, 1, 0), # Core Treatment = 1 if that OTU was present in every sample within a Treatment, 0 otherwise
#             detect = ifelse(plot_freq > 0, 1, 0)) %>%    # 1 if that OTU was detected at all within that Treatment and 0 if not
#   group_by(otu) %>% # group by OTU
#   summarise(sumF = sum(plot_freq), # frequency of an OTU across treatments, so with 3 treatments (grass, forb, and grassxforb), a 3 would indicate that this OTU occurs in every treatment and every subsample of that treatment
#             sumG = sum(coreTrt), # total number of Treatments where that OTU was present in every subsample
#             nS = length(FunGroup)*2, # total number of Treatments (times 2) = 6
#             Index = (sumF + sumG)/nS) # calculating weighting Index based on number of Treatments detected
# 
# otu_ranked2 <- occ_abun %>%
#   left_join(PresenceSum, by = 'otu') %>% # why even need occ_abun here?
#   transmute(otu = otu,
#             rank = Index) %>%
#   arrange(desc(rank))

otu_ranked <- occ_abun %>%
  #left_join(PresenceSum, by='otu') %>%
  transmute(otu=otu,
            rank = (otu_rel) + (otu_occ)) %>%
  arrange(desc(rank))

BCaddition <- NULL

otu_start <- otu_ranked$otu[1] # take the first ranked OTU
start_matrix <- as.matrix(otu[otu_start,]) # extract that OTU's abundance per sample, this should be a one column vector
#start_matrix <- t(start_matrix) # turns it into a one row vector, but mine already is because my original OTU table is a dataframe, theirs is already a matrix

x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) #Error in combn(ncol(start_matrix), 2) : n < m (fixed above); take every combination of samples, and for each combination give the absolute difference in those two OTU abundances, sum them all (sum what?) and divide by 2* number of reads (because 2 samples); i have no idea what the sum function is doing given that inside sum there's only one number, this is relevant later I believe?
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 # i dont fully understand what this column is 
BCaddition <- rbind(BCaddition,df_s)
  
for(i in 2:294){
  otu_add = otu_ranked$otu[i] # for the top 500 ranked OTUs, starting with the second one beacuse we computed the first previously
  add_matrix <- as.matrix(otu[otu_add,])
  #add_matrix <- t(add_matrix) # again dont need this because its a dataframe
  start_matrix <- rbind(start_matrix, add_matrix) # bind together the next ranked OTU with the previously ranked OTU
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads)) # and compute the difference again
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
  df_a <- data.frame(x_names, x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))  
}

# and this code does it for the last one, I think 
x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]] - otu[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse = ' - '))
df_full <- data.frame(x_names, x)
#names(df_full)[2] <- length(rownames(otu)) # this assumes we are looking at all OTUs, not just the top 500 ranked, right?
names(df_full)[2] <- i + 1 #I think
BCfull <- left_join(BCaddition, df_full, by = 'x_names') # each column is an OTU, each row a comparison of two samples, so the first (non-sample) column is the comparison of the first ranked OTU's abundance between the two samples specified, the second column is the difference in the second rank OTU's abundance between teh two samples specified and so on all the way until the 500th ranked OTU; it does this for each combination of samples
 
rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC = mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC = MeanBC/max(MeanBC))

Increase = BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase) + 1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)

BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

(elbow.all.p <- ggplot(BC_ranked[1:250,], aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
  geom_point(aes(y = proportionBC)) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_text(size = 7, angle = 45)) +
  geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
  geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])), 
             lty = 3, col = 'blue', cex = .5) +
  labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
  annotate(geom = "text", 
           x = elbow + 10, 
           y = .15, 
           label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +    
  annotate(geom = "text", 
           x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4, 
           y = .08, 
           label = paste("Last 2% increase (",last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))

occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))]] <- 'core'

# Models for the whole community
obs.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats = FALSE, pool = NULL)
#sta.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats = TRUE, pool = NULL)

(core.all.p <- ggplot() +
  geom_point(data = occ_abun[occ_abun$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun[occ_abun$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: All Samples", x = "log10(mean relative abundance)", y = "Occupancy"))

core.all <- occ_abun$otu[occ_abun$fill == 'core']

otu_relabun <- decostand(otu, method = "total", MARGIN = 2)

plotDF <-  data.frame(otu = as.factor(row.names(otu_relabun)), otu_relabun) %>% # this is changing my sample names for some reason, replacing the "-" with a "."
  gather(SampleID.occ, relabun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  left_join(otu_ranked, by = 'otu') %>%
  filter(otu %in% core.all) %>% 
  group_by(FunGroup, otu) %>%
  summarise(plot_freq = sum(relabun>0)/length(relabun), 
            coreSite = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0))

plotDF$otu <- factor(plotDF$otu, levels = otu_ranked$otu[1:213]) # 1: # OTUs before cutoff of last 2% increase, but there are only 96 in the core so why up to 205? does this just rank them

(all.plotDF <- ggplot(plotDF, aes(x = otu, plot_freq, group = FunGroup, fill = FunGroup)) +    
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(plotDF$otu %in% core.all))) +
  theme(axis.text = element_text(size = 6)) +
  labs(x = 'Ranked OTUs', y = 'Occupancy by site'))

#### Core Family (16s): Grass ####
grass.core <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "Grass")

otu <- as.data.frame(t(otu_table(grass.core)))

nReads <- 9434 

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 

# for some reason this code replaces the - with a . in our sample names, so creating a new column to make it run (this is used later in PlotDF)
map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

otu_ranked <- occ_abun %>%
  transmute(otu=otu,
            rank = (otu_rel) + (otu_occ)) %>%
  arrange(desc(rank))

BCaddition <- NULL

otu_start <- otu_ranked$otu[1] 
start_matrix <- as.matrix(otu[otu_start,])

x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) 
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
BCaddition <- rbind(BCaddition,df_s)
  
for(i in 2:500){
  otu_add=otu_ranked$otu[i] 
  add_matrix <- as.matrix(otu[otu_add,])
  start_matrix <- rbind(start_matrix, add_matrix) 
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads)) 
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))  
}


x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names, x)
names(df_full)[2] <- i + 1 
BCfull <- left_join(BCaddition,df_full, by = 'x_names') 
 
rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))

Increase <- BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC = c(0, (Increase)), rank = factor(c(1:(length(Increase) + 1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)

BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

(elbow.g.p <- ggplot(BC_ranked[1:350,], 
                       aes(x = factor(BC_ranked$rank[1:350], levels = BC_ranked$rank[1:450]))) +
  geom_point(aes(y = proportionBC)) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_text(size = 7, angle = 45)) +
  geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
  geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])), 
             lty = 3, col = 'blue', cex = .5) +
  labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
  annotate(geom = "text", 
           x = elbow + 10, 
           y = .15, 
           label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +    
  annotate(geom = "text", 
           x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4, 
           y = .08, 
           label = paste("Last 2% increase (", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))

occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)]))]] <- 'core'

obs.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats=FALSE, pool = NULL)
sta.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats = TRUE, pool = NULL)

(core.g.p <- ggplot() +
  geom_point(data = occ_abun[occ_abun$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun[occ_abun$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: All Samples", x = "log10(mean relative abundance)", y = "Occupancy"))

core.g <- occ_abun$otu[occ_abun$fill == 'core']

otu_relabun <- decostand(otu, method = "total", MARGIN = 2)

#### Core Family (16s): Forb ####
forb.core <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "Forb")

otu <- as.data.frame(t(otu_table(forb.core)))

nReads <- 9434 

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 

# for some reason this code replaces the - with a . in our sample names, so creating a new column to make it run (this is used later in PlotDF)
map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

otu_ranked <- occ_abun %>%
  transmute(otu=otu,
            rank = (otu_rel) + (otu_occ)) %>%
  arrange(desc(rank))

BCaddition <- NULL

otu_start <- otu_ranked$otu[1] 
start_matrix <- as.matrix(otu[otu_start,])

x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) 
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
BCaddition <- rbind(BCaddition,df_s)
  
for(i in 2:500){
  otu_add=otu_ranked$otu[i] 
  add_matrix <- as.matrix(otu[otu_add,])
  start_matrix <- rbind(start_matrix, add_matrix) 
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads)) 
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))  
}


x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names, x)
names(df_full)[2] <- i + 1 
BCfull <- left_join(BCaddition,df_full, by = 'x_names') 
 
rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))

Increase <- BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase)+1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)

BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

(elbow.f.p <- ggplot(BC_ranked[1:250,], 
                       aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
  geom_point(aes(y = proportionBC)) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_text(size = 7, angle = 45)) +
  geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
  geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])), 
             lty = 3, col = 'blue', cex = .5) +
  labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
  annotate(geom = "text", 
           x = elbow + 50, 
           y = .15, 
           label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +    
  annotate(geom = "text", 
           x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 50, 
           y = .08, 
           label = paste("Last 2% increase (", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))

occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)]))]] <- 'core'

obs.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats = FALSE, pool = NULL)
sta.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats = TRUE, pool = NULL)

(core.f.p <- ggplot() +
  geom_point(data = occ_abun[occ_abun$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun[occ_abun$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: All Samples", x = "log10(mean relative abundance)", y = "Occupancy"))

core.f <- occ_abun$otu[occ_abun$fill == 'core']

otu_relabun <- decostand(otu, method = "total", MARGIN = 2)

#### Core Family (16s): Grass x Forb ####
gf.core <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "grass_x_forb")

otu <- as.data.frame(t(otu_table(gf.core)))

nReads <- 9434 

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 

# for some reason this code replaces the - with a . in our sample names, so creating a new column to make it run (this is used later in PlotDF)
map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

otu_ranked <- occ_abun %>%
  transmute(otu=otu,
            rank = (otu_rel) + (otu_occ)) %>%
  arrange(desc(rank))

BCaddition <- NULL

otu_start <- otu_ranked$otu[1] 
start_matrix <- as.matrix(otu[otu_start,])

x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) 
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
BCaddition <- rbind(BCaddition,df_s)
  
for(i in 2:500){
  otu_add = otu_ranked$otu[i] 
  add_matrix <- as.matrix(otu[otu_add,])
  start_matrix <- rbind(start_matrix, add_matrix) 
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads)) 
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))  
}


x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names, x)
names(df_full)[2] <- i + 1 
BCfull <- left_join(BCaddition, df_full, by = 'x_names') 
 
rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC = mean(BC)) %>%
  arrange(-desc(MeanBC)) %>%
  mutate(proportionBC = MeanBC/max(MeanBC))

Increase = BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase) + 1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)

BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

(elbow.gf.p <- ggplot(BC_ranked[1:250,], 
                       aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
  geom_point(aes(y = proportionBC)) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_text(size = 7, angle = 45)) +
  geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
  geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])), 
             lty = 3, col = 'blue', cex = .5) +
  labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
  annotate(geom = "text", 
           x = elbow + 10, 
           y = .15, 
           label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +    
  annotate(geom = "text", 
           x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4, 
           y = .08, 
           label = paste("Last 2% increase (", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))

occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)]))]] <- 'core'

obs.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats=FALSE, pool = NULL)
sta.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats = TRUE, pool = NULL)

(core.gf.p <- ggplot() +
  geom_point(data = occ_abun[occ_abun$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun[occ_abun$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: All Samples", x = "log10(mean relative abundance)", y = "Occupancy"))

core.gf <- occ_abun$otu[occ_abun$fill == 'core']

otu_relabun <- decostand(otu, method = "total", MARGIN = 2)

#### Comparison of Core ####
library(VennDiagram)
core.g <- data.frame(FunGroup = "Grass", OTU = core.g)
core.f <- data.frame(FunGroup = "Forb", OTU = core.f)
core.gf <- data.frame(FunGroup = "grass_x_forb", OTU = core.gf)
core.all <- data.frame(FunGroup = "All", OTU = core.all)

# forb = 1
# grass = 2
# grass x forb = 3

nrow(core.f[core.f$OTU %in% core.g$OTU,]) # n12: 206 forb OTUs also in grass
nrow(core.g[core.g$OTU %in% core.gf$OTU,]) # n23: 201 grass OTUs in forbsxgrasses
nrow(core.f[core.f$OTU %in% core.gf$OTU,]) # n13: 178 forb OTUs in forbsxgrasses

test <- core.g[core.g$OTU %in% core.gf$OTU,]
nrow(core.f[core.f$OTU %in% test$OTU,]) # n123

# core <- merge(core.gf, core.g, by = "OTU", all = T)
# core <- merge(core, core.f, by = "OTU", all = T)
# core <- merge(core, core.all, by = "OTU", all = T)

# grid.newpage()  
# draw.triple.venn(area1 = 241,   # Forb core                  
#                  area2 = 309,   # Grass core
#                  area3 = 219,   # g_x_f core
#                  n12 = 206,
#                  n23 = 201,
#                  n13 = 178,
#                  n123 = 168,
#                  fill = c("pink", "green", "orange"),
#                  lty = "blank",
#                  category = c("Forbs", "Grasses", "Forbs & Grasses"))

grid.newpage()  
draw.triple.venn(area1 = 219,   # Forb core                  
                 area2 = 267,   # Grass core
                 area3 = 212,   # g_x_f core
                 n12 = 180,
                 n23 = 188,
                 n13 = 168,
                 n123 = 154,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Forbs", "Grasses", "Forbs & Grasses"))
# grasses have a bigger core, but forbs are more diverse
#### Network co-occurence ####

















#### Differential abundance analysis: 16s ####
## 16S: Soil v Rhiz  ##

treat.16s = phyloseq_to_deseq2(ps.16s.nocontrols, ~ SampleSubType)

dds.16s.treat = DESeq(treat.16s, test = "Wald", fitType = "parametric")

resultsNames(dds.16s.treat) #gives us comparisons for our contrasts, if needed

alpha = 0.01

#get results 
res <- results(dds.16s.treat, pAdjustMethod = "bonferroni")

#filter results by p-value
res.alpha <- res[which(res$padj < alpha), ]

#Bind taxonomy to results
res.alpha.tax = cbind(as(res.alpha, "data.frame"), as(tax_table(ps.16s.nocontrols)[rownames(res.alpha), ], "matrix"))

#tidy results 
res <- tidy(res.alpha)
  
#generate plot of significant ASVs for each contrast
# Order order
x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Order, function(x) max(x))
x = sort(x, TRUE)
res.alpha.tax$Order = factor(as.character(res.alpha.tax$Order), levels=names(x))
# Genus order
x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Genus, function(x) max(x))
x = sort(x, TRUE)
res.alpha.tax$Genus = factor(as.character(res.alpha.tax$Genus), levels=names(x))
p <- ggplot(res.alpha.tax, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("Plant Rhizosphere vs. Background soil")
  
#plot results

p

#tidy results into a table to save
df.res <- as.data.frame(res)
names(df.res)[1] <- "ASV"

write.csv(df.res, "Data/16s.ddseq.soilrhiz.csv")


## 16S:  G v SA v ST  ##
ps.16s.nocontrols.nobss <- subset_samples(ps.16s.nocontrols, SampleSubType !="Background_soil_sand_mix")

treat.16s = phyloseq_to_deseq2(ps.16s.nocontrols.nobss, ~ TreatmentName)

dds.16s.treat = DESeq(treat.16s, test = "Wald", fitType = "parametric")

resultsNames(dds.16s.treat) #gives us comparisons for our contrasts

contrasts = c("SA.G", "ST.G",
              "SAG.G", "STG.G")

contrast.list <- list(SA.G = c("TreatmentName", "Stress_avoiding_forb", "Invasive_grass"),
                      ST.G = c("TreatmentName", "Stress_tolerant_forb", "Invasive_grass"),
                      SAG.G = c("TreatmentName", "SA_forb_X_grass", "Invasive_grass"),
                      STG.G = c("TreatmentName", "ST_forb_X_grass", "Invasive_grass"))


plot.name.list <- list(SA.G = "SA forb vs. Grass",
                       ST.G = "ST forb vs. Grass",
                       SAG.G = "SA forb X grass vs. Grass",
                       STG.G = "ST forb X grass vs. Grass")

alpha = 0.01
res.list <- list()
plot.list <- list()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.16s.treat, contrast = contrast.list[[i]], pAdjustMethod = "bonferroni")
  #filter results by p-value
  res.alpha <- res[which(res$padj < alpha), ]
  #Bind taxonomy to results
  res.alpha.tax = cbind(as(res.alpha, "data.frame"), as(tax_table(ps.16s.nocontrols.nobss)[rownames(res.alpha), ], "matrix"))
  #tidy results 
  res.list[[paste(i,sep = ".")]] <- tidy(res.alpha)
  
  #generate plot of significant ASVs for each contrast
  # Order order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Order, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Order = factor(as.character(res.alpha.tax$Order), levels=names(x))
  # Genus order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Genus, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Genus = factor(as.character(res.alpha.tax$Genus), levels=names(x))
  p <- ggplot(res.alpha.tax, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(plot.name.list[[i]])
  plot.list[[i]] = p
}

#plot results

plot.list$SA.G + plot.list$ST.G + plot.list$SAG.G + plot.list$STG.G

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

write.csv(df.res, "Data/16s.ddseq.treatment.csv")

 
## 16S: G v F ##
#make "grass" the result to compare to

sample_data(ps.16s.nocontrols.nobss)$FunGroup <- relevel(sample_data(ps.16s.nocontrols.nobss)$FunGroup, "Grass")

treat.16s = phyloseq_to_deseq2(ps.16s.nocontrols.nobss, ~ FunGroup)

dds.16s.treat = DESeq(treat.16s, test = "Wald", fitType = "parametric")

resultsNames(dds.16s.treat) #gives us comparisons for our contrasts

contrasts = c("F.G", "FG.G")

contrast.list <- list(F.G = c("FunGroup", "Forb", "Grass"),
                      FG.G = c("FunGroup", "grass_x_forb", "Grass"))


plot.name.list <- list(F.G = "Forb vs. Grass",
                       FG.G = "Forb x grass vs. Grass")

alpha = 0.01
res.list <- list()
plot.list <- list()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.16s.treat, contrast = contrast.list[[i]], pAdjustMethod = "bonferroni")
  #filter results by p-value
  res.alpha <- res[which(res$padj < alpha), ]
  #Bind taxonomy to results
  res.alpha.tax = cbind(as(res.alpha, "data.frame"), as(tax_table(ps.16s.nocontrols.nobss)[rownames(res.alpha), ], "matrix"))
  #tidy results 
  res.list[[paste(i,sep = ".")]] <- tidy(res.alpha)
  
  #generate plot of significant ASVs for each contrast
  # Order order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Order, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Order = factor(as.character(res.alpha.tax$Order), levels=names(x))
  # Genus order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Genus, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Genus = factor(as.character(res.alpha.tax$Genus), levels=names(x))
  p <- ggplot(res.alpha.tax, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(plot.name.list[[i]])
  plot.list[[i]] = p
}

#plot results

plot.list$F.G+ plot.list$FG.G 

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

write.csv(df.res, "Data/16s.ddseq.fungroup.csv")


#same as above but now with "forb" as the comparison
sample_data(ps.16s.nocontrols.nobss)$FunGroup <- relevel(sample_data(ps.16s.nocontrols.nobss)$FunGroup, "Forb")

treat.16s = phyloseq_to_deseq2(ps.16s.nocontrols.nobss, ~ FunGroup)

dds.16s.treat = DESeq(treat.16s, test = "Wald", fitType = "parametric")

resultsNames(dds.16s.treat) #gives us comparisons for our contrasts

contrasts = c("G.F", "FG.F")

contrast.list <- list(G.F = c("FunGroup", "Grass", "Forb"),
                      FG.F = c("FunGroup", "grass_x_forb", "Forb"))


plot.name.list <- list(G.F = "Grass vs. Forb",
                       FG.F = "Forb x grass vs. Forb")

alpha = 0.01
res.list <- list()
plot.list <- list()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.16s.treat, contrast = contrast.list[[i]], pAdjustMethod = "bonferroni")
  #filter results by p-value
  res.alpha <- res[which(res$padj < alpha), ]
  #Bind taxonomy to results
  res.alpha.tax = cbind(as(res.alpha, "data.frame"), as(tax_table(ps.16s.nocontrols.nobss)[rownames(res.alpha), ], "matrix"))
  #tidy results 
  res.list[[paste(i,sep = ".")]] <- tidy(res.alpha)
  
  #generate plot of significant ASVs for each contrast
  # Order order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Order, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Order = factor(as.character(res.alpha.tax$Order), levels=names(x))
  # Genus order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Genus, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Genus = factor(as.character(res.alpha.tax$Genus), levels=names(x))
  p <- ggplot(res.alpha.tax, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(plot.name.list[[i]])
  plot.list[[i]] = p
}

#plot results

plot.list$G.F+ plot.list$FG.F 

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

write.csv(df.res, "Data/16s.ddseq.fungroup.v2.csv")

# More ASVs are high abundance in grass x forb vs. forbs alone (vs. grass x forb vs. grasses alone)
# Does this suggest that forbs shift more than grasses during competion?
# Or that grasses are stronger drivers of rhizosphere communities?


## 16S: Competition ##
ps.16s.nocontrols.comp <- subset_samples(ps.16s.nocontrols, !(is.na(Competion)))

treat.16s = phyloseq_to_deseq2(ps.16s.nocontrols.nobss, ~ Competion)

dds.16s.treat = DESeq(treat.16s, test = "Wald", fitType = "parametric")

resultsNames(dds.16s.treat) #gives us comparisons for our contrasts


alpha = 0.01

#get results 
res <- results(dds.16s.treat, pAdjustMethod = "bonferroni")

#filter results by p-value
res.alpha <- res[which(res$padj < alpha), ]

# SURPRISE - no differentially abundant ASVs between one vs. two species
# But we see differences when split by forb / grass / forb x grass
# maybe this is too high level and the differences are subtle? or single species driven?

# So does this make sense to also do by plant species specifically? 
# Can certainly do this with the contrasts, not sure if necesary?











#### Differential abundance analysis: ITS ####

# Write scripts to test for ASVs that differ with:
# Plant species
# Competition
# These should mirror ordination comparisons, etc

# Can be done at ASV level or at higher levels of taxonomic classification (e.g. family, order)
# Starting with only the ASV level

## ITS: Soil v Rhiz  ##

treat.its = phyloseq_to_deseq2(ps.ITS.nocontrols, ~ SampleSubType)

dds.its.treat = DESeq(treat.its, test = "Wald", fitType = "parametric")

resultsNames(dds.its.treat) #gives us comparisons for our contrasts, if needed

alpha = 0.01

#get results 
res <- results(dds.its.treat, pAdjustMethod = "bonferroni")

#filter results by p-value
res.alpha <- res[which(res$padj < alpha), ]

#Bind taxonomy to results
res.alpha.tax = cbind(as(res.alpha, "data.frame"), as(tax_table(ps.ITS.nocontrols)[rownames(res.alpha), ], "matrix"))

#tidy results 
res <- tidy(res.alpha)

#generate plot of significant ASVs for each contrast
# Order order
x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Order, function(x) max(x))
x = sort(x, TRUE)
res.alpha.tax$Order = factor(as.character(res.alpha.tax$Order), levels=names(x))
# Genus order
x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Genus, function(x) max(x))
x = sort(x, TRUE)
res.alpha.tax$Genus = factor(as.character(res.alpha.tax$Genus), levels=names(x))
p <- ggplot(res.alpha.tax, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("Plant Rhizosphere vs. Background soil")

#plot results

p

#tidy results into a table to save
df.res <- as.data.frame(res)
names(df.res)[1] <- "ASV"

write.csv(df.res, "Data/its.ddseq.soilrhiz.csv")


## 16S:  G v SA v ST  ##
ps.its.nocontrols.nobss <- subset_samples(ps.ITS.nocontrols, SampleSubType !="Background_soil_sand_mix")

treat.its = phyloseq_to_deseq2(ps.its.nocontrols.nobss, ~ TreatmentName)

dds.its.treat = DESeq(treat.its, test = "Wald", fitType = "parametric")

resultsNames(dds.its.treat) #gives us comparisons for our contrasts

contrasts = c("SA.G", "ST.G",
              "SAG.G", "STG.G")

contrast.list <- list(SA.G = c("TreatmentName", "Stress_avoiding_forb", "Invasive_grass"),
                      ST.G = c("TreatmentName", "Stress_tolerant_forb", "Invasive_grass"),
                      SAG.G = c("TreatmentName", "SA_forb_X_grass", "Invasive_grass"),
                      STG.G = c("TreatmentName", "ST_forb_X_grass", "Invasive_grass"))


plot.name.list <- list(SA.G = "SA forb vs. Grass",
                       ST.G = "ST forb vs. Grass",
                       SAG.G = "SA forb X grass vs. Grass",
                       STG.G = "ST forb X grass vs. Grass")

alpha = 0.01
res.list <- list()
plot.list <- list()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.its.treat, contrast = contrast.list[[i]], pAdjustMethod = "bonferroni")
  #filter results by p-value
  res.alpha <- res[which(res$padj < alpha), ]
  #Bind taxonomy to results
  res.alpha.tax = cbind(as(res.alpha, "data.frame"), as(tax_table(ps.its.nocontrols.nobss)[rownames(res.alpha), ], "matrix"))
  #tidy results 
  res.list[[paste(i,sep = ".")]] <- tidy(res.alpha)
  
  #generate plot of significant ASVs for each contrast
  # Order order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Order, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Order = factor(as.character(res.alpha.tax$Order), levels=names(x))
  # Genus order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Genus, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Genus = factor(as.character(res.alpha.tax$Genus), levels=names(x))
  p <- ggplot(res.alpha.tax, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(plot.name.list[[i]])
  plot.list[[i]] = p
}

#plot results

plot.list$SA.G + plot.list$ST.G + plot.list$SAG.G + plot.list$STG.G

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

write.csv(df.res, "Data/its.ddseq.treatment.csv")


## 16S: G v F ##
#make "grass" the result to compare to

sample_data(ps.its.nocontrols.nobss)$FunGroup <- relevel(sample_data(ps.its.nocontrols.nobss)$FunGroup, "Grass")

treat.its = phyloseq_to_deseq2(ps.its.nocontrols.nobss, ~ FunGroup)

dds.its.treat = DESeq(treat.its, test = "Wald", fitType = "parametric")

resultsNames(dds.its.treat) #gives us comparisons for our contrasts

contrasts = c("F.G", "FG.G")

contrast.list <- list(F.G = c("FunGroup", "Forb", "Grass"),
                      FG.G = c("FunGroup", "grass_x_forb", "Grass"))


plot.name.list <- list(F.G = "Forb vs. Grass",
                       FG.G = "Forb x grass vs. Grass")

alpha = 0.01
res.list <- list()
plot.list <- list()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.its.treat, contrast = contrast.list[[i]], pAdjustMethod = "bonferroni")
  #filter results by p-value
  res.alpha <- res[which(res$padj < alpha), ]
  #Bind taxonomy to results
  res.alpha.tax = cbind(as(res.alpha, "data.frame"), as(tax_table(ps.its.nocontrols.nobss)[rownames(res.alpha), ], "matrix"))
  #tidy results 
  res.list[[paste(i,sep = ".")]] <- tidy(res.alpha)
  
  #generate plot of significant ASVs for each contrast
  # Order order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Order, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Order = factor(as.character(res.alpha.tax$Order), levels=names(x))
  # Genus order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Genus, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Genus = factor(as.character(res.alpha.tax$Genus), levels=names(x))
  p <- ggplot(res.alpha.tax, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(plot.name.list[[i]])
  plot.list[[i]] = p
}

#plot results

plot.list$F.G+ plot.list$FG.G 

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

write.csv(df.res, "Data/its.ddseq.fungroup.csv")


#same as above but now with "forb" as the comparison
sample_data(ps.its.nocontrols.nobss)$FunGroup <- relevel(sample_data(ps.its.nocontrols.nobss)$FunGroup, "Forb")

treat.its = phyloseq_to_deseq2(ps.its.nocontrols.nobss, ~ FunGroup)

dds.its.treat = DESeq(treat.its, test = "Wald", fitType = "parametric")

resultsNames(dds.its.treat) #gives us comparisons for our contrasts

contrasts = c("G.F", "FG.F")

contrast.list <- list(G.F = c("FunGroup", "Grass", "Forb"),
                      FG.F = c("FunGroup", "grass_x_forb", "Forb"))


plot.name.list <- list(G.F = "Grass vs. Forb",
                       FG.F = "Forb x grass vs. Forb")

alpha = 0.01
res.list <- list()
plot.list <- list()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.its.treat, contrast = contrast.list[[i]], pAdjustMethod = "bonferroni")
  #filter results by p-value
  res.alpha <- res[which(res$padj < alpha), ]
  #Bind taxonomy to results
  res.alpha.tax = cbind(as(res.alpha, "data.frame"), as(tax_table(ps.its.nocontrols.nobss)[rownames(res.alpha), ], "matrix"))
  #tidy results 
  res.list[[paste(i,sep = ".")]] <- tidy(res.alpha)
  
  #generate plot of significant ASVs for each contrast
  # Order order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Order, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Order = factor(as.character(res.alpha.tax$Order), levels=names(x))
  # Genus order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Genus, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Genus = factor(as.character(res.alpha.tax$Genus), levels=names(x))
  p <- ggplot(res.alpha.tax, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(plot.name.list[[i]])
  plot.list[[i]] = p
}

#plot results

plot.list$G.F+ plot.list$FG.F 

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

write.csv(df.res, "Data/its.ddseq.fungroup.v2.csv")

# More ASVs are high abundance in grass x forb vs. forbs alone (vs. grass x forb vs. grasses alone)
# Does this suggest that forbs shift more than grasses during competion?
# Or that grasses are stronger drivers of rhizosphere communities?
# See same results in 16S communities


## 16S: Competition ##
ps.its.nocontrols.comp <- subset_samples(ps.ITS.nocontrols, !(is.na(Competion)))

treat.its = phyloseq_to_deseq2(ps.its.nocontrols.nobss, ~ Competion)

dds.its.treat = DESeq(treat.its, test = "Wald", fitType = "parametric")

resultsNames(dds.its.treat) #gives us comparisons for our contrasts


alpha = 0.01

#get results 
res <- results(dds.its.treat, pAdjustMethod = "bonferroni")

#filter results by p-value
res.alpha <- res[which(res$padj < alpha), ]

# SURPRISE - no differentially abundant ASVs between one vs. two species 
# Same as in the 16S!
# But we see differences when split by forb / grass / forb x grass
# maybe this is too high level and the differences are subtle? or single species driven?

# So does this make sense to also do by plant species specifically? 
# Can certainly do this with the contrasts, not sure if necesary?




#### SourceTracker #####
#Attempting to use sourcetracker to assess % of community originating from forb / grass during competion

#### SourceTracker: 16S ####

#Load 16s metadata
metadata_16s <- read.csv('Data/ps.nocontrols.16s.rare9434.mappingdata.csv', row.names = 1)

#Load 16s ASVs
asvs_16s <- read.csv('Data/ps.nocontrols.16s.rare9434.csv', row.names = 1)

# Manipulate metadata to add Env and SourceSink columns
#Add "Env" column, this could be rhiz-vs.-background soil or specific treatments (e.g. pairwise plants)

#Env column for treatment
#Here we will use "background soil" as the source
metadata_16s  %<>% 
  mutate(Env = fct_explicit_na(TreatmentName, "Background_soil")) 

#Here we are using 'FunGroup' for Env  
#Here we will use "background soil", "forb" and "grass" as sources
metadata_16s  %<>% 
  mutate(Env.v2 = ifelse(Env == "Background_soil" & FunGroup == "grass_x_forb", "Background_soil", as.character(FunGroup)))

#Add source-sink columns for each Env 
metadata_16s %<>%
  mutate (SourceSink = ifelse(Env == "Background_soil", "source", "sink"),
          SourceSink.v2 = ifelse(Env.v2 =="Background_soil" | Env.v2 == "Forb" | Env.v2 == "Grass", "source", "sink"))

#fix rownames
row.names(metadata_16s) <- metadata_16s$SampleID_Fix

# extract only those samples in common between the two tables 
# note, since we ported from phyloseq both the mapping data and the asvs
# this should always have teh same number of samples! 
# but let's just check anyway
common.sample.ids_16s <- intersect(rownames(metadata_16s), rownames(asvs_16s))
asvs_16s <- asvs_16s[common.sample.ids_16s,]
metadata_16s <- metadata_16s[common.sample.ids_16s,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids_16s) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}



## Run SourceTracker on 'Treatment' ##

# extract the source environments and source/sink indices
train.ix <- which(metadata_16s$SourceSink=='source')
test.ix <- which(metadata_16s$SourceSink=='sink')
envs <- metadata_16s$Env
if(is.element('Description',colnames(metadata_16s))) desc <- metadata_16s$Description


# load SourceTracker package
source('Scripts/SourceTracker.r')


## Couldn't get the tuning to work so skipped ##
# tune the alpha values using cross-validation (this is slow!)
#tune.results <- tune.st(asvs_16s[train.ix,], envs[train.ix])
#alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st <- sourcetracker(asvs_16s[train.ix,], envs[train.ix])

# the next command was too slow to finish running, even overnight on Cassie's old computer!
# However, it looked like it was working
# Results seemed to be 60% ASVs sourced from background soil, 40% from 'unknown' source
# Maybe Marina can try running??? Or Cassie can put the data on the cluster!

# Estimate source proportions in test data
results <- predict(st,asvs_16s[test.ix,], alpha1=alpha1, alpha2=alpha2)

# Estimate leave-one-out source proportions in training data 
results.train <- predict(st, alpha1=alpha1, alpha2=alpha2)

# plot results
labels <- sprintf('%s %s', envs,desc)
plot(results, labels[test.ix], type='pie')

# other plotting functions
# plot(results, labels[test.ix], type='bar')
# plot(results, labels[test.ix], type='dist')
# plot(results.train, labels[train.ix], type='pie')
# plot(results.train, labels[train.ix], type='bar')
# plot(results.train, labels[train.ix], type='dist')

# plot results with legend
# plot(results, labels[test.ix], type='pie', include.legend=TRUE, env.colors=c('#47697E','#5B7444','#CC6666','#79BEDB','#885588'))


## Run SourceTracker on 'FunGroup' ##

# extract the source environments and source/sink indices
train.ix.fg <- which(metadata_16s$SourceSink.v2=='source')
test.ix.fg <- which(metadata_16s$SourceSink.v2=='sink')
envs.fg <- metadata_16s$Env.v2
if(is.element('Description',colnames(metadata_16s))) desc <- metadata_16s$Description

## Couldn't get the tuning to work above so skipped ##
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st.fg <- sourcetracker(asvs_16s[train.ix.fg,], envs.fg[train.ix.fg])


# the next command was too slow to finish running, even overnight on Cassie's old computer!
# However, it looked like it was working
# Results seemed to be: 13% from background, 48% from Forb, 30% from grass, 9% unknown
# INTERESTING!

# Estimate source proportions in test data
results.fg <- predict(st.fg,asvs_16s[test.ix.fg,], alpha1=alpha1, alpha2=alpha2)

# Estimate leave-one-out source proportions in training data 
results.train.fg <- predict(st.fg, alpha1=alpha1, alpha2=alpha2)

# plot results
labels.fg <- sprintf('%s %s', envs.fg,desc)
plot(results.fg, labels.fg[test.ix.fg], type='pie')




#### Sourcetracker: ITS ####

#Load ITS metadata
metadata_its <- read.csv('Data/ps.nocontrols.its.rare7557.mappingdata.csv', row.names = 1)

#Load ITS ASVs
asvs_its <- read.csv('Data/ps.nocontrols.its.rare7557.csv', row.names = 1)

# Manipulate metadata to add Env and SourceSink columns
#Add "Env" column, this could be rhiz-vs.-background soil or specific treatments (e.g. pairwise plants)

#Env column for treatment
#Here we will use "background soil" as the source
metadata_its  %<>% 
  mutate(Env = fct_explicit_na(TreatmentName, "Background_soil")) 

#Here we are using 'FunGroup' for Env  
#Here we will use "background soil", "forb" and "grass" as sources
metadata_its  %<>% 
  mutate(Env.v2 = ifelse(Env == "Background_soil" & FunGroup == "grass_x_forb", "Background_soil", as.character(FunGroup)))

#Add source-sink columns for each Env 
metadata_its %<>%
  mutate (SourceSink = ifelse(Env == "Background_soil", "source", "sink"),
          SourceSink.v2 = ifelse(Env.v2 =="Background_soil" | Env.v2 == "Forb" | Env.v2 == "Grass", "source", "sink"))

#fix rownames
row.names(metadata_its) <- metadata_its$SampleID_Fix


## Run SourceTracker on 'Treatment' ##

# extract the source environments and source/sink indices
train.ix.its <- which(metadata_its$SourceSink=='source')
test.ix.its <- which(metadata_its$SourceSink=='sink')
envs.its <- metadata_its$Env
if(is.element('Description',colnames(metadata_its))) desc <- metadata_its$Description


# load SourceTracker package
source('Scripts/SourceTracker.r')


## Couldn't get the tuning to work so skipped ##
# tune the alpha values using cross-validation (this is slow!)
#tune.results <- tune.st(asvs_16s[train.ix,], envs[train.ix])
#alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st.its <- sourcetracker(asvs_its[train.ix.its,], envs.its[train.ix.its])

# the next command was too slow to finish running, even overnight on Cassie's old computer!
# However, it looked like it was working
# Results seemed to be: 20% background, 80% unknown - Much less than the 16S!

# Estimate source proportions in test data
results.its <- predict(st.its,asvs_its[test.ix.its,], alpha1=alpha1, alpha2=alpha2)

# Estimate leave-one-out source proportions in training data 
results.train.its <- predict(st.its, alpha1=alpha1, alpha2=alpha2)

# plot results
labels.its <- sprintf('%s %s', envs.its,desc)
plot(results.its, labels.its[test.ix.its], type='pie')

# other plotting functions
# plot(results, labels[test.ix], type='bar')
# plot(results, labels[test.ix], type='dist')
# plot(results.train, labels[train.ix], type='pie')
# plot(results.train, labels[train.ix], type='bar')
# plot(results.train, labels[train.ix], type='dist')

# plot results with legend
# plot(results, labels[test.ix], type='pie', include.legend=TRUE, env.colors=c('#47697E','#5B7444','#CC6666','#79BEDB','#885588'))


## Run SourceTracker on 'FunGroup' ##

# extract the source environments and source/sink indices
train.ix.fg.its <- which(metadata_its$SourceSink.v2=='source')
test.ix.fg.its <- which(metadata_its$SourceSink.v2=='sink')
envs.fg.its <- metadata_its$Env.v2
if(is.element('Description',colnames(metadata_its))) desc <- metadata_its$Description

## Couldn't get the tuning to work above so skipped ##
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st.fg.its <- sourcetracker(asvs_its[train.ix.fg.its,], envs.fg.its[train.ix.fg.its])

# the next command was too slow to finish running, even overnight on Cassie's old computer!
# However, it looked like it was working
# Results seemed to be: 4% background soil, 45% forb, 25% grass, 25% unknown

# Estimate source proportions in test data
results.fg.its <- predict(st.fg.its,asvs_its[test.ix.fg.its,], alpha1=alpha1, alpha2=alpha2)

# Estimate leave-one-out source proportions in training data 
results.train.fg.its <- predict(st.fg.its, alpha1=alpha1, alpha2=alpha2)

# plot results
labels.fg.its <- sprintf('%s %s', envs.fg.its,desc)
plot(results.fg.its, labels.fg.its[test.ix.fg.its], type='pie')