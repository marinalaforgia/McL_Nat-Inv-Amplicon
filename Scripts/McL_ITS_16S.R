### 16S and ITS Analysis April 2020 ###
rm(list=ls())
#### Setup ####

#setwd("/Users/Marina/Documents/UC-Davis/Projects/McL_Nat-Inv-Amplicon") # Marina's Working Directory
#setwd("") # Cassie's Working Directory

### Load Libraries ###

library(EcolUtils)
library(ggplot2)
library(plyr)
library(multcomp)
library(FSA)
library(Rmisc)
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

### Read in Data ###

## Remove background soil for analysis ##

## 16S ##
ps.16s.nocontrols.rare <- readRDS("Data/16S/Intermediates/ps-nocontrols-rare9434.RDS") # rarefied 16s data

## ITS ##
ps.ITS.nocontrols.rare <- readRDS("Data/ITS/greenhouse_its_itsx_decontam_controlsremoved_rare7557.rds") # rarefied ITS data

# Metadata mapping file including biomass and traits (version4)
mapping <- read.csv("/Users/Marina/Documents/UC-Davis/Projects/McL_Nat-Inv-Amplicon/Data/Setup/Grassland-Amplicon-Mapping-File4.csv") 
mapping$FunGroup <- ifelse(mapping$TreatmentName == "Invasive_grass", "Grass", NA)
mapping$FunGroup <- ifelse(mapping$Competion == "SingleSpecies" & mapping$TreatmentName != "Invasive_grass", "Forb", mapping$FunGroup)
mapping[is.na(mapping$FunGroup),]$FunGroup <- "grass_x_forb"
mapping$FunGroup <- as.factor(mapping$FunGroup)

# 16S mapping file does not match up with the mapping file for the ITS data, adjust here 
for(i in sample_names(ps.16s.nocontrols.rare)) {
    sample_names(ps.16s.nocontrols.rare)[which(sample_names(ps.16s.nocontrols.rare) == i)] <- as.character(mapping[mapping$SampleID.16S == i,]$SampleID_Fix)
  }

sample_names(ps.16s.nocontrols.rare)
sample_names(ps.ITS.nocontrols.rare) # they match up!

df.16s <- data.frame(sample_data(ps.16s.nocontrols.rare))
df.ITS <- data.frame(sample_data(ps.ITS.nocontrols.rare))

missing <- df.ITS[!(df.ITS$SampleID_Fix %in% df.16s$SampleID_Fix),]
missing2 <- df.16s[!(df.16s$SampleID_Fix %in% df.ITS$SampleID_Fix),]
row.names(mapping) <- mapping$SampleID_Fix
sample_data(ps.ITS.nocontrols.rare) <- mapping[which(mapping$SampleID_Fix %in% ps.ITS.nocontrols.rare@sam_data$SampleID_Fix),]
sample_data(ps.ITS.nocontrols.rare) # fixed!

sample_data(ps.16s.nocontrols.rare) <- mapping

rm(missing, missing2, df.16s, df.ITS, i)

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

mantel(dist_ITS, dist_16s, method = "spearman", permutations = 9999) # non parametric, they are highly correlated

