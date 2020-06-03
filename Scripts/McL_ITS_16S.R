### 16S and ITS Analysis April 2020 ###
rm(list=ls())

#### Setup ####

#setwd("/Users/Marina/Documents/UC-Davis/Projects/McL_Nat-Inv-Amplicon") # Marina's Working Directory
#setwd("/Users/Cassie/Documents/Dropbox/Greenhouse_ITS/McL_Nat-Inv-Amplicon") # Cassie's Working Directory

#### Load Libraries ####

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
library(reshape2)
library(betapart)
library(dada2)
library(extrafont)
library(tidyverse) # for tidy data packages, automatically loads dplyr
library(magrittr) # for piping
library(DESeq2)
library(biobroom)
library(patchwork)
library(VennDiagram)

#### Read in Data ####

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
mapping <- read.csv("Data/Grassland-Amplicon-Mapping-File4.csv") #Cassie's path
# 
# #### Mapping file matching NOT NEEDED ####
# #Adding columns to the mapping file
# mapping$FunGroup <- ifelse(mapping$TreatmentName == "Invasive_grass", "Grass", NA)
# mapping$FunGroup <- ifelse(mapping$Competion == "SingleSpecies" & mapping$TreatmentName != "Invasive_grass", "Forb", mapping$FunGroup)
# 
# mapping$FunGroup <- ifelse(mapping$Competion == "TwoSpecies", "grass_x_forb", as.character(mapping$FunGroup))
# 
# write.csv(mapping, "/Users/Marina/Documents/UC-Davis/Projects/McL_Nat-Inv-Amplicon/Data/Setup/Grassland-Amplicon-Mapping-File4.csv", row.names = F)

# #16S mapping file does not match up with the mapping file for the ITS data, adjust here
# for(i in sample_names(ps.16s.nocontrols)) {
#   sample_names(ps.16s.nocontrols)[which(sample_names(ps.16s.nocontrols) == i)] <- as.character(mapping[mapping$SampleID.16S == i,]$SampleID_Fix)
# }
# 
# sample_names(ps.16s.nocontrols)
# sample_names(ps.ITS.nocontrols) # they match up!
# 
# df.16s <- data.frame(sample_data(ps.16s.nocontrols))
# df.ITS <- data.frame(sample_data(ps.ITS.nocontrols))
# 
# missing <- df.ITS[!(df.ITS$SampleID_Fix %in% df.16s$SampleID_Fix),]
# missing2 <- df.16s[!(df.16s$SampleID_Fix %in% df.ITS$SampleID_Fix),]
# row.names(mapping) <- mapping$SampleID_Fix
# sample_data(ps.ITS.nocontrols) <- mapping[which(mapping$SampleID_Fix %in% ps.ITS.nocontrols@sam_data$SampleID_Fix),]
# sample_data(ps.ITS.nocontrols) # fixed!
# 
# sample_data(ps.16s.nocontrols) <- mapping
# 
# rm(missing, missing2, df.16s, df.ITS, i)
# 
# saveRDS(ps.16s.nocontrols, "Data/16S/Intermediates/ps-nocontrols-unrare.RDS")
# saveRDS(ps.ITS.nocontrols, "Data/ITS/greenhouse_its_itsx_decontam_controlsremoved.rds")

#### DESeq fix taxonomy ####

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

#### Rarefied Data ####

## 16S ##
# ps.nocontrols.16s.rare9434 <- rarefy_even_depth(ps.16s.nocontrols, 9434, replace = FALSE, rngseed = 5311)
# saveRDS(ps.nocontrols.16s.rare9434, "Data/16S/Intermediates/ps-nocontrols-rare9434.RDS")
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
# ps.nocontrols.its.rare7557 <- rarefy_even_depth(ps.ITS.nocontrols, 7557, replace = FALSE, rngseed = 5311)
# saveRDS(ps.nocontrols.its.rare7557, "Data/ITS/greenhouse_its_itsx_decontam_controlsremoved_rare7557.rds")
# # 132 samples

# Save for first pass at Halla (https://huttenhower.sph.harvard.edu/halla)
#write.table(t(otu_table(ps.nocontrols.its.rare7557)), "Data/ps.nocontrols.its.rare7557.transposed.txt", sep = "\t")



#Save for first pass sourcetracker (https://github.com/danknights/sourcetracker/)
#write.csv(otu_table(ps.nocontrols.its.rare7557), "Data/ps.nocontrols.its.rare7557.csv")
#write.csv(sample_data(ps.nocontrols.its.rare7557), "Data/ps.nocontrols.its.rare7557.mappingdata.csv")


# Load 
ps.ITS.nocontrols.rare <- readRDS("Data/ITS/greenhouse_its_itsx_decontam_controlsremoved_rare7557.rds") # rarefied ITS data
ps.ITS.nocontrols.rare <- readRDS("Data/greenhouse_its_itsx_decontam_controlsremoved_rare7557.rds") # rarefied ITS data, Cassie's path


## 16s ##

ps.16s.nocontrols.rare.spp <- subset_samples(ps.16s.nocontrols.rare, Competion == "SingleSpecies")

ps.16s.nocontrols.rare.spp.RA <- transform_sample_counts(ps.16s.nocontrols.rare.spp, function(x) x / sum(x)) # relative abundance (used later)

## ITS ##

ps.ITS.nocontrols.rare.spp <- subset_samples(ps.ITS.nocontrols.rare, Competion == "SingleSpecies") 

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

#### Ordination (16S): spp ####

###
# Forb Species
###

ps.16s.nocontrols.rare.spp.f <- subset_samples(ps.16s.nocontrols.rare.spp, grass == "none")

ps.16s.nc.rare.spp.ord <- ordinate(ps.16s.nocontrols.rare.spp.f, "PCoA", "wunifrac")

plot_ordination(ps.16s.nocontrols.rare.spp.f, ps.16s.nc.rare.spp.ord, shape = "Plant SpeciesSampled", color = "PlantSpeciesSampled") +
  theme_bw(base_size = 20) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = PlantSpeciesSampled)) 

DistWU.Spp <- phyloseq::distance(ps.16s.nocontrols.rare.spp.f, method = "wunifrac", type = "samples")

adonis(DistWU.Spp ~ PlantSpeciesSampled, as(sample_data(ps.16s.nocontrols.rare.spp.f), "data.frame"), permutations = 9999) # not significantly different among species... though really effing close

C <- betadisper(DistWU.Spp, as(sample_data(ps.16s.nocontrols.rare.spp.f), "data.frame")$PlantSpeciesSampled)
permutest(C, permutations = 9999) # spread doesnt vary between forbs species

### 
# Grass species
###
ps.16s.nocontrols.rare.spp.g <- subset_samples(ps.16s.nocontrols.rare.spp, forb == "none")

ps.16s.nc.rare.spp.ord <- ordinate(ps.16s.nocontrols.rare.spp.g, "PCoA", "wunifrac")

plot_ordination(ps.16s.nocontrols.rare.spp.g, ps.16s.nc.rare.spp.ord, shape = "Plant SpeciesSampled", color = "PlantSpeciesSampled") +
  theme_bw(base_size = 20) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = PlantSpeciesSampled)) 


DistWU.Spp <- phyloseq::distance(ps.16s.nocontrols.rare.spp.g, method = "wunifrac", type = "samples")

adonis(DistWU.Spp ~ PlantSpeciesSampled, as(sample_data(ps.16s.nocontrols.rare.spp.g), "data.frame"), permutations = 9999) # grass microbiomes are sig different, particularly medusahead

C <- betadisper(DistWU.Spp, as(sample_data(ps.16s.nocontrols.rare.spp.g), "data.frame")$PlantSpeciesSampled)
permutest(C, permutations = 9999)

pair.spp.g <- adonis.pair(DistWU.Spp, as(sample_data(ps.16s.nocontrols.rare.spp.g), "data.frame")$grass, nper = 9999, corr.method = "BH")
kable(pair.spp.g) # TACA sig diff but also probably dont have enough samples to really do this



#### Ordination (16S): Competition ####
ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion)))
#ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, forb != "none")

sample_data(ps.16s.nocontrols.rare)$FunGroup <- recode(sample_data(ps.16s.nocontrols.rare)$FunGroup, grass_x_forb = "Competition")

sample_data(ps.16s.nocontrols.rare)$FunGroup <- factor(sample_data(ps.16s.nocontrols.rare)$FunGroup, levels = c("Forb", "Grass", "Competition"))
  
ps.rare.ord.tr <- ordinate(ps.16s.nocontrols.rare, "PCoA", "wunifrac")

ord.plot <- plot_ordination(ps.16s.nocontrols.rare, ps.rare.ord.tr, color = "FunGroup", shape = "FunGroup") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = FunGroup)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(.82, .82),
    legend.text = element_text(size = 10),
    legend.margin = margin(c(.1,.1,.1,.1))
  ) +
  scale_color_manual(values = c("magenta4", "cyan4" ,"goldenrod3"))

ggsave("Figures/ord-plot.jpg", ord.plot, width = 5, height = 4, units = "in", dpi = 600)

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

pair.am2 <- adonis.pair(DistWU, as(sample_data(ps.16s.nocontrols.rare), "data.frame")$TreatmentName, corr.method = "none")
pair.am2 <- pair.am2[-c(8,7),1:6]
pair.am2 <- cbind(pair.am2, p.val.corr = p.adjust(pair.am2$P.value, method = "BH", n = length(pair.am2$P.value)))
kable(pair.am2)

#sample_data(ps.16s.nocontrols.rare)$FunGroup <- as.factor(sample_data(ps.16s.nocontrols.rare)$FunGroup)
pair.am3 <- adonis.pair(DistWU, as(sample_data(ps.16s.nocontrols.rare), "data.frame")$FunGroup, nper = 9999, corr.method = "BH")
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


sample_data(ps.ITS.nocontrols.rare)$FunGroup <- recode(sample_data(ps.ITS.nocontrols.rare)$FunGroup, grass_x_forb = "Competition")

sample_data(ps.ITS.nocontrols.rare)$FunGroup <- factor(sample_data(ps.ITS.nocontrols.rare)$FunGroup, levels = c("Forb", "Grass", "Competition"))

ord.plot.ITS <- plot_ordination(ps.ITS.nocontrols.rare, GL_pcoa, color = "FunGroup", shape = "FunGroup") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = FunGroup)) +
  theme(
    # legend.title = element_blank(),
    # legend.position = c(.82, .82),
    # legend.text = element_text(size = 10),
    # legend.margin = margin(c(.1,.1,.1,.1))
    legend.position = "none",
  ) +
  scale_color_manual(values = c("magenta4", "cyan4" ,"goldenrod3"))

ggsave("Figures/ord-plot-ITS.jpg", ord.plot.ITS, width = 5, height = 4, units = "in", dpi = 600)

DistBC = phyloseq::distance(ps.ITS.nocontrols.rare, method = "bray", type="samples")

adonis(DistBC ~ TreatmentName, as(sample_data(ps.ITS.nocontrols.rare), "data.frame"), permutations = 9999) # marginal treatment diffs

adonis(DistBC ~ Competion, as(sample_data(ps.ITS.nocontrols.rare), "data.frame"), permutations = 9999) # marginal competition shifts

adonis(DistBC ~ FunGroup, as(sample_data(ps.ITS.nocontrols.rare), "data.frame"), permutations = 9999) # marginal fungroup diffs

# because all these are marginal, we should not be looking at pairwise comparisons

sample_data(ps.ITS.nocontrols.rare)$FunGroup <- as.factor(sample_data(ps.ITS.nocontrols.rare)$FunGroup)
pair.its.com <- adonis.pair(DistBC, as(sample_data(ps.ITS.nocontrols.rare), "data.frame")$FunGroup, nper = 9999, corr.method = "BH")
kable(pair.its.com)

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
## Is this code correct? double check, ITS uses tax_glom; Also I changed this to look at family, that might make more sense however given that we probably arent using this anymore I'm not going through an changing them all
AvgRA99_Order <- filter_taxa(ps.16s.nocontrols.rare.RA, function(x) mean(x) > .007, TRUE)
df_Order <- psmelt(AvgRA99_Order)
grouped_Order <- group_by(df_Order, FunGroup, Family, Order, Class, Phylum, Kingdom)
avgs_Order <- dplyr::summarise(grouped_Order, mean=100*mean(Abundance), sd=100*sd(Abundance), se=100*se(Abundance))

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
taxa = unique(DataSet$Family)
df_taxa <- data.frame(taxa = taxa, chisq = rep(NA, length(taxa)), pvals = rep(NA, length(taxa)))

for (cat in df_taxa$taxa) { # for each Order that has a mean abundance above 1%
  new_df <- subset(DataSet, DataSet$Family == cat)
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
  new_df <- subset(DataSet, DataSet$Family == cat)
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
plot_richness(ps.16s.nocontrols.rare, measures = "Shannon", x = "FunGroup", color = "FunGroup") + 
  theme(text = element_text(size=24)) + 
  geom_boxplot() + 
  geom_jitter() +
  theme(axis.text.x = element_text(angle=-70, hjust=0, vjust=.5)) 

p.rich.16s <- plot_richness(ps.16s.nocontrols.rare, measures = "Shannon", x = "FunGroup", color = "FunGroup") + 
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    strip.text = element_blank()
  ) +
  scale_color_manual(values = c("magenta4", "cyan4" ,"goldenrod3")) +
  labs(y = "Shannon Diversity")
  
ggsave("Figures/rich-16S.jpg", p.rich.16s, width = 5, height = 4, units = "in", dpi = 600)

## Stats ##

GL_Alpha <- estimate_richness(ps.16s.nocontrols.rare, measures = c("Observed","Shannon", "InvSimpson", "Chao1"))
GL_Alpha2 <- cbind(GL_Alpha, sample_data(ps.16s.nocontrols.rare))

#kruskal_test(Observed ~ FunGroup, distribution = approximate(nresample = 9999), data = GL_Alpha2)

GL_Alpha2$FunGroup <- as.factor(GL_Alpha2$FunGroup)
kruskal_test(Shannon ~ FunGroup, distribution = approximate(nresample = 9999), data = GL_Alpha2)
# breaking news: shannon diversity does vary

dunnTest(Shannon ~ FunGroup, data = GL_Alpha2, method = "bh") # forbs have higher alpha diversity than grasses; competition treatments have marginally higher diversity than grasses, but not than forbs. forbs add diversity to competition

#### Alpha diversity: ITS ####
plot_richness(ps.ITS.nocontrols.rare, measures = c("Observed", "Shannon"), x = "FunGroup", color = "FunGroup") + 
  theme(text = element_text(size=24)) + geom_boxplot() + 
  geom_jitter() +
  theme(axis.text.x = element_text(angle=-70, hjust=0, vjust=.5)) 

p.rich.ITS <- plot_richness(ps.ITS.nocontrols.rare, measures = "Shannon", x = "FunGroup", color = "FunGroup") + 
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    strip.text = element_blank()
  ) +
  scale_color_manual(values = c("magenta4", "cyan4" ,"goldenrod3")) +
  labs(y = "Shannon Diversity")

ggsave("Figures/rich-ITS.jpg", p.rich.ITS, width = 5, height = 4, units = "in", dpi = 600)

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
### Note: just because an ASV/family/taxonomic unit does not appear in the CORE does not mean it does not appear in the SAMPLE, not sure how useful it is to break it up and look for changes in the core between these different groups #

ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion))) # these were likely removed earlier but just to be sure
otu <- as.data.frame(t(otu_table(ps.16s.nocontrols.rare)))

nReads <- 9434
map <- mapping

otu_PA <- 1*((otu>0)==1)  # presence-absence data (if OTU is present (greater than 0) assign a 1)
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  # total number of sites that OTU is present in, divided by the number of sites  (occupancy calculation)
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     # relative abundance: For each column divide every entry by the column total (this give relative abundance  of each OTU per sample), then give calculate the mean relative abundance per OTU by calculating the mean of each row
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance; occupancy of an OTU is average number of samples an OTU occurs in, abundance is the average relative abundance of that OTU within a sample

# for some reason this code replaces the - with a . in our sample names, so creating a new column to make it run (this is used later in PlotDF)
map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
  gather(SampleID.occ, abun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  group_by(otu, FunGroup) %>% # first try it by FunGroup?
  summarise(plot_freq = sum(abun > 0)/length(abun), # number of samples within a treatment where that OTU was present (greater than 0) divided by total number of samples within that treatment; so this is 1 if the OTU was present in every sample in that treatment and lower otherwise; so basically it's the percentage of subsamples within a treatment that OTU occurs in
            coreTrt = ifelse(plot_freq == 1, 1, 0), # Core Treatment = 1 if that OTU was present in every sample within a Treatment, 0 otherwise
            detect = ifelse(plot_freq > 0, 1, 0)) %>%    # 1 if that OTU was detected at all within that Treatment and 0 if not
  group_by(otu) %>% # group by OTU
  summarise(sumF = sum(plot_freq), # frequency of an OTU across treatments, so with 3 treatments (grass, forb, and grassxforb), a 3 would indicate that this OTU occurs in every treatment and every subsample of that treatment
            sumG = sum(coreTrt), # total number of Treatments where that OTU was present in every subsample
            nS = length(FunGroup)*2, # total number of Treatments (times 2) = 6
            Index = (sumF + sumG)/nS) # calculating weighting Index based on number of Treatments detected

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by = 'otu') %>% # why even need occ_abun here?
  transmute(otu = otu,
            rank = Index) %>%
  arrange(dplyr::desc(rank))
# 
# BCaddition <- NULL
# 
# otu_start <- otu_ranked$otu[1] # take the first ranked OTU
# start_matrix <- as.matrix(otu[otu_start,]) # extract that OTU's abundance per sample, this should be a one column vector
# #start_matrix <- t(start_matrix) # turns it into a one row vector, but mine already is because my original OTU table is a dataframe, theirs is already a matrix
# 
# x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) #Error in combn(ncol(start_matrix), 2) : n < m (fixed above); take every combination of samples, and for each combination give the absolute difference in those two OTU abundances, sum them all (sum what?) and divide by 2* number of reads (because 2 samples); i have no idea what the sum function is doing given that inside sum there's only one number, this is relevant later I believe?
# x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
# df_s <- data.frame(x_names,x)
# names(df_s)[2] <- 1 # i dont fully understand what this column is 
# BCaddition <- rbind(BCaddition,df_s)
#   
# for(i in 2:500){
#   otu_add = otu_ranked$otu[i] # for the top 500 ranked OTUs, starting with the second one beacuse we computed the first previously
#   add_matrix <- as.matrix(otu[otu_add,])
#   #add_matrix <- t(add_matrix) # again dont need this because its a dataframe
#   start_matrix <- rbind(start_matrix, add_matrix) # bind together the next ranked OTU with the previously ranked OTU
#   x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads)) # and compute the difference again
#   x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
#   df_a <- data.frame(x_names, x)
#   names(df_a)[2] <- i 
#   BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))  
# }
# 
# # and this code does it for the last one, I think 
# x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]] - otu[,x[2]]))/(2*nReads))
# x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse = ' - '))
# df_full <- data.frame(x_names, x)
# #names(df_full)[2] <- length(rownames(otu)) # this assumes we are looking at all OTUs, not just the top 500 ranked, right?
# names(df_full)[2] <- i + 1 #I think
# BCfull <- left_join(BCaddition, df_full, by = 'x_names') # each column is an OTU, each row a comparison of two samples, so the first (non-sample) column is the comparison of the first ranked OTU's abundance between the two samples specified, the second column is the difference in the second rank OTU's abundance between teh two samples specified and so on all the way until the 500th ranked OTU; it does this for each combination of samples
#  
# rownames(BCfull) <- BCfull$x_names
# temp_BC <- BCfull
# temp_BC$x_names <- NULL
# temp_BC_matrix <- as.matrix(temp_BC)
# 
# BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>% 
#   gather(comparison, BC, -rank) %>%
#   group_by(rank) %>%
#   summarise(MeanBC = mean(BC)) %>%
#   arrange(-dplyr::desc(MeanBC)) %>%
#   mutate(proportionBC = MeanBC/max(MeanBC))
# 
# Increase = BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
# increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase) + 1))))
# BC_ranked <- left_join(BC_ranked, increaseDF)
# BC_ranked <- BC_ranked[-nrow(BC_ranked),]
# 
# rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)
# 
# BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
# 
# elbow <- which.max(BC_ranked$fo_diffs)
# 
# (elbow.all.p <- ggplot(BC_ranked[1:250,], aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
#   geom_point(aes(y = proportionBC)) +
#   theme_classic() + 
#   theme(strip.background = element_blank(), 
#         axis.text.x = element_text(size = 7, angle = 45)) +
#   geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
#   geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])), 
#              lty = 3, col = 'blue', cex = .5) +
#   labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
#   annotate(geom = "text", 
#            x = elbow + 10, 
#            y = .15, 
#            label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +    
#   annotate(geom = "text", 
#            x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4, 
#            y = .08, 
#            label = paste("Last 2% increase (",last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))
# 
# occ_abun$fill <- 'no'
# occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))]] <- 'core'
# 
# saveRDS(occ_abun, "Data/16S/Core/occ_abun_asv_all.RDS")
# saveRDS(BC_ranked, "Data/16S/Core/BC_ranked_asv_all.RDS")

occ_abun.asv.all <- readRDS("Data/16S/Core/occ_abun_asv_all.RDS")
BC_ranked.asv.all <- readRDS("Data/16S/Core/BC_ranked_asv_all.RDS")

# Models for the whole community
taxon <- as.matrix(rownames(otu))
rownames(taxon) <- taxon
spp <- t(otu)
obs.np <- sncm.fit(spp = spp, taxon = taxon, stats = FALSE, pool = NULL)

ggplot() +
  geom_point(data = occ_abun.asv.all[occ_abun.asv.all$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun.asv.all[occ_abun.asv.all$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: All Samples", x = "log10(mean relative abundance)", y = "Occupancy")

core.all <- occ_abun.asv.all$otu[occ_abun.asv.all$fill == 'core']

otu_relabun.asv.all <- decostand(otu, method = "total", MARGIN = 2)

plotDF <-  data.frame(otu = as.factor(row.names(otu_relabun.asv.all)), otu_relabun.asv.all) %>%
  gather(SampleID.occ, relabun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  left_join(otu_ranked, by = 'otu') %>%
  filter(otu %in% core.all) %>% 
  group_by(FunGroup, otu) %>%
  summarise(plot_freq = sum(relabun>0)/length(relabun), 
            coreSite = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0))

plotDF$otu <- factor(plotDF$otu, levels = otu_ranked$otu[1:205])

(all.plotDF <- ggplot(plotDF, aes(x = otu, plot_freq, group = FunGroup, fill = FunGroup)) +    
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(plotDF$otu %in% core.all))) +
  theme(axis.text = element_text(size = 6)) +
  labs(x = 'Ranked OTUs', y = 'Occupancy by sample'))

## divide core into neutral vs deterministic vs dispersal limited ##
otu.mod.all <- data.frame(otu = occ_abun.asv.all[occ_abun.asv.all$fill == "core",]$otu, mod = rep(NA, 205))

for(i in otu.mod.all$otu){
otu.mod.all[otu.mod.all$otu == i,]$mod <- ifelse(log10(occ_abun.asv.all[occ_abun.asv.all$otu == i,]$otu_occ) > log10(obs.np[obs.np$V1 == i,]$pred.upr), "det", "neu")

otu.mod.all[otu.mod.all$otu == i,]$mod <- ifelse(log10(occ_abun.asv.all[occ_abun.asv.all$otu == i,]$otu_occ) < log10(obs.np[obs.np$V1 == i,]$pred.lwr), "dis", otu.mod.all[otu.mod.all$otu == i,]$mod)

}

#rm(occ_abun, obs.np, otu, otu_PA, otu_ranked, spp, taxon, otu_occ, otu_rel, i, nReads, PresenceSum)

# look at differences in positively selected for and negatively selected for ASVS
plotDF.s <- filter(plotDF, otu %in% otu.mod.all[otu.mod.all$mod != "neu",]$otu)
plotDF.dt <- filter(plotDF, otu %in% otu.mod.all[otu.mod.all$mod == "det",]$otu)

ggplot(plotDF.dt, aes(x = otu, plot_freq, group = FunGroup, fill = FunGroup)) +    
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(plotDF$otu %in% core.all))) +
  theme(axis.text = element_text(size = 6)) +
  labs(x = 'Ranked OTUs', y = 'Occupancy by sample')


# ok now I have core ASVs that are positively and negatively selected for, and I can look at how those vary across functional groups... but how to choose there are so many!

#### Core ASV (16s): Grass ####
grass.core <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "Grass")

otu <- as.data.frame(t(otu_table(grass.core)))

nReads <- 9434

otu_PA <- 1*((otu>0)==1)
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu')


map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
  gather(SampleID.occ, abun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  group_by(otu, FunGroup) %>% 
  summarise(plot_freq = sum(abun > 0)/length(abun), 
            coreTrt = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0)) %>%    
  group_by(otu) %>% 
  summarise(sumF = sum(plot_freq), 
            sumG = sum(coreTrt), 
            nS = length(FunGroup)*2, 
            Index = (sumF + sumG)/nS) 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by = 'otu') %>% 
  transmute(otu = otu,
            rank = Index) %>%
  arrange(desc(rank))

# 
# BCaddition <- NULL
# 
# otu_start <- otu_ranked$otu[1]
# start_matrix <- as.matrix(otu[otu_start,])
# 
# x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads))
# x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
# df_s <- data.frame(x_names,x)
# names(df_s)[2] <- 1
# BCaddition <- rbind(BCaddition,df_s)
# 
# for(i in 2:500){
#   otu_add=otu_ranked$otu[i]
#   add_matrix <- as.matrix(otu[otu_add,])
#   start_matrix <- rbind(start_matrix, add_matrix)
#   x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
#   x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
#   df_a <- data.frame(x_names,x)
#   names(df_a)[2] <- i
#   BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))
# }
# 
# 
# x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
# x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
# df_full <- data.frame(x_names, x)
# names(df_full)[2] <- i + 1
# BCfull <- left_join(BCaddition,df_full, by = 'x_names')
# 
# rownames(BCfull) <- BCfull$x_names
# temp_BC <- BCfull
# temp_BC$x_names <- NULL
# temp_BC_matrix <- as.matrix(temp_BC)
# 
# BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>%
#   gather(comparison, BC, -rank) %>%
#   group_by(rank) %>%
#   summarise(MeanBC=mean(BC)) %>%
#   arrange(-dplyr::desc(MeanBC)) %>%
#   mutate(proportionBC=MeanBC/max(MeanBC))
# 
# Increase <- BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
# increaseDF <- data.frame(IncreaseBC = c(0, (Increase)), rank = factor(c(1:(length(Increase) + 1))))
# BC_ranked <- left_join(BC_ranked, increaseDF)
# BC_ranked <- BC_ranked[-nrow(BC_ranked),]
# 
# rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)
# 
# BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
# 
# elbow <- which.max(BC_ranked$fo_diffs)
# 
# (elbow.g.p <- ggplot(BC_ranked[1:350,],
#                        aes(x = factor(BC_ranked$rank[1:350], levels = BC_ranked$rank[1:450]))) +
#   geom_point(aes(y = proportionBC)) +
#   theme_classic() +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(size = 7, angle = 45)) +
#   geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
#   geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),
#              lty = 3, col = 'blue', cex = .5) +
#   labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
#   annotate(geom = "text",
#            x = elbow + 10,
#            y = .15,
#            label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +
#   annotate(geom = "text",
#            x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4,
#            y = .08,
#            label = paste("Last 2% increase (", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))
# 
# occ_abun$fill <- 'no'
# occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)]))]] <- 'core'
# 
# saveRDS(occ_abun, "Data/16S/Core/occ_abun_asv_g.RDS")
# saveRDS(BC_ranked, "Data/16S/Core/BC_ranked_asv_g.RDS")

occ_abun.asv.g <- readRDS("Data/16S/Core/occ_abun_asv_g.RDS")
BC_ranked.asv.g <- readRDS("Data/16S/Core/BC_ranked_asv_g.RDS")

taxon <- as.matrix(rownames(otu))
rownames(taxon) <- taxon
spp <- t(otu)

obs.np = sncm.fit(spp, taxon, stats = FALSE, pool = NULL)

(core.g.p <- ggplot() +
  geom_point(data = occ_abun.asv.g[occ_abun.asv.g$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun.asv.g[occ_abun.asv.g$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: Grass Samples", x = "log10(mean relative abundance)", y = "Occupancy"))

core.g <- occ_abun.asv.g$otu[occ_abun.asv.g$fill == 'core']

otu_relabun.asv.g <- decostand(otu, method = "total", MARGIN = 2)

## divide core into neutral vs deterministic vs dispersal limited ##
otu.mod.g <- data.frame(otu = occ_abun.asv.g[occ_abun.asv.g$fill == "core",]$otu, mod = rep(NA, length(unique(occ_abun.asv.g[occ_abun.asv.g$fill == "core",]$otu))))

for(i in otu.mod.g$otu){
otu.mod.g[otu.mod.g$otu == i,]$mod <- ifelse(log10(occ_abun.asv.g[occ_abun.asv.g$otu == i,]$otu_occ) > log10(obs.np[obs.np$V1 == i,]$pred.upr), "det", "neu")

otu.mod.g[otu.mod.g$otu == i,]$mod <- ifelse(log10(occ_abun.asv.g[occ_abun.asv.g$otu == i,]$otu_occ) < log10(obs.np[obs.np$V1 == i,]$pred.lwr), "dis", otu.mod.g[otu.mod.g$otu == i,]$mod)

}

rm(occ_abun, obs.np, otu, otu_PA, otu_ranked, spp, taxon, otu_occ, otu_rel, i, nReads, PresenceSum, grass.core)

#### Core ASV (16s): Forb ####
forb.core <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "Forb")

otu <- as.data.frame(t(otu_table(forb.core)))

nReads <- 9434 

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 


map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
  gather(SampleID.occ, abun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  group_by(otu, FunGroup) %>%
  summarise(plot_freq = sum(abun > 0)/length(abun), 
            coreTrt = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0)) %>%    
  group_by(otu) %>% # group by OTU
  summarise(sumF = sum(plot_freq), 
            sumG = sum(coreTrt), 
            nS = length(FunGroup)*2, 
            Index = (sumF + sumG)/nS) 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by = 'otu') %>% 
  transmute(otu = otu,
            rank = Index) %>%
  arrange(dplyr::desc(rank))

# 
# BCaddition <- NULL
# 
# otu_start <- otu_ranked$otu[1] 
# start_matrix <- as.matrix(otu[otu_start,])
# 
# x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) 
# x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
# df_s <- data.frame(x_names,x)
# names(df_s)[2] <- 1 
# BCaddition <- rbind(BCaddition,df_s)
#   
# for(i in 2:500){
#   otu_add=otu_ranked$otu[i] 
#   add_matrix <- as.matrix(otu[otu_add,])
#   start_matrix <- rbind(start_matrix, add_matrix) 
#   x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads)) 
#   x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
#   df_a <- data.frame(x_names,x)
#   names(df_a)[2] <- i 
#   BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))  
# }
# 
# 
# x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
# x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
# df_full <- data.frame(x_names, x)
# names(df_full)[2] <- i + 1 
# BCfull <- left_join(BCaddition,df_full, by = 'x_names') 
#  
# rownames(BCfull) <- BCfull$x_names
# temp_BC <- BCfull
# temp_BC$x_names <- NULL
# temp_BC_matrix <- as.matrix(temp_BC)
# 
# BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
#   gather(comparison, BC, -rank) %>%
#   group_by(rank) %>%
#   summarise(MeanBC=mean(BC)) %>%
#   arrange(-dplyr::desc(MeanBC)) %>%
#   mutate(proportionBC=MeanBC/max(MeanBC))
# 
# Increase <- BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
# increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase)+1))))
# BC_ranked <- left_join(BC_ranked, increaseDF)
# BC_ranked <- BC_ranked[-nrow(BC_ranked),]
# 
# rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)
# 
# BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
# 
# elbow <- which.max(BC_ranked$fo_diffs)
# 
# (elbow.f.p <- ggplot(BC_ranked[1:250,], 
#                        aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
#   geom_point(aes(y = proportionBC)) +
#   theme_classic() + 
#   theme(strip.background = element_blank(), 
#         axis.text.x = element_text(size = 7, angle = 45)) +
#   geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
#   geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])), 
#              lty = 3, col = 'blue', cex = .5) +
#   labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
#   annotate(geom = "text", 
#            x = elbow + 50, 
#            y = .15, 
#            label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +    
#   annotate(geom = "text", 
#            x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 50, 
#            y = .08, 
#            label = paste("Last 2% increase (", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))
# 
# occ_abun$fill <- 'no'
# occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)]))]] <- 'core'
# 
# saveRDS(occ_abun, "Data/16S/Core/occ_abun_asv_f.RDS")
# saveRDS(BC_ranked, "Data/16S/Core/BC_ranked_asv_f.RDS")

occ_abun.asv.f <- readRDS("Data/16S/Core/occ_abun_asv_f.RDS")
BC_ranked.asv.f <- readRDS("Data/16S/Core/BC_ranked_asv_f.RDS")

taxon <- as.matrix(rownames(otu))
rownames(taxon) <- taxon
spp <- t(otu)
obs.np = sncm.fit(spp, taxon, stats = FALSE, pool = NULL)

(core.f.p <- ggplot() +
  geom_point(data = occ_abun.asv.f[occ_abun.asv.f$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun.asv.f[occ_abun.asv.f$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: Forb Samples", x = "log10(mean relative abundance)", y = "Occupancy"))

core.f <- occ_abun.asv.f$otu[occ_abun.asv.f$fill == 'core']

otu_relabun.asv.f <- decostand(otu, method = "total", MARGIN = 2)


## divide core into neutral vs deterministic vs dispersal limited ##
otu.mod.f <- data.frame(otu = occ_abun.asv.f[occ_abun.asv.f$fill == "core",]$otu, mod = rep(NA, length(unique(occ_abun.asv.f[occ_abun.asv.f$fill == "core",]$otu))))

for(i in otu.mod.f$otu){
otu.mod.f[otu.mod.f$otu == i,]$mod <- ifelse(log10(occ_abun.asv.f[occ_abun.asv.f$otu == i,]$otu_occ) > log10(obs.np[obs.np$V1 == i,]$pred.upr), "det", "neu")

otu.mod.f[otu.mod.f$otu == i,]$mod <- ifelse(log10(occ_abun.asv.f[occ_abun.asv.f$otu == i,]$otu_occ) < log10(obs.np[obs.np$V1 == i,]$pred.lwr), "dis", otu.mod.f[otu.mod.f$otu == i,]$mod)

}

rm(occ_abun, obs.np, otu, otu_PA, otu_ranked, spp, taxon, otu_occ, otu_rel, i, nReads, PresenceSum, forb.core)

#### Core ASV (16s): Grass x Forb ####
gf.core <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "grass_x_forb")

otu <- as.data.frame(t(otu_table(gf.core)))

nReads <- 9434 

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 

map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
  gather(SampleID.occ, abun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  group_by(otu, FunGroup) %>% 
  summarise(plot_freq = sum(abun > 0)/length(abun), 
            coreTrt = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0)) %>%    
  group_by(otu) %>% 
  summarise(sumF = sum(plot_freq), 
            sumG = sum(coreTrt), 
            nS = length(FunGroup)*2, 
            Index = (sumF + sumG)/nS) 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by = 'otu') %>% 
  transmute(otu = otu,
            rank = Index) %>%
  arrange(dplyr::desc(rank))

# BCaddition <- NULL
# 
# otu_start <- otu_ranked$otu[1] 
# start_matrix <- as.matrix(otu[otu_start,])
# 
# x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) 
# x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
# df_s <- data.frame(x_names,x)
# names(df_s)[2] <- 1 
# BCaddition <- rbind(BCaddition,df_s)
#   
# for(i in 2:500){
#   otu_add = otu_ranked$otu[i] 
#   add_matrix <- as.matrix(otu[otu_add,])
#   start_matrix <- rbind(start_matrix, add_matrix) 
#   x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads)) 
#   x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
#   df_a <- data.frame(x_names,x)
#   names(df_a)[2] <- i 
#   BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))  
# }
# 
# 
# x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
# x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
# df_full <- data.frame(x_names, x)
# names(df_full)[2] <- i + 1 
# BCfull <- left_join(BCaddition, df_full, by = 'x_names') 
#  
# rownames(BCfull) <- BCfull$x_names
# temp_BC <- BCfull
# temp_BC$x_names <- NULL
# temp_BC_matrix <- as.matrix(temp_BC)
# 
# BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>% 
#   gather(comparison, BC, -rank) %>%
#   group_by(rank) %>%
#   summarise(MeanBC = mean(BC)) %>%
#   arrange(-dplyr::desc(MeanBC)) %>%
#   mutate(proportionBC = MeanBC/max(MeanBC))
# 
# Increase = BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
# increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase) + 1))))
# BC_ranked <- left_join(BC_ranked, increaseDF)
# BC_ranked <- BC_ranked[-nrow(BC_ranked),]
# 
# rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)
# 
# BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
# 
# elbow <- which.max(BC_ranked$fo_diffs)
# 
# (elbow.gf.p <- ggplot(BC_ranked[1:250,], 
#                        aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
#   geom_point(aes(y = proportionBC)) +
#   theme_classic() + 
#   theme(strip.background = element_blank(), 
#         axis.text.x = element_text(size = 7, angle = 45)) +
#   geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
#   geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])), 
#              lty = 3, col = 'blue', cex = .5) +
#   labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
#   annotate(geom = "text", 
#            x = elbow + 10, 
#            y = .15, 
#            label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +    
#   annotate(geom = "text", 
#            x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4, 
#            y = .08, 
#            label = paste("Last 2% increase (", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))
# 
# occ_abun$fill <- 'no'
# occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)]))]] <- 'core'
# 
# saveRDS(occ_abun, "Data/16S/Core/occ_abun_asv_gf.RDS")
# saveRDS(BC_ranked, "Data/16S/Core/BC_ranked_asv_gf.RDS")

occ_abun.asv.gf <- readRDS("Data/16S/Core/occ_abun_asv_gf.RDS")
BC_ranked.asv.gf <- readRDS("Data/16S/Core/BC_ranked_asv_gf.RDS")

taxon <- as.matrix(rownames(otu))
rownames(taxon) <- taxon
spp <- t(otu)

obs.np = sncm.fit(spp, taxon, stats = FALSE, pool = NULL)

(core.gf.p <- ggplot() +
  geom_point(data = occ_abun.asv.gf[occ_abun.asv.gf$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun.asv.gf[occ_abun.asv.gf$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: Grass + Forb", x = "log10(mean relative abundance)", y = "Occupancy"))

core.gf <- occ_abun.asv.gf$otu[occ_abun.asv.gf$fill == 'core']

otu_relabun.asv.gf <- decostand(otu, method = "total", MARGIN = 2)


## divide core into neutral vs deterministic vs dispersal limited ##
otu.mod.fg <- data.frame(otu = occ_abun.asv.gf[occ_abun.asv.gf$fill == "core",]$otu, mod = rep(NA, length(unique(occ_abun.asv.gf[occ_abun.asv.gf$fill == "core",]$otu))))

for(i in otu.mod.fg$otu){
otu.mod.fg[otu.mod.fg$otu == i,]$mod <- ifelse(log10(occ_abun.asv.gf[occ_abun.asv.gf$otu == i,]$otu_occ) > log10(obs.np[obs.np$V1 == i,]$pred.upr), "det", "neu")

otu.mod.fg[otu.mod.fg$otu == i,]$mod <- ifelse(log10(occ_abun.asv.gf[occ_abun.asv.gf$otu == i,]$otu_occ) < log10(obs.np[obs.np$V1 == i,]$pred.lwr), "dis", otu.mod.fg[otu.mod.fg$otu == i,]$mod)

}

rm(occ_abun, obs.np, otu, otu_PA, otu_ranked, spp, taxon, otu_occ, otu_rel, i, nReads, PresenceSum, gf.core)

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

# what about those that are deterministically chosen by the plant?
otu.mod.f.det <- filter(otu.mod.f, mod == "det")
otu.mod.g.det <- filter(otu.mod.g, mod == "det")
otu.mod.fg.det <- filter(otu.mod.fg, mod == "det")

nrow(otu.mod.f.det[otu.mod.f.det$otu %in% otu.mod.g.det$otu,]) # n12: 67 forb OTUs also in grass
nrow(otu.mod.g.det[otu.mod.g.det$otu %in% otu.mod.fg.det$otu,]) # n23: 70 grass OTUs in forbsxgrasses
nrow(otu.mod.f.det[otu.mod.f.det$otu %in% otu.mod.fg.det$otu,]) # n13: 59 forb OTUs in forbsxgrasses

test <- otu.mod.g.det[otu.mod.g.det$otu %in% otu.mod.fg.det$otu,]
nrow(otu.mod.f.det[otu.mod.f.det$otu %in% test$otu,]) # n123 45

nrow(otu.mod.f[otu.mod.f$mod == "det",]) # 93
nrow(otu.mod.g[otu.mod.g$mod == "det",]) #171
nrow(otu.mod.fg[otu.mod.fg$mod == "det",]) #106

grid.newpage()  
draw.triple.venn(area1 = 93,   # Forb core                  
                 area2 = 171,   # Grass core
                 area3 = 106,   # g_x_f core
                 n12 = 67,
                 n23 = 70,
                 n13 = 59,
                 n123 = 45,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Forbs", "Grasses", "Forbs & Grasses"))


#### Core Family (16s): All ####

# There is probably a much cleaner, fast to way to manipulate these files into the right format, but whatever 
ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion))) # these were likely removed earlier but just to be sure
core.all.fam <- tax_glom(ps.16s.nocontrols.rare, taxrank = "Family", NArm = FALSE)

otu <- data.frame(otu = rownames(t(otu_table(core.all.fam))), t(otu_table(core.all.fam)))
tax <- t(data.frame(t(tax_table(core.all.fam))))
tax <- data.frame(otu = rownames(tax), tax)

tax$ord.fam <- paste(tax$Order, tax$Family, sep = ".")

otu <- merge(tax[,c(1,9)], otu, by = "otu")

rownames(otu) <- otu$ord.fam
otu <- otu[ , c(3:ncol(otu))]

nReads <- 9434 
map <- mapping

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)   
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 

map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
  gather(SampleID.occ, abun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  group_by(otu, FunGroup) %>% 
  summarise(plot_freq = sum(abun > 0)/length(abun), 
            coreTrt = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0)) %>%    
  group_by(otu) %>% 
  summarise(sumF = sum(plot_freq), 
            sumG = sum(coreTrt), 
            nS = length(FunGroup)*2,
            Index = (sumF + sumG)/nS) 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by = 'otu') %>% 
  transmute(otu = otu,
            rank = Index) %>%
  arrange(dplyr::desc(rank))
# 
# BCaddition <- NULL
# 
# otu_start <- otu_ranked$otu[1] 
# start_matrix <- as.matrix(otu[otu_start,]) 
# 
# x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) 
# x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
# df_s <- data.frame(x_names,x)
# names(df_s)[2] <- 1 
# BCaddition <- rbind(BCaddition,df_s)
#   
# for(i in 2:299){
#   otu_add = otu_ranked$otu[i] 
#   add_matrix <- as.matrix(otu[otu_add,])
#   start_matrix <- rbind(start_matrix, add_matrix) 
#   x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads)) 
#   x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
#   df_a <- data.frame(x_names, x)
#   names(df_a)[2] <- i 
#   BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))  
# }
# 
#  
# x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]] - otu[,x[2]]))/(2*nReads))
# x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse = ' - '))
# df_full <- data.frame(x_names, x)
# 
# names(df_full)[2] <- i + 1 
# BCfull <- left_join(BCaddition, df_full, by = 'x_names') 
#  
# rownames(BCfull) <- BCfull$x_names
# temp_BC <- BCfull
# temp_BC$x_names <- NULL
# temp_BC_matrix <- as.matrix(temp_BC)
# 
# BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>% 
#   gather(comparison, BC, -rank) %>%
#   group_by(rank) %>%
#   summarise(MeanBC = mean(BC)) %>%
#   arrange(-dplyr::desc(MeanBC)) %>%
#   mutate(proportionBC = MeanBC/max(MeanBC))
# 
# Increase = BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
# increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase) + 1))))
# BC_ranked <- left_join(BC_ranked, increaseDF)
# BC_ranked <- BC_ranked[-nrow(BC_ranked),]
# 
# rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)
# 
# BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
# 
# elbow <- which.max(BC_ranked$fo_diffs)
# 
# (elbow.all.p <- ggplot(BC_ranked[1:250,], aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
#   geom_point(aes(y = proportionBC)) +
#   theme_classic() + 
#   theme(strip.background = element_blank(), 
#         axis.text.x = element_text(size = 7, angle = 45)) +
#   geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
#   geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])), 
#              lty = 3, col = 'blue', cex = .5) +
#   labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
#   annotate(geom = "text", 
#            x = elbow + 10, 
#            y = .15, 
#            label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +    
#   annotate(geom = "text", 
#            x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4, 
#            y = .08, 
#            label = paste("Last 2% increase (",last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))
# 
# occ_abun$fill <- 'no'
# occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))]] <- 'core'
# 
# 
# saveRDS(occ_abun, "Data/16S/Core/occ_abun_fam_all.RDS")
# saveRDS(BC_ranked, "Data/16S/Core/BC_ranked_fam_all.RDS")

occ_abun.fam.all <- readRDS("Data/16S/Core/occ_abun_fam_all.RDS")
BC_ranked.fam.all <- readRDS("Data/16S/Core/BC_ranked_fam_all.RDS")

taxon <- as.matrix(rownames(otu))
rownames(taxon) <- taxon
spp <- t(otu)

obs.np <- sncm.fit(spp, taxon, stats = FALSE, pool = NULL)

(core.all.p <- ggplot() +
  geom_point(data = occ_abun.fam.all[occ_abun.fam.all$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun.fam.all[occ_abun.fam.all$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: All Samples", x = "log10(mean relative abundance)", y = "Occupancy"))

core.all <- occ_abun.fam.all$otu[occ_abun.fam.all$fill == 'core']

otu_relabun.fam.all <- decostand(otu, method = "total", MARGIN = 2)

plotDF <-  data.frame(otu = as.factor(row.names(otu_relabun.fam.all)), otu_relabun.fam.all) %>% 
  gather(SampleID.occ, relabun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  left_join(otu_ranked, by = 'otu') %>%
  filter(otu %in% core.all) %>% 
  group_by(FunGroup, otu) %>%
  summarise(plot_freq = sum(relabun>0)/length(relabun), 
            coreSite = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0))

plotDF$otu <- factor(plotDF$otu, levels = otu_ranked$otu[1:78]) 

(all.plotDF <- ggplot(plotDF, aes(x = otu, plot_freq, group = FunGroup, fill = FunGroup)) +    
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(plotDF$otu %in% core.all))) +
  theme(axis.text = element_text(size = 6)) +
  labs(x = 'Ranked OTUs', y = 'Occupancy by site')) # everything is in like 90% of samples

otu_RA_core.all <- data.frame(otu = as.factor(row.names(otu_relabun.fam.all)), otu_relabun.fam.all) %>%
  filter(otu %in% plotDF$otu)

## divide core into neutral vs deterministic vs dispersal limited ##
otu.mod.all <- data.frame(otu = occ_abun.fam.all[occ_abun.fam.all$fill == "core",]$otu, mod = rep(NA, length(unique(occ_abun.fam.all[occ_abun.fam.all$fill == "core",]$otu))))

for(i in otu.mod.all$otu){
otu.mod.all[otu.mod.all$otu == i,]$mod <- ifelse(log10(occ_abun.fam.all[occ_abun.fam.all$otu == i,]$otu_occ) > log10(obs.np[obs.np$V1 == i,]$pred.upr), "det", "neu")

otu.mod.all[otu.mod.all$otu == i,]$mod <- ifelse(log10(occ_abun.fam.all[occ_abun.fam.all$otu == i,]$otu_occ) < log10(obs.np[obs.np$V1 == i,]$pred.lwr), "dis", otu.mod.all[otu.mod.all$otu == i,]$mod)

}

#### Core Family (16s): Grass ####
grass.core <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "Grass")

core.g.ord <- tax_glom(grass.core, taxrank = "Family", NArm = FALSE)

otu <- data.frame(otu = rownames(t(otu_table(core.g.ord))), t(otu_table(core.g.ord)))
tax <- t(data.frame(t(tax_table(core.g.ord))))
tax <- data.frame(otu = rownames(tax), tax)

tax$ord.fam <- paste(tax$Order, tax$Family, sep = ".")

otu <- merge(tax[,c(1,9)], otu, by = "otu")

rownames(otu) <- otu$ord.fam
otu <- otu[ , c(3:ncol(otu))]

nReads <- 9434 

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 

map <- mapping
map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)


PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
  gather(SampleID.occ, abun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  group_by(otu, FunGroup) %>% 
  summarise(plot_freq = sum(abun > 0)/length(abun), 
            coreTrt = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0)) %>%    
  group_by(otu) %>% 
  summarise(sumF = sum(plot_freq), 
            sumG = sum(coreTrt), 
            nS = length(FunGroup)*2, 
            Index = (sumF + sumG)/nS) 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by = 'otu') %>% # 
  transmute(otu = otu,
            rank = Index) %>%
  arrange(dplyr::desc(rank))


# BCaddition <- NULL
# 
# otu_start <- otu_ranked$otu[1]
# start_matrix <- as.matrix(otu[otu_start,])
# 
# x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads))
# x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
# df_s <- data.frame(x_names,x)
# names(df_s)[2] <- 1
# BCaddition <- rbind(BCaddition,df_s)
# 
# for(i in 2:299){
#   otu_add=otu_ranked$otu[i]
#   add_matrix <- as.matrix(otu[otu_add,])
#   start_matrix <- rbind(start_matrix, add_matrix)
#   x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
#   x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
#   df_a <- data.frame(x_names,x)
#   names(df_a)[2] <- i
#   BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))
# }
# 
# 
# x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
# x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
# df_full <- data.frame(x_names, x)
# names(df_full)[2] <- i + 1
# BCfull <- left_join(BCaddition,df_full, by = 'x_names')
# 
# rownames(BCfull) <- BCfull$x_names
# temp_BC <- BCfull
# temp_BC$x_names <- NULL
# temp_BC_matrix <- as.matrix(temp_BC)
# 
# BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>%
#   gather(comparison, BC, -rank) %>%
#   group_by(rank) %>%
#   summarise(MeanBC=mean(BC)) %>%
#   arrange(-dplyr::desc(MeanBC)) %>%
#   mutate(proportionBC=MeanBC/max(MeanBC))
# 
# Increase <- BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
# increaseDF <- data.frame(IncreaseBC = c(0, (Increase)), rank = factor(c(1:(length(Increase) + 1))))
# BC_ranked <- left_join(BC_ranked, increaseDF)
# BC_ranked <- BC_ranked[-nrow(BC_ranked),]
# 
# rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)
# 
# BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
# 
# elbow <- which.max(BC_ranked$fo_diffs)
# 
# (elbow.g.p <- ggplot(BC_ranked[1:350,],
#                        aes(x = factor(BC_ranked$rank[1:350], levels = BC_ranked$rank[1:450]))) +
#   geom_point(aes(y = proportionBC)) +
#   theme_classic() +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(size = 7, angle = 45)) +
#   geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
#   geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),
#              lty = 3, col = 'blue', cex = .5) +
#   labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
#   annotate(geom = "text",
#            x = elbow + 10,
#            y = .15,
#            label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +
#   annotate(geom = "text",
#            x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4,
#            y = .08,
#            label = paste("Last 2% increase (", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))
# 
# occ_abun$fill <- 'no'
# occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)]))]] <- 'core'
# 
# saveRDS(occ_abun, "Data/16S/Core/occ_abun_fam_g.RDS")
# saveRDS(BC_ranked, "Data/16S/Core/BC_ranked_fam_g.RDS")

occ_abun.fam.g <- readRDS("Data/16S/Core/occ_abun_fam_g.RDS")
BC_ranked.fam.g <- readRDS("Data/16S/Core/BC_ranked_fam_g.RDS")

spp <- t(otu)
taxon <- as.matrix(rownames(otu))
rownames(taxon) <- taxon
obs.np = sncm.fit(spp, taxon, stats=FALSE, pool = NULL)
#sta.np = sncm.fit(t(otu), as.vector(rownames(otu)), stats = TRUE, pool = NULL)

(core.g.p <- ggplot() +
  geom_point(data = occ_abun.fam.g[occ_abun.fam.g$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun.fam.g[occ_abun.fam.g$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Fam Core: Grasses", x = "log10(mean relative abundance)", y = "Occupancy"))

core.g <- occ_abun.fam.g$otu[occ_abun.fam.g$fill == 'core']

otu_relabun.fam.g <- decostand(otu, method = "total", MARGIN = 2)

#### Core Family (16s): Forb ####
forb.core <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "Forb")
core.f.fam <- tax_glom(forb.core, taxrank = "Family", NArm = FALSE)

otu <- data.frame(otu = rownames(t(otu_table(core.f.fam))), t(otu_table(core.f.fam)))
tax <- t(data.frame(t(tax_table(core.f.fam))))
tax <- data.frame(otu = rownames(tax), tax)

tax$ord.fam <- paste(tax$Order, tax$Family, sep = ".")

otu <- merge(tax[,c(1,9)], otu, by = "otu")

rownames(otu) <- otu$ord.fam
otu <- otu[ , c(3:ncol(otu))]

nReads <- 9434 

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 

map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
  gather(SampleID.occ, abun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  group_by(otu, FunGroup) %>% 
  summarise(plot_freq = sum(abun > 0)/length(abun), 
            coreTrt = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0)) %>%
  group_by(otu) %>% # group by OTU
  summarise(sumF = sum(plot_freq),
            sumG = sum(coreTrt), 
            nS = length(FunGroup)*2, 
            Index = (sumF + sumG)/nS) 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by = 'otu') %>% 
  transmute(otu = otu,
            rank = Index) %>%
  arrange(dplyr::desc(rank))


# BCaddition <- NULL
# 
# otu_start <- otu_ranked$otu[1] 
# start_matrix <- as.matrix(otu[otu_start,])
# 
# x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) 
# x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
# df_s <- data.frame(x_names,x)
# names(df_s)[2] <- 1 
# BCaddition <- rbind(BCaddition,df_s)
#   
# for(i in 2:299){
#   otu_add=otu_ranked$otu[i] 
#   add_matrix <- as.matrix(otu[otu_add,])
#   start_matrix <- rbind(start_matrix, add_matrix) 
#   x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads)) 
#   x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
#   df_a <- data.frame(x_names,x)
#   names(df_a)[2] <- i 
#   BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))  
# }
# 
# 
# x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
# x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
# df_full <- data.frame(x_names, x)
# names(df_full)[2] <- i + 1 
# BCfull <- left_join(BCaddition,df_full, by = 'x_names') 
#  
# rownames(BCfull) <- BCfull$x_names
# temp_BC <- BCfull
# temp_BC$x_names <- NULL
# temp_BC_matrix <- as.matrix(temp_BC)
# 
# BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
#   gather(comparison, BC, -rank) %>%
#   group_by(rank) %>%
#   summarise(MeanBC=mean(BC)) %>%
#   arrange(-dplyr::desc(MeanBC)) %>%
#   mutate(proportionBC=MeanBC/max(MeanBC))
# 
# Increase <- BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
# increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase)+1))))
# BC_ranked <- left_join(BC_ranked, increaseDF)
# BC_ranked <- BC_ranked[-nrow(BC_ranked),]
# 
# rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)
# 
# BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
# 
# elbow <- which.max(BC_ranked$fo_diffs)
# 
# (elbow.f.p <- ggplot(BC_ranked[1:250,], 
#                        aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
#   geom_point(aes(y = proportionBC)) +
#   theme_classic() + 
#   theme(strip.background = element_blank(), 
#         axis.text.x = element_text(size = 7, angle = 45)) +
#   geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
#   geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])), 
#              lty = 3, col = 'blue', cex = .5) +
#   labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
#   annotate(geom = "text", 
#            x = elbow + 50, 
#            y = .15, 
#            label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +    
#   annotate(geom = "text", 
#            x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 50, 
#            y = .08, 
#            label = paste("Last 2% increase (", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))
# 
# occ_abun$fill <- 'no'
# occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)]))]] <- 'core'
# 
# saveRDS(occ_abun, "Data/16S/Core/occ_abun_fam_f.RDS")
# saveRDS(BC_ranked, "Data/16S/Core/BC_ranked_fam_f.RDS")

occ_abun.fam.f <- readRDS("Data/16S/Core/occ_abun_fam_f.RDS")
BC_ranked.fam.f <- readRDS("Data/16S/Core/BC_ranked_fam_f.RDS")

spp <- t(otu)
taxon <- as.matrix(rownames(otu))
rownames(taxon) <- taxon
obs.np <- sncm.fit(spp, taxon, stats = FALSE, pool = NULL)

(core.f.p <- ggplot() +
  geom_point(data = occ_abun.fam.f[occ_abun.fam.f$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun.fam.f[occ_abun.fam.f$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Family Core: Forbs", x = "log10(mean relative abundance)", y = "Occupancy"))

core.f <- occ_abun.fam.f$otu[occ_abun.fam.f$fill == 'core']

otu_relabun.fam.f <- decostand(otu, method = "total", MARGIN = 2)

otu.mod.f <- data.frame(otu = occ_abun.fam.f[occ_abun.fam.f$fill == "core",]$otu, mod = rep(NA, 82))

for(i in otu.mod.f$otu){
otu.mod.f[otu.mod.f$otu == i,]$mod <- ifelse(log10(occ_abun.fam.f[occ_abun.fam.f$otu == i,]$otu_occ) > log10(obs.np[obs.np$V1 == i,]$pred.upr), "det", "neu")

otu.mod.f[otu.mod.f$otu == i,]$mod <- ifelse(log10(occ_abun.fam.f[occ_abun.fam.f$otu == i,]$otu_occ) < log10(obs.np[obs.np$V1 == i,]$pred.lwr), "dis", otu.mod.f[otu.mod.f$otu == i,]$mod)

}

#### Core Family (16s): Grass x Forb ####
gf.core <- subset_samples(ps.16s.nocontrols.rare, FunGroup == "grass_x_forb")
core.gf.fam <- tax_glom(gf.core, taxrank = "Family", NArm = FALSE)

otu <- data.frame(otu = rownames(t(otu_table(core.gf.fam))), t(otu_table(core.gf.fam)))
tax <- t(data.frame(t(tax_table(core.gf.fam))))
tax <- data.frame(otu = rownames(tax), tax)

tax$ord.fam <- paste(tax$Order, tax$Family, sep = ".")

otu <- merge(tax[,c(1,9)], otu, by = "otu")

rownames(otu) <- otu$ord.fam
otu <- otu[ , c(3:ncol(otu))]

nReads <- 9434 

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 

map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
  gather(SampleID.occ, abun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  group_by(otu, FunGroup) %>% 
  summarise(plot_freq = sum(abun > 0)/length(abun), 
            coreTrt = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0)) %>%    
  group_by(otu) %>% 
  summarise(sumF = sum(plot_freq), 
            sumG = sum(coreTrt), 
            nS = length(FunGroup)*2, 
            Index = (sumF + sumG)/nS) 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by = 'otu') %>% 
  transmute(otu = otu,
            rank = Index) %>%
  arrange(dplyr::desc(rank))


# BCaddition <- NULL
# 
# otu_start <- otu_ranked$otu[1]
# start_matrix <- as.matrix(otu[otu_start,])
# 
# x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads))
# x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
# df_s <- data.frame(x_names,x)
# names(df_s)[2] <- 1
# BCaddition <- rbind(BCaddition,df_s)
# 
# for(i in 2:299){
#   otu_add = otu_ranked$otu[i]
#   add_matrix <- as.matrix(otu[otu_add,])
#   start_matrix <- rbind(start_matrix, add_matrix)
#   x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
#   x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
#   df_a <- data.frame(x_names,x)
#   names(df_a)[2] <- i
#   BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))
# }
# 
# 
# x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
# x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
# df_full <- data.frame(x_names, x)
# names(df_full)[2] <- i + 1
# BCfull <- left_join(BCaddition, df_full, by = 'x_names')
# 
# rownames(BCfull) <- BCfull$x_names
# temp_BC <- BCfull
# temp_BC$x_names <- NULL
# temp_BC_matrix <- as.matrix(temp_BC)
# 
# BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>%
#   gather(comparison, BC, -rank) %>%
#   group_by(rank) %>%
#   summarise(MeanBC = mean(BC)) %>%
#   arrange(-dplyr::desc(MeanBC)) %>%
#   mutate(proportionBC = MeanBC/max(MeanBC))
# 
# Increase = BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
# increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase) + 1))))
# BC_ranked <- left_join(BC_ranked, increaseDF)
# BC_ranked <- BC_ranked[-nrow(BC_ranked),]
# 
# rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)
# 
# BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
# 
# elbow <- which.max(BC_ranked$fo_diffs)
# 
# (elbow.gf.p <- ggplot(BC_ranked[1:250,],
#                        aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
#   geom_point(aes(y = proportionBC)) +
#   theme_classic() +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(size = 7, angle = 45)) +
#   geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
#   geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),
#              lty = 3, col = 'blue', cex = .5) +
#   labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
#   annotate(geom = "text",
#            x = elbow + 10,
#            y = .15,
#            label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +
#   annotate(geom = "text",
#            x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4,
#            y = .08,
#            label = paste("Last 2% increase (", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))
# 
# occ_abun$fill <- 'no'
# occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)]))]] <- 'core'
# 
# saveRDS(occ_abun, "Data/16S/Core/occ_abun_fam_gf.RDS")
# saveRDS(BC_ranked, "Data/16S/Core/BC_ranked_fam_gf.RDS")

occ_abun.fam.gf <- readRDS("Data/16S/Core/occ_abun_fam_gf.RDS")
BC_ranked.fam.gf <- readRDS("Data/16S/Core/BC_ranked_fam_gf.RDS")

spp <- t(otu)
taxon <- as.matrix(rownames(otu))
rownames(taxon) <- taxon
obs.np <- sncm.fit(spp, taxon, stats = FALSE, pool = NULL)

(core.gf.p <- ggplot() +
  geom_point(data = occ_abun.fam.gf[occ_abun.fam.gf$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun.fam.gf[occ_abun.fam.gf$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Fam Core: grass x forb", x = "log10(mean relative abundance)", y = "Occupancy"))

core.gf <- occ_abun.fam.gf$otu[occ_abun.fam.gf$fill == 'core']

otu_relabun.fam.gf <- decostand(otu, method = "total", MARGIN = 2)

#### Comparison of Core ####
library(VennDiagram)
core.g <- data.frame(FunGroup = "Grass", OTU = core.g)
core.f <- data.frame(FunGroup = "Forb", OTU = core.f)
core.gf <- data.frame(FunGroup = "grass_x_forb", OTU = core.gf)
core.all <- data.frame(FunGroup = "All", OTU = core.all)

# forb = 1
# grass = 2
# grass x forb = 3

nrow(core.f[core.f$OTU %in% core.g$OTU,]) # n12: 61 forb OTUs also in grass
nrow(core.g[core.g$OTU %in% core.gf$OTU,]) # n23: 48 grass OTUs in forbsxgrasses
nrow(core.f[core.f$OTU %in% core.gf$OTU,]) # n13: 52 forb OTUs in forbsxgrasses

test <- core.g[core.g$OTU %in% core.gf$OTU,]
nrow(core.f[core.f$OTU %in% test$OTU,]) # n123 46

grid.newpage()  
draw.triple.venn(area1 = 82,   # Forb core                  
                 area2 = 80,   # Grass core
                 area3 = 55,   # g_x_f core
                 n12 = 61,
                 n23 = 48,
                 n13 = 52,
                 n123 = 46,
                 fill = c("pink", "green", "orange"),
                 lty = "blank",
                 category = c("Forbs", "Grasses", "Forbs & Grasses"))


# what are the new families that don't show up in grasses or forbs?
core.gf[!core.gf$OTU %in% core.f$OTU,]

#### Core ASV (ITS): All ####
ps.ITS.nocontrols.rare <- subset_samples(ps.ITS.nocontrols.rare, !(is.na(Competion))) # these were likely removed earlier but just to be sure
otu <- as.data.frame(t(otu_table(ps.ITS.nocontrols.rare)))

nReads <- 7557
map <- mapping

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)   
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 
map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
  gather(SampleID.occ, abun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  group_by(otu, FunGroup) %>% 
  summarise(plot_freq = sum(abun > 0)/length(abun), 
            coreTrt = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0)) %>%   
  group_by(otu) %>% 
  summarise(sumF = sum(plot_freq), 
            sumG = sum(coreTrt), 
            nS = length(FunGroup)*2, 
            Index = (sumF + sumG)/nS) 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by = 'otu') %>% 
  transmute(otu = otu,
            rank = Index) %>%
  arrange(dplyr::desc(rank))
# 
# BCaddition <- NULL
# 
# otu_start <- otu_ranked$otu[1] 
# start_matrix <- as.matrix(otu[otu_start,]) 
# 
# x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads)) 
# x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
# df_s <- data.frame(x_names,x)
# names(df_s)[2] <- 1 
# BCaddition <- rbind(BCaddition,df_s)
# 
# for(i in 2:500){
#   otu_add = otu_ranked$otu[i] 
#   add_matrix <- as.matrix(otu[otu_add,])
#   start_matrix <- rbind(start_matrix, add_matrix) 
#   x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
#   x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
#   df_a <- data.frame(x_names, x)
#   names(df_a)[2] <- i
#   BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))
# }
# 
# 
# x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]] - otu[,x[2]]))/(2*nReads))
# x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse = ' - '))
# df_full <- data.frame(x_names, x)
# 
# names(df_full)[2] <- i + 1
# BCfull <- left_join(BCaddition, df_full, by = 'x_names') 
# 
# rownames(BCfull) <- BCfull$x_names
# temp_BC <- BCfull
# temp_BC$x_names <- NULL
# temp_BC_matrix <- as.matrix(temp_BC)
# 
# BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>%
#   gather(comparison, BC, -rank) %>%
#   group_by(rank) %>%
#   summarise(MeanBC = mean(BC)) %>%
#   arrange(-dplyr::desc(MeanBC)) %>%
#   mutate(proportionBC = MeanBC/max(MeanBC))
# 
# Increase = BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
# increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase) + 1))))
# BC_ranked <- left_join(BC_ranked, increaseDF)
# BC_ranked <- BC_ranked[-nrow(BC_ranked),]
# 
# rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)
# 
# BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
# 
# elbow <- which.max(BC_ranked$fo_diffs)
# 
# (elbow.all.p <- ggplot(BC_ranked[1:250,], aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
#   geom_point(aes(y = proportionBC)) +
#   theme_classic() +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(size = 7, angle = 45)) +
#   geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
#   geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),
#              lty = 3, col = 'blue', cex = .5) +
#   labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
#   annotate(geom = "text",
#            x = elbow + 10,
#            y = .15,
#            label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +
#   annotate(geom = "text",
#            x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4,
#            y = .08,
#            label = paste("Last 2% increase (",last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))
# 
# occ_abun$fill <- 'no'
# occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))]] <- 'core'
# 
# saveRDS(occ_abun, "Data/ITS/Core/occ_abun_asv_all.RDS")
# saveRDS(BC_ranked, "Data/ITS/Core/BC_ranked_asv_all.RDS")

occ_abun.asv.all <- readRDS("Data/ITS/Core/occ_abun_asv_all.RDS")
BC_ranked.asv.all <- readRDS("Data/ITS/Core/BC_ranked_asv_all.RDS")

# Models for the whole community
taxon <- as.matrix(rownames(otu))
rownames(taxon) <- taxon
spp <- t(otu)
obs.np <- sncm.fit(spp = spp, taxon = taxon, stats = FALSE, pool = NULL)

ggplot() +
  geom_point(data = occ_abun.asv.all[occ_abun.asv.all$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun.asv.all[occ_abun.asv.all$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "ITS Core: All Samples", x = "log10(mean relative abundance)", y = "Occupancy")

core.all <- occ_abun.asv.all$otu[occ_abun.asv.all$fill == 'core']

otu_relabun.asv.all <- decostand(otu, method = "total", MARGIN = 2)

plotDF <-  data.frame(otu = as.factor(row.names(otu_relabun.asv.all)), otu_relabun.asv.all) %>%
  gather(SampleID.occ, relabun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  left_join(otu_ranked, by = 'otu') %>%
  filter(otu %in% core.all) %>% 
  group_by(FunGroup, otu) %>%
  summarise(plot_freq = sum(relabun>0)/length(relabun), 
            coreSite = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0))

plotDF$otu <- factor(plotDF$otu, levels = otu_ranked$otu[1:187])

(all.plotDF <- ggplot(plotDF, aes(x = otu, plot_freq, group = FunGroup, fill = FunGroup)) +    
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(plotDF$otu %in% core.all))) +
  theme(axis.text = element_text(size = 6)) +
  labs(x = 'Ranked OTUs', y = 'Occupancy by sample'))

## divide core into neutral vs deterministic vs dispersal limited ##
otu.mod.all <- data.frame(otu = occ_abun.asv.all[occ_abun.asv.all$fill == "core",]$otu, mod = rep(NA, 187))

for(i in otu.mod.all$otu){
otu.mod.all[otu.mod.all$otu == i,]$mod <- ifelse(log10(occ_abun.asv.all[occ_abun.asv.all$otu == i,]$otu_occ) > log10(obs.np[obs.np$V1 == i,]$pred.upr), "det", "neu")

otu.mod.all[otu.mod.all$otu == i,]$mod <- ifelse(log10(occ_abun.asv.all[occ_abun.asv.all$otu == i,]$otu_occ) < log10(obs.np[obs.np$V1 == i,]$pred.lwr), "dis", otu.mod.all[otu.mod.all$otu == i,]$mod)

}

rm(occ_abun, obs.np, otu, otu_PA, otu_ranked, spp, taxon, otu_occ, otu_rel, i, nReads, PresenceSum)

#### Core Family (ITS): All #### 
ps.ITS.nocontrols.rare <- subset_samples(ps.ITS.nocontrols.rare, !(is.na(Competion))) # these were likely removed earlier but just to be sure
core.all.fam <- tax_glom(ps.ITS.nocontrols.rare, taxrank = "Family", NArm = FALSE)

otu <- data.frame(otu = rownames(t(otu_table(core.all.fam))), t(otu_table(core.all.fam)))
tax <- t(data.frame(t(tax_table(core.all.fam))))
tax <- data.frame(otu = rownames(tax), tax)

tax$ord.fam <- paste(tax$Order, tax$Family, sep = ".")

otu <- merge(tax[,c(1,9)], otu, by = "otu")
otu <- filter(otu, ord.fam != "NA.NA")
rownames(otu) <- otu$ord.fam
otu <- otu[ , c(3:ncol(otu))]

nReads <- 7557 
map <- mapping

otu_PA <- 1*((otu>0)==1)  
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)  
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)   
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 

map$SampleID.occ <- gsub("-", "\\.", map$SampleID_Fix)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
  gather(SampleID.occ, abun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  group_by(otu, FunGroup) %>% 
  dplyr::summarise(plot_freq = sum(abun > 0)/length(abun), 
            coreTrt = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0)) %>%    
  group_by(otu) %>% 
  dplyr::summarise(sumF = sum(plot_freq), 
            sumG = sum(coreTrt), 
            nS = length(FunGroup)*2,
            Index = (sumF + sumG)/nS) 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by = 'otu') %>% 
  transmute(otu = otu,
            rank = Index) %>%
  arrange(dplyr::desc(rank))

# BCaddition <- NULL
# 
# otu_start <- otu_ranked$otu[1]
# start_matrix <- as.matrix(otu[otu_start,])
# 
# x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[, x[2]]))/(2 * nReads))
# x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
# df_s <- data.frame(x_names,x)
# names(df_s)[2] <- 1
# BCaddition <- rbind(BCaddition,df_s)
# 
# for(i in 2:255){
#   otu_add = otu_ranked$otu[i]
#   add_matrix <- as.matrix(otu[otu_add,])
#   start_matrix <- rbind(start_matrix, add_matrix)
#   x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
#   x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
#   df_a <- data.frame(x_names, x)
#   names(df_a)[2] <- i
#   BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))
# }
# 
# 
# x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]] - otu[,x[2]]))/(2*nReads))
# x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse = ' - '))
# df_full <- data.frame(x_names, x)
# 
# names(df_full)[2] <- i + 1
# BCfull <- left_join(BCaddition, df_full, by = 'x_names')
# 
# rownames(BCfull) <- BCfull$x_names
# temp_BC <- BCfull
# temp_BC$x_names <- NULL
# temp_BC_matrix <- as.matrix(temp_BC)
# 
# BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))), t(temp_BC_matrix)) %>%
#   gather(comparison, BC, -rank) %>%
#   group_by(rank) %>%
#   dplyr::summarise(MeanBC = mean(BC)) %>%
#   arrange(-dplyr::desc(MeanBC)) %>%
#   mutate(proportionBC = MeanBC/max(MeanBC))
# 
# Increase = BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
# increaseDF <- data.frame(IncreaseBC = c(0,(Increase)), rank = factor(c(1:(length(Increase) + 1))))
# BC_ranked <- left_join(BC_ranked, increaseDF)
# BC_ranked <- BC_ranked[-nrow(BC_ranked),]
# 
# rm(BCaddition, BCfull, x, x_names, i, temp_BC, temp_BC_matrix, start_matrix, increaseDF, otu_PA, df_s, df_full, df_a, add_matrix, otu_add, otu_start, Increase)
# 
# BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
# 
# elbow <- which.max(BC_ranked$fo_diffs)
# 
# (elbow.all.p <- ggplot(BC_ranked[1:250,], aes(x = factor(BC_ranked$rank[1:250], levels = BC_ranked$rank[1:250]))) +
#   geom_point(aes(y = proportionBC)) +
#   theme_classic() +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(size = 7, angle = 45)) +
#   geom_vline(xintercept = elbow, lty = 3, col = 'red', cex = .5) +
#   geom_vline(xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),
#              lty = 3, col = 'blue', cex = .5) +
#   labs(x = 'ranked OTUs', y = 'Bray-Curtis similarity') +
#   annotate(geom = "text",
#            x = elbow + 10,
#            y = .15,
#            label = paste("Elbow method"," (",elbow,")", sep = ''), color = "red") +
#   annotate(geom = "text",
#            x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])) - 4,
#            y = .08,
#            label = paste("Last 2% increase (",last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),')', sep = ''), color = "blue"))
# 
# occ_abun$fill <- 'no'
# occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))]] <- 'core'
# 
# 
# saveRDS(occ_abun, "Data/ITS/Core/ITS_occ_abun_fam_all.RDS")
# saveRDS(BC_ranked, "Data/ITS/Core/ITS_BC_ranked_fam_all.RDS")

occ_abun.fam.all <- readRDS("Data/ITS/Core/ITS_occ_abun_fam_all.RDS")
BC_ranked.fam.all <- readRDS("Data/ITS/Core/ITS_BC_ranked_fam_all.RDS")

taxon <- as.matrix(rownames(otu))
rownames(taxon) <- taxon
spp <- t(otu)

obs.np <- sncm.fit(spp, taxon, stats = FALSE, pool = NULL) # this doesnt work because of NAs in teh data

(core.all.p <- ggplot() +
  geom_point(data = occ_abun.fam.all[occ_abun.fam.all$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'white', alpha = .2) +
  geom_point(data = occ_abun.fam.all[occ_abun.fam.all$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), 
             pch = 21, fill = 'blue', size = 1.8) +
  geom_line(data = obs.np, aes(y = freq.pred, x = log10(p)), 
            size = 1, color = 'black', alpha = .25) +
  geom_line(data = obs.np, aes(y = pred.upr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25)+
  geom_line(data = obs.np, aes(y = obs.np$pred.lwr, x = log10(p)), 
            color = 'black', lty = 'twodash', size = 1, alpha = .25) +
  labs(title = "16S Core: All Samples", x = "log10(mean relative abundance)", y = "Occupancy"))

core.all <- occ_abun.fam.all$otu[occ_abun.fam.all$fill == 'core']

otu_relabun.fam.all <- decostand(otu, method = "total", MARGIN = 2)

plotDF <-  data.frame(otu = as.factor(row.names(otu_relabun.fam.all)), otu_relabun.fam.all) %>% 
  gather(SampleID.occ, relabun, -otu) %>%
  left_join(map, by = 'SampleID.occ') %>%
  left_join(otu_ranked, by = 'otu') %>%
  filter(otu %in% core.all) %>% 
  group_by(FunGroup, otu) %>%
  dplyr::summarise(plot_freq = sum(relabun>0)/length(relabun), 
            coreSite = ifelse(plot_freq == 1, 1, 0), 
            detect = ifelse(plot_freq > 0, 1, 0))

plotDF$otu <- factor(plotDF$otu, levels = otu_ranked$otu[1:78]) 

(all.plotDF <- ggplot(plotDF, aes(x = otu, plot_freq, group = FunGroup, fill = FunGroup)) +    
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(plotDF$otu %in% core.all))) +
  theme(axis.text = element_text(size = 6)) +
  labs(x = 'Ranked OTUs', y = 'Occupancy by site')) # everything is in like 90% of samples

otu_RA_core.all <- data.frame(otu = as.factor(row.names(otu_relabun.fam.all)), otu_relabun.fam.all) %>%
  filter(otu %in% plotDF$otu)

## divide core into neutral vs deterministic vs dispersal limited ##
otu.mod.all <- data.frame(otu = occ_abun.fam.all[occ_abun.fam.all$fill == "core",]$otu, mod = rep(NA, length(unique(occ_abun.fam.all[occ_abun.fam.all$fill == "core",]$otu))))

for(i in otu.mod.all$otu){
otu.mod.all[otu.mod.all$otu == i,]$mod <- ifelse(log10(occ_abun.fam.all[occ_abun.fam.all$otu == i,]$otu_occ) > log10(obs.np[obs.np$V1 == i,]$pred.upr), "det", "neu")

otu.mod.all[otu.mod.all$otu == i,]$mod <- ifelse(log10(occ_abun.fam.all[occ_abun.fam.all$otu == i,]$otu_occ) < log10(obs.np[obs.np$V1 == i,]$pred.lwr), "dis", otu.mod.all[otu.mod.all$otu == i,]$mod)

}
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

sample_data(ps.16s.nocontrols.nobss)$FunGroup <- relevel(as.factor(sample_data(ps.16s.nocontrols.nobss)$FunGroup), "Grass")

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

plot.list$F.G + plot.list$FG.G 

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

write.csv(df.res, "Data/16s.ddseq.fungroup.csv")


#same as above but now with "forb" as the comparison
sample_data(ps.16s.nocontrols.nobss)$FunGroup <- relevel(as.factor(sample_data(ps.16s.nocontrols.nobss)$FunGroup), "Forb")

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

# try to extract normalized data for use in models
dds.16s.treat = DESeq(treat.16s, test = "Wald", fitType = "parametric")
vsd.blind <- varianceStabilizingTransformation(dds.16s.treat , blind = TRUE)
norm.blind <- assay(vsd.blind)
norm.count <- counts(dds.16s.treat, normalized = TRUE)

norm.count <- data.frame(cbind(SampleID_Fix = rownames(norm.count), norm.count))

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



#### Differential abundance analysis: 16s - Family ####
## 16S: Soil v Rhiz  ##
ps.16s.nocontrols.fam <- tax_glom(ps.16s.nocontrols, taxrank = "Family", NArm = FALSE )

treat.16s = phyloseq_to_deseq2(ps.16s.nocontrols.fam, ~ SampleSubType)

dds.16s.treat = DESeq(treat.16s, test = "Wald", fitType = "parametric")

resultsNames(dds.16s.treat) #gives us comparisons for our contrasts, if needed

alpha = 0.01

#get results 
res <- results(dds.16s.treat, pAdjustMethod = "bonferroni")

#filter results by p-value
res.alpha <- res[which(res$padj < alpha), ]

#Bind taxonomy to results
res.alpha.tax = cbind(as(res.alpha, "data.frame"), as(tax_table(ps.16s.nocontrols.fam)[rownames(res.alpha), ], "matrix"))

#tidy results 
res <- tidy(res.alpha)

#generate plot of significant ASVs for each contrast
# Order order
x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Order, function(x) max(x))
x = sort(x, TRUE)
res.alpha.tax$Order = factor(as.character(res.alpha.tax$Order), levels=names(x))
# Family order
x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Family, function(x) max(x))
x = sort(x, TRUE)
res.alpha.tax$Family = factor(as.character(res.alpha.tax$Family), levels=names(x))
p <- ggplot(res.alpha.tax, aes(x=Family, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("Plant Rhizosphere vs. Background soil")

#plot results

p

#tidy results into a table to save
df.res <- as.data.frame(res)
names(df.res)[1] <- "ASV"

df.res.tax <- cbind(df.res, as(tax_table(ps.16s.nocontrols.fam)[df.res$ASV, ], "matrix"))

write.csv(df.res.tax, "Data/16s.ddseq.soilrhiz.fam.csv")


## 16S:  G v SA v ST  ##
ps.16s.nocontrols.nobss <- subset_samples(ps.16s.nocontrols.fam, SampleSubType !="Background_soil_sand_mix")

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

alpha = 0.05
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
  # Family order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Family, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Family = factor(as.character(res.alpha.tax$Family), levels=names(x))
  p <- ggplot(res.alpha.tax, aes(x=Family, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(plot.name.list[[i]])
  plot.list[[i]] = p
}

#plot results

plot.list$SA.G + plot.list$ST.G + plot.list$SAG.G + plot.list$STG.G

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

df.res.tax <- cbind(df.res, as(tax_table(ps.16s.nocontrols.fam)[df.res$ASV, ], "matrix"))

write.csv(df.res, "Data/16s.ddseq.treatment.fam.csv")


## 16S: G v F ##
#make "grass" the result to compare to

sample_data(ps.16s.nocontrols.nobss)$FunGroup <- relevel(as.factor(sample_data(ps.16s.nocontrols.nobss)$FunGroup), "Grass")

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
  # Family order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Family, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Family = factor(as.character(res.alpha.tax$Family), levels=names(x))
  p <- ggplot(res.alpha.tax, aes(x=Family, y=log2FoldChange, color=Family)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(plot.name.list[[i]])
  plot.list[[i]] = p
}

#plot results

plot.list$F.G+ plot.list$FG.G 

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

df.res.tax <- cbind(df.res, as(tax_table(ps.16s.nocontrols.fam)[df.res$ASV, ], "matrix"))

write.csv(df.res, "Data/16s.ddseq.fungroup.fam.csv")


#same as above but now with "forb" as the comparison
ps.16s.nocontrols.nobss <- subset_samples(ps.16s.nocontrols.fam, SampleSubType !="Background_soil_sand_mix")

sample_data(ps.16s.nocontrols.nobss)$FunGroup <- relevel(as.factor(sample_data(ps.16s.nocontrols.nobss)$FunGroup), "Forb")

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
  # Family order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Family, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Family = factor(as.character(res.alpha.tax$Family), levels=names(x))
  p <- ggplot(res.alpha.tax, aes(x=Family, y=log2FoldChange, color=Family)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(plot.name.list[[i]])
  plot.list[[i]] = p
}


#plot results

plot.list$G.F+ plot.list$FG.F 

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

df.res.tax <- cbind(df.res, as(tax_table(ps.16s.nocontrols.fam)[df.res$ASV, ], "matrix"))

write.csv(df.res, "Data/16s.ddseq.fungroup.v2.fam.csv")

# More families are high abundance in grass x forb vs. forbs alone (vs. grass x forb vs. grasses alone)
# Does this suggest that forbs shift more than grasses during competion?
# Or that grasses are stronger drivers of rhizosphere communities?


## 16S: Competition ##
treat.16s = phyloseq_to_deseq2(ps.16s.nocontrols.nobss, ~ Competion)

dds.16s.treat = DESeq(treat.16s, test = "Wald", fitType = "parametric")

resultsNames(dds.16s.treat) #gives us comparisons for our contrasts


alpha = 0.01

#get results 
res <- results(dds.16s.treat, pAdjustMethod = "bonferroni")

#filter results by p-value
res.alpha <- res[which(res$padj < alpha), ]

#Bind taxonomy to results
res.alpha.tax = cbind(as(res.alpha, "data.frame"), as(tax_table(ps.16s.nocontrols.nobss)[rownames(res.alpha), ], "matrix"))

#tidy results 
res <- tidy(res.alpha)

#generate plot of significant ASVs for each contrast
# Order order
x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Order, function(x) max(x))
x = sort(x, TRUE)
res.alpha.tax$Order = factor(as.character(res.alpha.tax$Order), levels=names(x))
# Family order
x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Family, function(x) max(x))
x = sort(x, TRUE)
res.alpha.tax$Family = factor(as.character(res.alpha.tax$Family), levels=names(x))
p <- ggplot(res.alpha.tax, aes(x=Family, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("Two species vs. One species")

#plot results

p

#tidy results into a table to save
df.res <- as.data.frame(res)
names(df.res)[1] <- "ASV"

df.res.tax <- cbind(df.res, as(tax_table(ps.16s.nocontrols.fam)[df.res$ASV, ], "matrix"))

write.csv(df.res.tax, "Data/16s.ddseq.competition.fam.csv")

# Some families differ during competition (where as ASVs did not)
# however the log2fold change is very low, so may not be informative








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


## ITS:  G v SA v ST  ##
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


## ITS: G v F ##
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


## ITS: Competition ##
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





#### Differential abundance analysis: ITS - Family ####
ps.its.nocontrols.fam <- tax_glom(ps.ITS.nocontrols, taxrank = "Family", NArm = FALSE )


## ITS: Soil v Rhiz  ##

treat.its = phyloseq_to_deseq2(ps.its.nocontrols.fam, ~ SampleSubType)

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
# Family order
x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Family, function(x) max(x))
x = sort(x, TRUE)
res.alpha.tax$Family = factor(as.character(res.alpha.tax$Family), levels=names(x))
p <- ggplot(res.alpha.tax, aes(x=Family, y=log2FoldChange, color=Family)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("Plant Rhizosphere vs. Background soil")

#plot results

p

#tidy results into a table to save
df.res <- as.data.frame(res)
names(df.res)[1] <- "ASV"

df.res.tax <- cbind(df.res, as(tax_table(ps.its.nocontrols.fam)[df.res$ASV, ], "matrix"))
write.csv(df.res, "Data/its.ddseq.soilrhiz.fam.csv")


## ITS:  G v SA v ST  ##
ps.its.nocontrols.nobss <- subset_samples(ps.its.nocontrols.fam, SampleSubType !="Background_soil_sand_mix")

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

alpha = 0.1
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
  # Family order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Family, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Family = factor(as.character(res.alpha.tax$Family), levels=names(x))
  p <- ggplot(res.alpha.tax, aes(x=Family, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(plot.name.list[[i]])
  plot.list[[i]] = p
}

#plot results

plot.list$SA.G + plot.list$ST.G + plot.list$SAG.G + plot.list$STG.G

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

df.res.tax <- cbind(df.res, as(tax_table(ps.its.nocontrols.fam)[df.res$ASV, ], "matrix"))

write.csv(df.res, "Data/its.ddseq.treatment.fam.csv")

#results disappear at family level, means ASV level differences are important here


## ITS: G v F ##
#make "grass" the result to compare to

sample_data(ps.its.nocontrols.nobss)$FunGroup <- relevel(as.factor(sample_data(ps.its.nocontrols.nobss)$FunGroup), "Grass")

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
  # Family order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Family, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Family = factor(as.character(res.alpha.tax$Family), levels=names(x))
  p <- ggplot(res.alpha.tax, aes(x=Family, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(plot.name.list[[i]])
  plot.list[[i]] = p
}

#plot results

plot.list$F.G+ plot.list$FG.G 

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

df.res.tax <- cbind(df.res, as(tax_table(ps.its.nocontrols.nobss)[df.res$ASV, ], "matrix"))

write.csv(df.res, "Data/its.ddseq.fungroup.fam.csv")


#same as above but now with "forb" as the comparison
sample_data(ps.its.nocontrols.nobss)$FunGroup <- relevel(as.factor(sample_data(ps.its.nocontrols.nobss)$FunGroup), "Forb")

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
  # Family order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Family, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Family = factor(as.character(res.alpha.tax$Family), levels=names(x))
  p <- ggplot(res.alpha.tax, aes(x=Family, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(plot.name.list[[i]])
  plot.list[[i]] = p
}

#plot results

plot.list$G.F+ plot.list$FG.F 

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

df.res.tax <- cbind(df.res, as(tax_table(ps.its.nocontrols.nobss)[df.res$ASV, ], "matrix"))

write.csv(df.res, "Data/its.ddseq.fungroup.v2.fam.csv")




## ITS: Competition ##
treat.its = phyloseq_to_deseq2(ps.its.nocontrols.nobss, ~ Competion)

dds.its.treat = DESeq(treat.its, test = "Wald", fitType = "parametric")

resultsNames(dds.its.treat) #gives us comparisons for our contrasts


alpha = 0.01

#get results 
res <- results(dds.its.treat, pAdjustMethod = "bonferroni")

#filter results by p-value
res.alpha <- res[which(res$padj < alpha), ]


#Bind taxonomy to results
res.alpha.tax = cbind(as(res.alpha, "data.frame"), as(tax_table(ps.its.nocontrols.nobss)[rownames(res.alpha), ], "matrix"))

#tidy results 
res <- tidy(res.alpha)

#generate plot of significant ASVs for each contrast
# Order order
x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Order, function(x) max(x))
x = sort(x, TRUE)
res.alpha.tax$Order = factor(as.character(res.alpha.tax$Order), levels=names(x))
# Family order
x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Family, function(x) max(x))
x = sort(x, TRUE)
res.alpha.tax$Family = factor(as.character(res.alpha.tax$Family), levels=names(x))
p <- ggplot(res.alpha.tax, aes(x=Family, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("Two species vs. One species")

#plot results

p

#tidy results into a table to save
df.res <- as.data.frame(res)
names(df.res)[1] <- "ASV"

df.res.tax <- cbind(df.res, as(tax_table(ps.its.nocontrols.nobss)[df.res$ASV, ], "matrix"))

write.csv(df.res.tax, "Data/its.ddseq.competition.fam.csv")

# Some families differ during competition (where as ASVs did not)
# however the log2fold change is very low, so may not be informative



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
#results <- predict(st,asvs_16s[test.ix,], alpha1=alpha1, alpha2=alpha2)

results <- readRDS("Data/source_tracker/16s.treatment.ST.RDS")

# Estimate leave-one-out source proportions in training data 
# results.train <- predict(st, alpha1=alpha1, alpha2=alpha2)

results.train <- readRDS("Data/source_tracker/16s.treatment.ST.RDS")



# plot results
labels <- sprintf('%s %s', envs,desc)
plot(results, labels[test.ix], type='pie')

#Restructure results for plotting
results.df <- as_tibble(results$proportions)
results.df$SampleID_Fix <- row.names(results$proportions)

results.df.meta <- inner_join(results.df, metadata_16s)

grouped_res <- group_by(results.df.meta, TreatmentName)
#calculate mean of means + sd (sd(x)/sqrt(length(x)))
avgs_res <- summarise(grouped_res, mean=100*mean(Background_soil), sd=100*sd(Background_soil)/sqrt(length(Background_soil)))

#Plot % community from background soil (rest of community is thus 'unknown')
ggplot(avgs_res, aes(x = TreatmentName, y = mean, fill = TreatmentName)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
                width = .4, position = position_dodge(.9)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 20)) + 
  ylab("Mean % Community from Background Soil")

kable(avgs_res, caption = "% of community from background soil")

# other plotting functions
# plot(results, labels[test.ix], type='bar')
#plot(results, labels[test.ix], type='dist')
#plot(results.train, labels[train.ix], type='pie')
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
#results.fg <- predict(st.fg,asvs_16s[test.ix.fg,], alpha1=alpha1, alpha2=alpha2)

results.fg <- readRDS("Data/source_tracker/16s.fg.ST.RDS")

# Estimate leave-one-out source proportions in training data 
#results.train.fg <- predict(st.fg, alpha1=alpha1, alpha2=alpha2)

results.train.fg <- readRDS("Data/source_tracker/16s.fg.ST.loo.RDS")
  
# plot results
labels.fg <- sprintf('%s %s', envs.fg,desc)
plot(results.fg, labels.fg[test.ix.fg], type='pie')


#Restructure results for plotting
results.df.fg <- as_tibble(results.fg$proportions)
results.df.fg$SampleID_Fix <- row.names(results.fg$proportions)

results.df.fg.meta <- inner_join(results.df.fg, metadata_16s)

grouped_res.fg <- group_by(results.df.fg.meta, FunGroup)
#calculate mean of means + sd (sd(x)/sqrt(length(x)))
avgs_res.fg <- summarise(grouped_res.fg, 
                         mean_grass=100*mean(Grass), 
                         sd_grass=100*sd(Grass)/sqrt(length(Grass)),
                         mean_forb=100*mean(Forb), 
                         sd_forb=100*sd(Forb)/sqrt(length(Forb)),
                         mean_soil=100*mean(Background_soil), 
                         sd_soil=100*sd(Background_soil)/sqrt(length(Background_soil)),
                         mean_unknown=100*mean(Unknown), 
                         sd_unknown=100*sd(Unknown)/sqrt(length(Unknown)))

avgs_res.fg.df <- data.frame("Source" = c('grass', 'forb', 'background soil', 'unknown'), "mean" = c(40.3419,42.13405,8.56619,8.957857), "sd" = c(1.325134,1.142347, 0.4818445,0.5978808))

kable(avgs_res.fg.df, caption = "% of community from source")
# 
# |Source          |      mean|        sd|
#   |:---------------|---------:|---------:|
#   |grass           | 40.341900| 1.3251340|
#   |forb            | 42.134050| 1.1423470|
#   |background soil |  8.566190| 0.4818445|
#   |unknown         |  8.957857| 0.5978808|

ggplot(avgs_res.fg.df, aes(x = Source, y = mean, fill = Source)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
                width = .4, position = position_dodge(.9)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 20)) + 
  ylab("Mean % Community from Source during Competition")

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
#results.its <- predict(st.its,asvs_its[test.ix.its,], alpha1=alpha1, alpha2=alpha2)

results.its <- readRDS("Data/source_tracker/its.treatment.ST.RDS")

# Estimate leave-one-out source proportions in training data 
#results.train.its <- predict(st.its, alpha1=alpha1, alpha2=alpha2)

results.train.its <- readRDS("Data/source_tracker/its.treatment.ST.loo.RDS")


# plot results
labels.its <- sprintf('%s %s', envs.its,desc)
plot(results.its, labels.its[test.ix.its], type='pie')


#Restructure results for plotting
results.df.its <- as_tibble(results.its$proportions)
results.df.its$SampleID_Fix <- row.names(results.its$proportions)

results.df.its.meta <- inner_join(results.df.its, metadata_its)

grouped_res.its <- group_by(results.df.its.meta, TreatmentName)
#calculate mean of means + sd (sd(x)/sqrt(length(x)))
avgs_res.its <- summarise(grouped_res.its, mean=100*mean(Background_soil), sd=100*sd(Background_soil)/sqrt(length(Background_soil)))

#Plot % community from background soil (rest of community is thus 'unknown')
ggplot(avgs_res.its, aes(x = TreatmentName, y = mean, fill = TreatmentName)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
                width = .4, position = position_dodge(.9)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 20)) + 
  ylab("Mean % Community from Background Soil")

kable(avgs_res.its, caption = "% of community from background soil")


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
#results.fg.its <- predict(st.fg.its,asvs_its[test.ix.fg.its,], alpha1=alpha1, alpha2=alpha2)

results.fg.its <- readRDS("Data/source_tracker/its.fg.ST.RDS")

# Estimate leave-one-out source proportions in training data 
#results.train.fg.its <- predict(st.fg.its, alpha1=alpha1, alpha2=alpha2)

results.train.fg.its <- readRDS("Data/source_tracker/its.fg.ST.loo.RDS")

# plot results
labels.fg.its <- sprintf('%s %s', envs.fg.its,desc)
plot(results.fg.its, labels.fg.its[test.ix.fg.its], type='pie')



#Restructure results for plotting
results.df.fg.its <- as_tibble(results.fg.its$proportions)
results.df.fg.its$SampleID_Fix <- row.names(results.fg.its$proportions)

results.df.fg.meta.its <- inner_join(results.df.fg.its, metadata_its)

grouped_res.fg.its <- group_by(results.df.fg.meta.its, FunGroup)
#calculate mean of means + sd (sd(x)/sqrt(length(x)))
avgs_res.fg.its <- summarise(grouped_res.fg.its, 
                         mean_grass=100*mean(Grass), 
                         sd_grass=100*sd(Grass)/sqrt(length(Grass)),
                         mean_forb=100*mean(Forb), 
                         sd_forb=100*sd(Forb)/sqrt(length(Forb)),
                         mean_soil=100*mean(Background_soil), 
                         sd_soil=100*sd(Background_soil)/sqrt(length(Background_soil)),
                         mean_unknown=100*mean(Unknown), 
                         sd_unknown=100*sd(Unknown)/sqrt(length(Unknown)))

avgs_res.fg.its.df <- data.frame("Source" = c('grass', 'forb', 'background soil', 'unknown'), "mean" = c(38.87845,32.64786,5.160238,23.31345), "sd" = c(2.263107,2.005945, 0.6747386,2.107941))

kable(avgs_res.fg.its.df, caption = "% of community from source")
# 
# |Source          |      mean|        sd|
#   |:---------------|---------:|---------:|
#   |grass           | 38.878450| 2.2631070|
#   |forb            | 32.647860| 2.0059450|
#   |background soil |  5.160238| 0.6747386|
#   |unknown         | 23.313450| 2.1079410|

ggplot(avgs_res.fg.its.df, aes(x = Source, y = mean, fill = Source)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
                width = .4, position = position_dodge(.9)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -70, hjust = 0, vjust = .5),
        text = element_text(size = 20)) + 
  ylab("Mean % Community from Source during Competition")


### Stats on sourcetracker  ###

## 16S ##

#The null hypothesis of these tests is that “sample distribution is normal”. 
#If the test is significant, the distribution is non-normal.
shapiro.test(grouped_res.fg$Grass) #normal
shapiro.test(grouped_res.fg$Forb) #normal
shapiro.test(grouped_res.fg$Background_soil) #non-normal
shapiro.test(grouped_res.fg$Unknown) #non-normal

#Probably should use a KW Test? But would need to reformat data so groupings work
#In mean time....

#Do grasses have the same mean as Forbs? 
wilcox.test(grouped_res.fg$Grass, grouped_res.fg$Forb) # p-value = 0.3413
wilcox.test(grouped_res.fg$Unknown, grouped_res.fg$Background_soil) #p-value = 0.6604
#rest are significant (can run if neeed)

## ITS ##

#The null hypothesis of these tests is that “sample distribution is normal”. 
#If the test is significant, the distribution is non-normal.
shapiro.test(grouped_res.fg.its$Grass) #non-normal
shapiro.test(grouped_res.fg.its$Forb) #non-normal
shapiro.test(grouped_res.fg.its$Background_soil) #non-normal
shapiro.test(grouped_res.fg.its$Unknown) #non-normal

#Probably should use a KW Test? But would need to reformat data so groupings work
#In mean time....

#Do grasses have the same mean as Forbs? 
wilcox.test(grouped_res.fg.its$Grass, grouped_res.fg.its$Forb) # p-value = 0.06625
#rest are significant (can run if neeed)

## Comparions ##

#ITS vs. 16S comparisons
wilcox.test(grouped_res.fg$Background_soil, grouped_res.fg.its$Background_soil) #p-value = 2.192e-09
wilcox.test(grouped_res.fg$Unknown, grouped_res.fg.its$Unknown) #p-value = 5.575e-09
wilcox.test(grouped_res.fg$Forb, grouped_res.fg.its$Forb) #p-value = 0.0001254
wilcox.test(grouped_res.fg$Grass, grouped_res.fg.its$Grass) #p-value = 0.3317






#### Biomass effects ####
hist(log(mapping$Weight.g))

mapping$FunGroup <- as.factor(mapping$FunGroup)

# Forb weight
m.bio.f <- lmer(log(forb.wt) ~ FunGroup + (1|forb), mapping[mapping$FunGroup != "Grass",])
plot(fitted(m.bio.f), resid(m.bio.f))
qqnorm(resid(m.bio.f))
qqline(resid(m.bio.f))
summary(m.bio.f)

# Grass weight
m.bio.g <- lmer(log(grass.wt) ~ FunGroup + (1|grass), mapping[mapping$FunGroup != "Forb",])
plot(fitted(m.bio.g), resid(m.bio.g))
qqnorm(resid(m.bio.g))
qqline(resid(m.bio.g))
summary(m.bio.g)

# Log response ratio
map.comp <- mapping %>%
  filter(Competion == "TwoSpecies") %>%
  pivot_longer(c(grass, forb), names_to = "FunGroup2", values_to = "Species")

biomass <- mapping %>%
  filter(Competion == "SingleSpecies") %>%
  dplyr::group_by(PlantSpeciesSampled) %>%
  dplyr::summarise(avg.wt = mean(Weight.g, na.rm = T)) %>%
  dplyr::rename(Species = PlantSpeciesSampled) %>%
  left_join(map.comp, by = "Species") %>%
  filter(Competion == "TwoSpecies")

biomass$weight.d <- ifelse(biomass$FunGroup2 == "forb", log(biomass$forb.wt/biomass$avg.wt), log(biomass$grass.wt/biomass$avg.wt))


m.bio <- lmer(weight.d ~ FunGroup2 + (1|Treatment), biomass)
plot(fitted(m.bio), resid(m.bio))
qqnorm(resid(m.bio))
qqline(resid(m.bio))
summary(m.bio)

library(Rmisc)
biomass.lrr.sum <- summarySE(biomass, groupvars = "FunGroup2", measurevar = "weight.d", na.rm = T)
biomass.lrr.sum$FunGroup2 <- recode(biomass.lrr.sum$FunGroup2, grass = "Grass", forb = "Forb")

bio.plot <- ggplot(biomass.lrr.sum, aes(y = weight.d, x = FunGroup2, fill = FunGroup2)) + 
  geom_bar(stat = "identity", position = "dodge", col = "black") + 
  geom_errorbar(aes(ymin = weight.d - se, ymax = weight.d + se), width = 0.2, position = position_dodge(1)) + 
  theme_bw(base_size = 15) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  labs(y = "Log Response Ratio") +
  scale_fill_manual(values = c("magenta4", "cyan4"))

ggsave("Figures/bio-lrr-plot.jpg", bio.plot, width = 5, height = 4, units = "in", dpi = 600)

biomass$Species <- factor(biomass$Species, levels = c("Clarkia purpurea", "Agoseris heterophylla", "Calycadenia pauciflora", "Hemizonia congesta", "Lasthenia californica", "Plantago erecta", "Elymus caput-medusae", "Bromus hordeaceus", "Avena fatua"))  

# biomass$Species.sh <- revalue(biomass$Species, c("Bromus hordeaceus" = "BRHO", "Clarkia purpurea" =  "CLPU", "Agoseris heterophylla" = "AGHE", "AVFA" = "Avena fatua", "CAPA" = "Calycadenia pauciflora", "HECO" = "Hemizonia congesta", "LACA" = "Lasthenia californica", "PLER" = "Plantago erecta", "TACA" = "Taeniatherum caput-medusae"))


# ggplot(biomass, aes(y = weight.d, x = Species, col = FunGroup2)) +
#   geom_boxplot()
# 
# ggplot(biomass, aes(y = weight.d, x = FunGroup2)) +
#   geom_point()

### looking at target and neighbor ##
biomass <- separate(biomass, PlantSpeciesSampled, c("species1", "species2"), sep = ",")
biomass <- biomass[!is.na(biomass$weight.d),]
biomass$neighbor <- ifelse(biomass$FunGroup2 == "forb", biomass$species1, biomass$species2)
biomass$neighbor2 <- ifelse(biomass$FunGroup2 == "forb", "grass", "forb")

bio.nt <- lmer(weight.d ~ neighbor2 + (1|Treatment), biomass) 
plot(fitted(bio.nt), resid(bio.nt))
qqnorm(resid(bio.nt))
qqline(resid(bio.nt))
summary(bio.nt)# competitive effects dont vary between grasses and forbs...

bio.nt <- lmer(weight.d ~ FunGroup2 + (1|Treatment), biomass) 
plot(fitted(bio.nt), resid(bio.nt))
qqnorm(resid(bio.nt))
qqline(resid(bio.nt))
summary(bio.nt) # competitive response doesnt vary between grasses and forbs

biomass.f <- filter(mapping[,c(1, 6, 8, 9, 13, 28, 30, 35)], forb != "none", !is.na(forb.wt))
biomass.g <- filter(mapping[,c(1, 6, 8, 9, 13, 27, 31, 35)], grass != "none")
biomass.f$g.f <- "Forb"
biomass.g$g.f <- "Grass"
colnames(biomass.g)[6:7] <- c("Species", "weight.g")
colnames(biomass.f)[6:7] <- c("Species", "weight.g")
biomass2 <- rbind(biomass.g, biomass.f)

biomass2$FunGroup <- ifelse(biomass2$FunGroup == "grass_x_forb", paste(biomass2$FunGroup, biomass2$g.f, sep = "."), as.character(biomass2$FunGroup))

hist(log(biomass2$weight.g))
m.bio2 <- lmer(log(weight.g) ~ FunGroup - 1 + (1|Species), biomass2)
plot(fitted(m.bio2), resid(m.bio2))
qqnorm(resid(m.bio2))
qqline(resid(m.bio2))
summary(m.bio2)

pairs(emmeans(m.bio2, ~FunGroup))

# contrasts

Fb <- c(1,0,0,0)
G <- c(0,1,0,0)
fxgF <- c(0,0,1,0)
fxgG <- c(0,0,0,1)

K <- rbind("Fb - G" = Fb - G,
           "Fb - fxgF" = Fb - fxgF,
           "G - fxgG" = G - fxgG,
           "fxgF  - fxgG" = fxgF  - fxgG)

library(multcomp)
summary(glht(m.bio2, linfct = K), test = adjusted("BH")) # forbs declined more than grasses in competition

hist(log(biomass2$weight.g))
m.bio2 <- lmer(log(weight.g) ~ Competion * g.f + (1|Species), biomass2)
plot(fitted(m.bio2), resid(m.bio2))
qqnorm(resid(m.bio2))
qqline(resid(m.bio2))
summary(m.bio2)

library(Rmisc)

biomass.sum <- summarySE(biomass2, measurevar = "weight.g", groupvars = c("Competion","g.f", "Species"))

biomass.sum$Competion <- as.numeric(biomass.sum$Competion)

ggplot(biomass.sum, aes(y = weight.g, x = Competion, col = Species, group = Species)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = weight.g - se, ymax = weight.g + se), width = 0.1) +
  facet_wrap(~g.f) +
  theme_classic()

biomass.sum <- summarySE(biomass2, measurevar = "weight.g", groupvars = c("Competion","g.f"))

ggplot(biomass.sum, aes(y = weight.g, x = Competion)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = weight.g - se, ymax = weight.g + se), width = 0.2, position = position_dodge(1)) +
  facet_wrap(~g.f) +
  theme_classic()

#### Biomass v Diversity ####
GL_Alpha <- estimate_richness(ps.16s.nocontrols.rare, measures = c("Observed","Shannon", "InvSimpson"))
GL_Alpha2 <- cbind(GL_Alpha, sample_data(ps.16s.nocontrols.rare))

GL_Alpha2 <- filter(GL_Alpha2, !is.na(Competion))

hist(log(GL_Alpha2$Shannon + 1))
m.sha <- lmer(log(Shannon + 1) ~ FunGroup + (1|Treatment), data = GL_Alpha2[GL_Alpha2$Shannon > 5.7,])
plot(fitted(m.sha), resid(m.sha))
qqnorm(resid(m.sha))
qqline(resid(m.sha))
summary(m.sha)


m.b <- lmer(log(forb.wt) ~ Shannon + FunGroup + (1 | forb), data = GL_Alpha2[GL_Alpha2$Shannon > 5.7 & GL_Alpha2$forb != "none",])
plot(fitted(m.b), resid(m.b))
qqnorm(resid(m.b))
qqline(resid(m.b))
summary(m.b)

GL_Alpha2$shannon.log <- log(GL_Alpha2$Shannon)
GL_Alpha2$forb.wt.log <- log(GL_Alpha2$forb.wt)

psem.16S <- psem(
  
 lme(shannon.log ~ FunGroup, random = ~ 1 | Treatment, data = GL_Alpha2[GL_Alpha2$forb != "none",], na.action = na.omit),
 
 lme(forb.wt.log ~ shannon.log + FunGroup, random = ~ 1 | Treatment, data = GL_Alpha2[GL_Alpha2$forb != "none",], na.action = na.omit)
  
  )

summary(psem.16S, .progressBar = FALSE) 


# reformatting 
# GL_Alpha2 <- GL_Alpha2[,c(2, 4, 9, 16, 12, 30, 31, 33, 34, 35, 38)]
# 
# GL_Alpha2.spp.f <- filter(GL_Alpha2[,-c(6,8,9)], Competion == "SingleSpecies", FunGroup == "Forb")
# names(GL_Alpha2.spp.f)[6] <- "Species_Name"
# 
# GL_Alpha2.spp.g <- filter(GL_Alpha2[,-c(7,8,9)], Competion == "SingleSpecies", FunGroup == "Grass")
# names(GL_Alpha2.spp.g)[6] <- "Species_Name"
# 
# GL_Alpha2.com.f <- filter(GL_Alpha2[,-c(6, 9, 10)], Competion == "TwoSpecies")
# names(GL_Alpha2.com.f)[6] <- "Species_Name"
# names(GL_Alpha2.com.f)[7] <- "Weight.g"
# 
# GL_Alpha2.com.g <- filter(GL_Alpha2[,-c(7, 8, 10)], Competion == "TwoSpecies")
# names(GL_Alpha2.com.g)[6] <- "Species_Name"
# names(GL_Alpha2.com.g)[7] <- "Weight.g"
# 
# GL_Alpha2 <- rbind(GL_Alpha2.spp.f, GL_Alpha2.spp.g, GL_Alpha2.com.f, GL_Alpha2.com.g)
# 
# GL_Alpha2$shannon.log <- log(GL_Alpha2$Shannon)
# GL_Alpha2$weight.log <- log(GL_Alpha2$Weight.g)
# 
# GL_Alpha2 <- filter(GL_Alpha2, !is.na(Weight.g))
# 
# psem.16S <- psem(
#   
#  lme(shannon.log ~ FunGroup, random = ~ 1 | Species_Name, data = GL_Alpha2),
#  
#  lme(weight.log ~ shannon.log + FunGroup, random = ~ 1 | Species_Name, data = GL_Alpha2)
#   
#   )
# 
# summary(psem.16S, .progressBar = FALSE) 
# 
# m.wt <- lme(weight.log ~ Shannon * FunGroup, random = ~ 1 | Treatment, data = GL_Alpha2)
# plot(fitted(m.wt), resid(m.wt))
# qqnorm(resid(m.wt))
# qqline(resid(m.wt))
# summary(m.wt)

#### Biomass v Ordination ####
ord.ax <- data.frame(cbind(SampleID_Fix = rownames(ps.rare.ord.tr$vectors), ax1 = ps.rare.ord.tr$vectors[,1]))
ord.ax$ax1 <- as.numeric(as.character(ord.ax$ax1))
ord.ax <- merge(ord.ax, mapping, by = "SampleID_Fix")

hist(ord.ax$ax1)
m.b <- lmer(log(forb.wt) ~ FunGroup + ax1 + (1|forb), data = ord.ax[ord.ax$FunGroup != "Grass",])
plot(fitted(m.b), resid(m.b))
qqnorm(resid(m.b))
qqline(resid(m.b))
summary(m.b)

m.b <- lmer(log(grass.wt) ~ FunGroup + ax1 + (1|grass), data = ord.ax[ord.ax$FunGroup != "Forb",])
plot(fitted(m.b), resid(m.b))
qqnorm(resid(m.b))
qqline(resid(m.b))
summary(m.b)

psem.16S <- psem(
  
 lmer(ax1 ~ FunGroup + (1 | forb), data = ord.ax[ord.ax$FunGroup != "Grass",], na.action = na.omit),
 
 lmer(forb.wt.log ~ ax1 + FunGroup + (1 | forb), data = df_ASV[df_ASV$forb != "none" & df_ASV$OTU == "SV645",], na.action = na.omit)
  
  )

summary(psem.16S, .progressBar = FALSE) 

#### Biomass v Order ####
ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion)))
ps.16s.nocontrols.rare.RA <- transform_sample_counts(ps.16s.nocontrols.rare, function(x) x / sum(x))
ps.16s.nocontrols.RA.ord <- tax_glom(ps.16s.nocontrols.rare.RA, taxrank = "Order", NArm = FALSE )
SEM.ord.df <- psmelt(ps.16s.nocontrols.RA.ord)

SEM.ord.df2 <- SEM.ord.df %>%
  filter(Competion == "SingleSpecies", forb != "none") %>%
  group_by(forb, Order) %>%
  summarise(avg.wt = mean(Weight.g, na.rm = T), avg.RA = mean(Abundance)) %>%
  left_join(SEM.ord.df, by = c("forb", "Order")) %>%
  filter(Competion == "TwoSpecies") %>%
  mutate(weight.d = log(forb.wt/avg.wt), Order.d = (Abundance - avg.RA)/avg.RA)

###
# Selenomonadales
###

# abundance v weight
ggplot(SEM.ord.df[SEM.ord.df$Order == "Selenomonadales",], aes(x = log(Abundance + 0.0001), y = log(forb.wt + 0.01), col = FunGroup, group = FunGroup)) +
  geom_point() +
  geom_smooth(method = "lm")

m.b <- lmer(log(forb.wt + 0.01) ~ log(Abundance + 0.0001) + FunGroup + (1|forb), data = SEM.ord.df[SEM.ord.df$Order == "Selenomonadales",])
plot(fitted(m.b), resid(m.b))
qqnorm(resid(m.b))
qqline(resid(m.b))
summary(m.b)

### 
# Flavobacteriales
###

ggplot(SEM.ord.df[SEM.ord.df$Order == "Flavobacteriales" & SEM.ord.df$FunGroup != "Grass",], aes(x = Abundance, y = log(forb.wt + 0.01), col = FunGroup, group = FunGroup)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(SEM.ord.df2[SEM.ord.df2$Order == "Flavobacteriales",], aes(x = Abundance, y = weight.d)) + 
  geom_point() +
  geom_smooth(method = "lm")

m.b <- lmer(log(forb.wt + 0.01) ~ Abundance + FunGroup + (1|forb), data = SEM.ord.df[SEM.ord.df$Order == "Flavobacteriales",])
plot(fitted(m.b), resid(m.b))
qqnorm(resid(m.b))
qqline(resid(m.b))
summary(m.b)

SEM.ord.df <- filter(SEM.ord.df, forb != "none", !is.na(forb.wt))
SEM.ord.df$Abun.log <- log(SEM.ord.df$Abundance + 0.0001)
SEM.ord.df$forb.wt.log <- log(SEM.ord.df$forb.wt + 0.01)

psem.16S <- psem(
 
 lme(forb.wt.log ~ Abundance + FunGroup, random =  ~ 1 | forb, na.action = na.omit, SEM.ord.df[SEM.ord.df$Order == "Flavobacteriales",]),
 
 lme(Abundance ~ FunGroup, random =  ~ 1 | forb, na.action = na.omit, data = SEM.ord.df[SEM.ord.df$Order == "Flavobacteriales",])
  
  )

summary(psem.16S, .progressBar = FALSE) 

###
# Betaproteobacteriales
###
ggplot(SEM.ord.df[SEM.ord.df$Order == "Betaproteobacteriales" & SEM.ord.df$FunGroup != "Grass",], aes(x = Abundance, y = log(forb.wt + 0.01), col = FunGroup, group = FunGroup)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(SEM.ord.df2[SEM.ord.df2$Order == "Betaproteobacteriales",], aes(x = Abundance, y = weight.d)) + 
  geom_point() +
  geom_smooth(method = "lm")

m.b <- lmer(log(forb.wt + 0.01) ~ Abundance + FunGroup + (1|forb), data = SEM.ord.df[SEM.ord.df$Order == "Betaproteobacteriales",])
plot(fitted(m.b), resid(m.b))
qqnorm(resid(m.b))
qqline(resid(m.b))
summary(m.b)

psem.16S <- psem(
 
 lme(forb.wt.log ~ Abundance + FunGroup, random =  ~ 1 | forb, na.action = na.omit, SEM.ord.df[SEM.ord.df$Order == "Betaproteobacteriales",]),
 
 lme(Abundance ~ FunGroup, random =  ~ 1 | forb, na.action = na.omit, data = SEM.ord.df[SEM.ord.df$Order == "Betaproteobacteriales",])
  
  )

summary(psem.16S, .progressBar = FALSE) 

###
# Clostridiales
###
ggplot(SEM.ord.df[SEM.ord.df$Order == "Clostridiales" & SEM.ord.df$FunGroup != "Grass",], aes(x = Abundance, y = log(forb.wt + 0.01), col = FunGroup, group = FunGroup)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(SEM.ord.df2[SEM.ord.df2$Order == "Clostridiales",], aes(x = Abundance, y = weight.d)) + 
  geom_point() +
  geom_smooth(method = "lm")

m.b <- lmer(log(forb.wt + 0.01) ~ Abundance * FunGroup + (1|forb), data = SEM.ord.df[SEM.ord.df$Order == "Clostridiales",])
plot(fitted(m.b), resid(m.b))
qqnorm(resid(m.b))
qqline(resid(m.b))
summary(m.b)

psem.16S <- psem(
 
 lme(forb.wt.log ~ Abundance + FunGroup, random =  ~ 1 | forb, na.action = na.omit, SEM.ord.df[SEM.ord.df$Order == "Clostridiales",]),
 
 lme(Abundance ~ FunGroup, random =  ~ 1 | forb, na.action = na.omit, data = SEM.ord.df[SEM.ord.df$Order == "Clostridiales",])
  
  )

summary(psem.16S, .progressBar = FALSE) 

#### Biomass v Family ####

ps.16s.nocontrols.RA.fam <- tax_glom(ps.16s.nocontrols.rare.RA, taxrank = "Family", NArm = FALSE )
SEM.fam.df <- psmelt(ps.16s.nocontrols.RA.fam)


fam.df.f <- filter(SEM.fam.df[,-c(30,34,35)], forb != "none", !is.na(forb.wt))
fam.df.g <- filter(SEM.fam.df[,-c(31,33,35)], grass != "none")
fam.df.f$g.f <- "Forb"
fam.df.g$g.f <- "Grass"
colnames(fam.df.g)[c(30,32)] <- c("Species", "weight.g")
colnames(fam.df.f)[c(30,32)] <- c("Species", "weight.g")
fam.df <- rbind(fam.df.g, fam.df.f)

fam.df$FunGroup <- ifelse(fam.df$FunGroup == "grass_x_forb", paste(fam.df$FunGroup, fam.df$g.f, sep = "."), as.character(fam.df$FunGroup))

fam.df$Competion <- revalue(fam.df$Competion, c("SingleSpecies" = "Alone", "TwoSpecies" = "Together"))

# SEM.fam.df2 <- SEM.fam.df %>%
#   filter(Competion == "SingleSpecies") %>%
#   group_by(forb, Family) %>%
#   summarise(avg.wt = mean(Weight.g, na.rm = T), avg.RA = mean(Abundance)) %>%
#   left_join(SEM.fam.df, by = c("forb", "Family")) %>%
#   filter(Competion == "TwoSpecies") %>%
#   mutate(weight.d = log(forb.wt/avg.wt), Family.d = (Abundance - avg.RA)/avg.RA)

SEM.fam.comp <- SEM.fam.df %>%
  filter(Competion == "TwoSpecies") %>%
  pivot_longer(c(grass, forb), names_to = "FunGroup2", values_to = "Species")

SEM.fam.df3 <- SEM.fam.df %>%
  dplyr::filter(Competion == "SingleSpecies") %>%
  dplyr::group_by(PlantSpeciesSampled, Family) %>%
  dplyr::summarise(avg.wt = mean(Weight.g, na.rm = T), avg.RA = mean(Abundance)) %>%
  dplyr::rename(Species = PlantSpeciesSampled) %>%
  left_join(SEM.fam.comp, by = c("Species", "Family")) %>%
  filter(Competion == "TwoSpecies") 
# %>%
#   mutate(weight.d = log(forb.wt/avg.wt), Family.d = (Abundance - avg.RA)/avg.RA)
  
SEM.fam.df3$weight.d <- ifelse(SEM.fam.df3$FunGroup2 == "forb", log(SEM.fam.df3$forb.wt/SEM.fam.df3$avg.wt), log(SEM.fam.df3$grass.wt/SEM.fam.df3$avg.wt))

SEM.fam.df3 <- separate(SEM.fam.df3, PlantSpeciesSampled, c("species1", "species2"), sep = ",")
SEM.fam.df3$neighbor <- ifelse(SEM.fam.df3$FunGroup2 == "forb", SEM.fam.df3$species1, SEM.fam.df3$species2)
SEM.fam.df3$neighbor2 <- ifelse(SEM.fam.df3$FunGroup2 == "forb", "grass", "forb")

#### Fibrobacteraceae ####

## Abundance vs weight
ggplot(fam.df[fam.df$Family == "Fibrobacteraceae" & fam.df$g.f == "Forb",], aes(x = Abundance, y = log(weight.g + 0.01), col = Species, group = Species)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Fibrobacteraceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance + FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Fibrobacteraceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Fibrobacteraceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Fibrobacteraceae", y = "Log Response Ratio")

m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Fibrobacteraceae",])
plot(fitted(m.f.lrr), resid(m.f.lrr))
qqnorm(resid(m.f.lrr))
qqline(resid(m.f.lrr))
summary(m.f.lrr)

m.f.lrr <- lmer(weight.d ~ Abundance * neighbor2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Fibrobacteraceae",])
plot(fitted(m.f.lrr), resid(m.f.lrr))
qqnorm(resid(m.f.lrr))
qqline(resid(m.f.lrr))
summary(m.f.lrr)

#### Veillonellaceae ####

## Abundance vs weight
ggplot(fam.df[fam.df$Family == "Veillonellaceae",], aes(x = Abundance, y = log(weight.g + 0.01), col = g.f, group = g.f)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Veillonellaceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance * FunGroup2 + (1|Treatment), data = fam.df[fam.df$Family == "Veillonellaceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Veillonellaceae" & SEM.fam.df3$Abundance < 0.01,], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Veillonellaceae", y = "Log Response Ratio")

m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Veillonellaceae" & SEM.fam.df3$Abundance < 0.01,])
plot(fitted(m.f.lrr), resid(m.f.lrr))
qqnorm(resid(m.f.lrr))
qqline(resid(m.f.lrr))
summary(m.f.lrr)

#### Clostridiaceae_1 ####

## Abundance vs weight
ggplot(fam.df[fam.df$Family == "Clostridiaceae_1",], aes(x = Abundance, y = log(weight.g + 0.01), col = g.f, group = g.f)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Clostridiaceae_1", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance * FunGroup2 + (1|Treatment), data = fam.df[fam.df$Family == "Clostridiaceae_1",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Clostridiaceae_1" & SEM.fam.df3$Abundance < 0.01,], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Clostridiaceae_1", y = "Log Response Ratio")

m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Clostridiaceae_1",])
plot(fitted(m.f.lrr), resid(m.f.lrr))
qqnorm(resid(m.f.lrr))
qqline(resid(m.f.lrr))
summary(m.f.lrr)

# outliers removed
m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Clostridiaceae_1" & SEM.fam.df3$Abundance < 0.01,])
plot(fitted(m.f.lrr), resid(m.f.lrr))
qqnorm(resid(m.f.lrr))
qqline(resid(m.f.lrr))
summary(m.f.lrr)

#### Methylophilaceae ####


## Abundance vs weight
ggplot(fam.df[fam.df$Family == "Methylophilaceae",], aes(x = Abundance, y = log(weight.g + 0.01), col = g.f, group = g.f)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Methylophilaceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance + FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Methylophilaceae" & fam.df$g.f != "Forb",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Methylophilaceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Methylophilaceae", y = "Log Response Ratio")

m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Methylophilaceae",])
plot(fitted(m.f.lrr), resid(m.f.lrr))
qqnorm(resid(m.f.lrr))
qqline(resid(m.f.lrr))
summary(m.f.lrr)

#### **- Burkholderiaceae ####

## Abundance vs weight
ggplot(fam.df[fam.df$Family == "Burkholderiaceae",], aes(x = Abundance, y = log(weight.g + 0.01), col = g.f, group = g.f)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Burkholderiaceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance * FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Burkholderiaceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
# ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Burkholderiaceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x) +
#   theme_bw() +
#   theme(
#     strip.background = element_rect(fill = "white"),
#     legend.title = element_blank(),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   labs(title = "Burkholderiaceae", y = "Log Response Ratio")
# 
# m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Burkholderiaceae",])
# plot(fitted(m.f.lrr), resid(m.f.lrr))
# qqnorm(resid(m.f.lrr))
# qqline(resid(m.f.lrr))
# summary(m.f.lrr)

#### Cyclobacteriaceae ####

## Abundance vs weight
ggplot(fam.df[fam.df$Family == "Cyclobacteriaceae",], aes(x = Abundance, y = log(weight.g + 0.01), col = g.f, group = g.f)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Cyclobacteriaceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance + FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Cyclobacteriaceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Cyclobacteriaceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Cyclobacteriaceae", y = "Log Response Ratio")

m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Cyclobacteriaceae",])
plot(fitted(m.f.lrr), resid(m.f.lrr))
qqnorm(resid(m.f.lrr))
qqline(resid(m.f.lrr))
summary(m.f.lrr)



#### Nitrospiraceae ####


## Abundance vs weight
ggplot(fam.df[fam.df$Family == "Nitrospiraceae",], aes(x = Abundance, y = log(weight.g + 0.01), col = g.f, group = g.f)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Nitrospiraceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance * FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Nitrospiraceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
# ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Nitrospiraceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x) +
#   theme_bw() +
#   theme(
#     strip.background = element_rect(fill = "white"),
#     legend.title = element_blank(),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   labs(title = "Nitrospiraceae", y = "Log Response Ratio")
# 
# m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Nitrospiraceae",])
# plot(fitted(m.f.lrr), resid(m.f.lrr))
# qqnorm(resid(m.f.lrr))
# qqline(resid(m.f.lrr))
# summary(m.f.lrr)

#### Steroidobacteraceae ####


## Abundance vs weight
ggplot(fam.df[fam.df$Family == "Steroidobacteraceae",], aes(x = Abundance, y = log(weight.g + 0.01), col = g.f, group = g.f)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Steroidobacteraceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance * FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Steroidobacteraceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Steroidobacteraceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Steroidobacteraceae", y = "Log Response Ratio")

m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Steroidobacteraceae",])
plot(fitted(m.f.lrr), resid(m.f.lrr))
qqnorm(resid(m.f.lrr))
qqline(resid(m.f.lrr))
summary(m.f.lrr)

#### Rhodocyclaceae ####

## Abundance vs weight
ggplot(fam.df[fam.df$Family == "Rhodocyclaceae",], aes(x = Abundance, y = log(weight.g + 0.01), col = g.f, group = g.f)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Rhodocyclaceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance + FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Rhodocyclaceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Rhodocyclaceae" & SEM.fam.df3$Abundance < 0.04,], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Rhodocyclaceae", y = "Log Response Ratio")

m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Rhodocyclaceae",])
plot(fitted(m.f.lrr), resid(m.f.lrr))
qqnorm(resid(m.f.lrr))
qqline(resid(m.f.lrr))
summary(m.f.lrr)

#### Weeksellaceae ####

## Abundance vs weight
ggplot(fam.df[fam.df$Family == "Weeksellaceae",], aes(x = Abundance, y = log(weight.g + 0.01), col = g.f, group = g.f)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Weeksellaceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance + FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Weeksellaceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Weeksellaceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Weeksellaceae", y = "Log Response Ratio")

m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Weeksellaceae",])
plot(fitted(m.f.lrr), resid(m.f.lrr))
qqnorm(resid(m.f.lrr))
qqline(resid(m.f.lrr))
summary(m.f.lrr)

#### Streptomycetaceae ####

## Abundance vs weight
ggplot(fam.df[fam.df$Family == "Streptomycetaceae" & fam.df$Abundance < 0.01,], aes(x = Abundance, y = log(weight.g + 0.01), col = g.f, group = g.f)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Streptomycetaceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance + FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Streptomycetaceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Streptomycetaceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Streptomycetaceae", y = "Log Response Ratio")

m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Streptomycetaceae",])
plot(fitted(m.f.lrr), resid(m.f.lrr))
qqnorm(resid(m.f.lrr))
qqline(resid(m.f.lrr))
summary(m.f.lrr)

#### Microbacteriaceae ####

## Abundance vs weight
ggplot(fam.df[fam.df$Family == "Microbacteriaceae",], aes(x = Abundance, y = log(weight.g + 0.01), col = g.f, group = g.f)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Microbacteriaceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance + FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Microbacteriaceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
# ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Microbacteriaceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x) +
#   theme_bw() +
#   theme(
#     strip.background = element_rect(fill = "white"),
#     legend.title = element_blank(),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   labs(title = "Microbacteriaceae", y = "Log Response Ratio")
# 
# m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Microbacteriaceae",])
# plot(fitted(m.f.lrr), resid(m.f.lrr))
# qqnorm(resid(m.f.lrr))
# qqline(resid(m.f.lrr))
# summary(m.f.lrr)

#### **+ Beijerinckiaceae ####

## Abundance vs weight ** driven by Agoseris
ggplot(fam.df[fam.df$Family == "Beijerinckiaceae",], aes(x = Abundance, y = log(weight.g + 0.01), col = g.f, group = g.f)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Beijerinckiaceae", y = "Log Weight (g)")

m.f.wt <- lm(log(weight.g + 0.05) ~ Abundance * FunGroup, data = fam.df[fam.df$Family == "Beijerinckiaceae" & fam.df$g.f == "Forb",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
# ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Beijerinckiaceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x) +
#   theme_bw() +
#   theme(
#     strip.background = element_rect(fill = "white"),
#     legend.title = element_blank(),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   labs(title = "Beijerinckiaceae", y = "Log Response Ratio")
# 
# m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Beijerinckiaceae",])
# plot(fitted(m.f.lrr), resid(m.f.lrr))
# qqnorm(resid(m.f.lrr))
# qqline(resid(m.f.lrr))
# summary(m.f.lrr)

#### **+ Acetobacteraceae ####

## Abundance vs weight ** driven by plantago
ggplot(fam.df[fam.df$Family == "Acetobacteraceae",], aes(x = Abundance, y = log(weight.g + 0.01), col = Species, group = Species)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Acetobacteraceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.03) ~ Abundance * FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Acetobacteraceae" & fam.df$g.f == "Forb",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
# ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Acetobacteraceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x) +
#   theme_bw() +
#   theme(
#     strip.background = element_rect(fill = "white"),
#     legend.title = element_blank(),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   labs(title = "Acetobacteraceae", y = "Log Response Ratio")
# 
# m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Acetobacteraceae",])
# plot(fitted(m.f.lrr), resid(m.f.lrr))
# qqnorm(resid(m.f.lrr))
# qqline(resid(m.f.lrr))
# summary(m.f.lrr)

#### **+ Azospirillaceae ####

## Abundance vs weight ** driven by plantago
ggplot(fam.df[fam.df$Family == "Azospirillaceae",], aes(x = Abundance, y = log(weight.g + 0.01), col = Species, group = Species)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Azospirillaceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance * FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Azospirillaceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
# ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Azospirillaceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x) +
#   theme_bw() +
#   theme(
#     strip.background = element_rect(fill = "white"),
#     legend.title = element_blank(),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   labs(title = "Azospirillaceae", y = "Log Response Ratio")
# 
# m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Azospirillaceae",])
# plot(fitted(m.f.lrr), resid(m.f.lrr))
# qqnorm(resid(m.f.lrr))
# qqline(resid(m.f.lrr))
# summary(m.f.lrr)

#### **+ Geodermatophilaceae ####

## Abundance vs weight ** driven by plantago
ggplot(fam.df[fam.df$Family == "Geodermatophilaceae" & fam.df$g.f == "Forb",], aes(x = Abundance, y = log(weight.g + 0.01), col = Species, group = Species)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Geodermatophilaceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance * FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Geodermatophilaceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
# ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Geodermatophilaceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x) +
#   theme_bw() +
#   theme(
#     strip.background = element_rect(fill = "white"),
#     legend.title = element_blank(),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   labs(title = "Geodermatophilaceae", y = "Log Response Ratio")
# 
# m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Geodermatophilaceae",])
# plot(fitted(m.f.lrr), resid(m.f.lrr))
# qqnorm(resid(m.f.lrr))
# qqline(resid(m.f.lrr))
# summary(m.f.lrr)

#### **+ Propionibacteriaceae ####

## Abundance vs weight ## driven by plantago mostly, some hemizonia
ggplot(fam.df[fam.df$Family == "Propionibacteriaceae" & fam.df$g.f == "Forb",], aes(x = Abundance, y = log(weight.g + 0.01), col = Species, group = Species)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Propionibacteriaceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance * FunGroup + (1|Treatment), data = fam.df[fam.df$Family == "Propionibacteriaceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
# ggplot(SEM.fam.df3[SEM.fam.df3$Family == "Propionibacteriaceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x) +
#   theme_bw() +
#   theme(
#     strip.background = element_rect(fill = "white"),
#     legend.title = element_blank(),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   labs(title = "Propionibacteriaceae", y = "Log Response Ratio")
# 
# m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Propionibacteriaceae",])
# plot(fitted(m.f.lrr), resid(m.f.lrr))
# qqnorm(resid(m.f.lrr))
# qqline(resid(m.f.lrr))
# summary(m.f.lrr)

#### Pairs ####
# fam <- res.alpha.tax$Family
# fam.df <- fam.df[fam.df$Family %in% fam,]
# group <- unique(fam.df$FunGroup)
# 
# fam.df.w <- pivot_wider(fam.df[,c(3,32,35)], names_from = FunGroup, values_from = weight.g)
# 
# plot_list = list()
# model_list = list()
# 
# for(i in fam){
#   for(j in group) {
#     p <- ggplot(fam.df[fam.df$Family == i & fam.df$FunGroup == j,], aes(y = log(weight.g + 0.01), x = Abundance)) +
#       geom_point() +
#       geom_smooth(method = "lm", se = F, formula = y ~ x, col = "red") +
#       labs(title = paste(i,j)) +
#       theme_bw()
#       
#    plot_list[[i]][j] <- p 
#    
#    model_list[[i]][j] <- summary(lmer(log(weight.g + 0.01) ~ Abundance + (1|Species), fam.df[fam.df$Family == i & fam.df$FunGroup == j,]))
#   }
# }
# 
# print(plot_list[[1]][1])

#### ITS vs biomass ####
ps.ITS.nocontrols.rare.RA <- transform_sample_counts(ps.ITS.nocontrols.rare, function(x) x / sum(x))
ps.ITS.nocontrols.RA.fam <- tax_glom(ps.ITS.nocontrols.rare.RA, taxrank = "Family", NArm = FALSE )
ITS.fam <- psmelt(ps.ITS.nocontrols.RA.fam)

ITS.fam.df.f <- filter(ITS.fam[,-c(30,34,35)], forb != "none", !is.na(forb.wt))
ITS.fam.df.g <- filter(ITS.fam[,-c(31,33,35)], grass != "none")
ITS.fam.df.f$g.f <- "Forb"
ITS.fam.df.g$g.f <- "Grass"
colnames(ITS.fam.df.g)[c(30,32)] <- c("Species", "weight.g")
colnames(ITS.fam.df.f)[c(30,32)] <- c("Species", "weight.g")
ITS.fam.df <- rbind(ITS.fam.df.g, ITS.fam.df.f)

ITS.fam.df$FunGroup <- ifelse(ITS.fam.df$FunGroup == "grass_x_forb", paste(ITS.fam.df$FunGroup, ITS.fam.df$g.f, sep = "."), as.character(ITS.fam.df$FunGroup))

ITS.fam.df$Competion <- revalue(ITS.fam.df$Competion, c("SingleSpecies" = "Alone", "TwoSpecies" = "Together"))

ITS.fam.comp <- ITS.fam %>%
  dplyr::filter(Competion == "TwoSpecies") %>%
  pivot_longer(c(grass, forb), names_to = "FunGroup2", values_to = "Species")


ITS.fam.df3 <- ITS.fam %>%
  dplyr::filter(Competion == "SingleSpecies") %>%
  dplyr::group_by(PlantSpeciesSampled, Family) %>%
  dplyr::summarise(avg.wt = mean(Weight.g, na.rm = T), avg.RA = mean(Abundance)) %>%
  dplyr::rename(Species = PlantSpeciesSampled) %>%
  left_join(ITS.fam.comp, by = c("Species", "Family")) %>%
  filter(Competion == "TwoSpecies") 
# %>%
#   mutate(weight.d = log(forb.wt/avg.wt), Family.d = (Abundance - avg.RA)/avg.RA)
  
ITS.fam.df3$weight.d <- ifelse(ITS.fam.df3$FunGroup2 == "forb", log(ITS.fam.df3$forb.wt/ITS.fam.df3$avg.wt), log(ITS.fam.df3$grass.wt/ITS.fam.df3$avg.wt))

ITS.fam.df3 <- separate(ITS.fam.df3, PlantSpeciesSampled, c("species1", "species2"), sep = ",")
ITS.fam.df3$neighbor <- ifelse(ITS.fam.df3$FunGroup2 == "forb", ITS.fam.df3$species1, ITS.fam.df3$species2)
ITS.fam.df3$neighbor2 <- ifelse(ITS.fam.df3$FunGroup2 == "forb", "grass", "forb")

# Tubeufiaceae #

## Abundance vs weight
ggplot(ITS.fam.df[ITS.fam.df$Family == "f__Tubeufiaceae",], aes(x = Abundance, y = log(weight.g + 0.01), col = g.f, group = g.f)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~ Competion) + 
  labs(title = "Tubeufiaceae", y = "Log Weight (g)")

m.f.wt <- lmer(log(weight.g + 0.01) ~ Abundance * FunGroup + (1|Treatment), data = ITS.fam.df[ITS.fam.df$Family == "f__Tubeufiaceae",])
plot(fitted(m.f.wt), resid(m.f.wt))
qqnorm(resid(m.f.wt))
qqline(resid(m.f.wt))
summary(m.f.wt)

## log response ratio 
ggplot(ITS.fam.df3[ITS.fam.df3$Family == "f__Tubeufiaceae",], aes(x = Abundance, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Tubeufiaceae", y = "Log Response Ratio")

# m.f.lrr <- lmer(weight.d ~ Abundance * FunGroup2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Tubeufiaceae",])
# plot(fitted(m.f.lrr), resid(m.f.lrr))
# qqnorm(resid(m.f.lrr))
# qqline(resid(m.f.lrr))
# summary(m.f.lrr)
# 
# m.f.lrr <- lmer(weight.d ~ Abundance * neighbor2 + (1|Treatment), data = SEM.fam.df3[SEM.fam.df3$Family == "Tubeufiaceae",])
# plot(fitted(m.f.lrr), resid(m.f.lrr))
# qqnorm(resid(m.f.lrr))
# qqline(resid(m.f.lrr))
# summary(m.f.lrr)
