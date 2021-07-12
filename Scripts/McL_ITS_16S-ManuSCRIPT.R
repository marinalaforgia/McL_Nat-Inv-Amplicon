### 16S and ITS Analysis August 2020 ###
rm(list=ls())

#### Set Dir ####

#setwd("/Users/Marina/Documents/UC-Davis/Projects/McL_Nat-Inv-Amplicon") # Marina's Working Directory
#setwd("/Users/Cassie/Documents/Dropbox/Greenhouse_ITS/McL_Nat-Inv-Amplicon") # Cassie's Working Directory

#### Load Libraries ####

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
library(gridExtra)
library(gtable)
library(sjstats)
library(MuMIn)

#### Load Data ####

### 
# Unrarefied
###

## 16S ##
ps.16s.nocontrols <- readRDS("Data/16S/Intermediates/ps-nocontrols-unrare.RDS") # 16s data, Marina's path
ps.16s.nocontrols <- readRDS("Data/ps-nocontrols-unrare.RDS") #16s data, Cassie's path
#132 samples

## ITS ##
ps.ITS.nocontrols <- readRDS("Data/ITS/greenhouse_its_itsx_decontam_controlsremoved.rds") # ITS data
ps.ITS.nocontrols <- readRDS("Data/greenhouse_its_itsx_decontam_controlsremoved.rds") # ITS data, Cassie's path
#132 samples

###
# Rarefied Data 
###

## 16S ##

# Load 
ps.16s.nocontrols.rare <- readRDS("Data/16S/Intermediates/ps-nocontrols-rare9434.RDS") # rarefied 16s data
ps.16s.nocontrols.rare <- readRDS("Data/ps-nocontrols-rare9434.RDS") # rarefied 16s data, Cassie's path

sample_data(ps.16s.nocontrols.rare)$PlantSpeciesSampled <- recode(sample_data(ps.16s.nocontrols.rare)$PlantSpeciesSampled, `Elymus caput-medusae` = "Taeniatherum caput-medusae")

# Load 
ps.ITS.nocontrols.rare <- readRDS("Data/ITS/greenhouse_its_itsx_decontam_controlsremoved_rare7557.rds") # rarefied ITS data
ps.ITS.nocontrols.rare <- readRDS("Data/greenhouse_its_itsx_decontam_controlsremoved_rare7557.rds") # rarefied ITS data, Cassie's path

sample_data(ps.ITS.nocontrols.rare)$PlantSpeciesSampled <- recode(sample_data(ps.ITS.nocontrols.rare)$PlantSpeciesSampled, `Elymus caput-medusae` = "Taeniatherum caput-medusae")

###
# Metadata mapping file including biomass and traits (version4)
###

mapping <- read.csv("/Users/Marina/Documents/UC-Davis/Projects/McL_Nat-Inv-Amplicon/Data/Setup/Grassland-Amplicon-Mapping-File4.csv")
mapping <- read.csv("Data/Grassland-Amplicon-Mapping-File4.csv") #Cassie's path

###
# Single species and RA files
###

## 16s ##

ps.16s.nocontrols.rare.spp <- subset_samples(ps.16s.nocontrols.rare, Competion == "SingleSpecies")

ps.16s.nocontrols.rare.spp.RA <- transform_sample_counts(ps.16s.nocontrols.rare.spp, function(x) x / sum(x)) # relative abundance (used later)

## ITS ##

ps.ITS.nocontrols.rare.spp <- subset_samples(ps.ITS.nocontrols.rare, Competion == "SingleSpecies") 

ps.ITS.nocontrols.rare.spp.RA <- transform_sample_counts(ps.ITS.nocontrols.rare.spp, function(x) x / sum(x)) 

#### DESeq prep ####

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

#### Ordination (16S) - Soil ####
ps.16s.nocontrols.rare.soil <- subset_samples(ps.16s.nocontrols.rare, SampleType == "Soil")

sample_data(ps.16s.nocontrols.rare.soil)$SampleSubType <- recode(sample_data(ps.16s.nocontrols.rare.soil)$SampleSubType, Rhizosphere_soil = "Rhizosphere", Background_soil_sand_mix = "Background soil")
  
ps.rare.ord.tr.soil <- ordinate(ps.16s.nocontrols.rare.soil, "PCoA", "wunifrac")

DistWU <- phyloseq::distance(ps.16s.nocontrols.rare.soil, method = "wunifrac", type = "samples")

# adonis compares centroids and is sensitive to differences in spread
set.seed(50)
adonis(DistWU ~ SampleSubType, as(sample_data(ps.16s.nocontrols.rare.soil), "data.frame"), permutations = 9999)

C <- betadisper(DistWU, as(sample_data(ps.16s.nocontrols.rare.soil), "data.frame")$SampleSubType)
set.seed(50)
permutest(C, permutations = 9999) # spread varies between treatments
plot(C, label = F) # centroids


#### Ordination (16S) - FG ####
ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion)))

sample_data(ps.16s.nocontrols.rare)$FunGroup <- recode(sample_data(ps.16s.nocontrols.rare)$FunGroup, grass_x_forb = "Competition")

sample_data(ps.16s.nocontrols.rare)$FunGroup <- factor(sample_data(ps.16s.nocontrols.rare)$FunGroup, levels = c("Forb", "Grass", "Competition"))
  
ps.rare.ord.tr <- ordinate(ps.16s.nocontrols.rare, "PCoA", "wunifrac")

DistWU <- phyloseq::distance(ps.16s.nocontrols.rare, method = "wunifrac", type = "samples")

# adonis compares centroids and is sensitive to differences in spread
set.seed(50)
adonis(DistWU ~ FunGroup, as(sample_data(ps.16s.nocontrols.rare), "data.frame"), permutations = 9999)

C <- betadisper(DistWU, as(sample_data(ps.16s.nocontrols.rare), "data.frame")$FunGroup)
set.seed(50)
permutest(C, permutations = 9999) # spread doesnt vary between treatments
#plot(C, label = F) # centroids

set.seed(50)
pair.16s.com <- adonis.pair(DistWU, as(sample_data(ps.16s.nocontrols.rare), "data.frame")$FunGroup, nper = 9999, corr.method = "BH")
kable(pair.16s.com)


#### Ordination (16S) - Spp ####

###
# Forb Species
###

ps.16s.nocontrols.rare.spp.f <- subset_samples(ps.16s.nocontrols.rare.spp, grass == "none")

ps.16s.nc.rare.spp.ord <- ordinate(ps.16s.nocontrols.rare.spp.f, "PCoA", "wunifrac")

ord.plot.spp.f <- plot_ordination(ps.16s.nocontrols.rare.spp.f, ps.16s.nc.rare.spp.ord, color = "PlantSpeciesSampled", shape = "PlantSpeciesSampled") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = PlantSpeciesSampled)) +
  labs(x = "PCoA1 (25.2%)", y = "PCoA2 (19.0%)") +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.margin = margin(c(.1,.1,.1,.1))
  )

ggsave("Figures/Final/Fig-S10a.pdf", ord.plot.spp.f, width = 7, height = 4, units = "in", dpi = 600)


DistWU.Spp <- phyloseq::distance(ps.16s.nocontrols.rare.spp.f, method = "wunifrac", type = "samples")

set.seed(50)
adonis(DistWU.Spp ~ PlantSpeciesSampled, as(sample_data(ps.16s.nocontrols.rare.spp.f), "data.frame"), permutations = 9999)

C <- betadisper(DistWU.Spp, as(sample_data(ps.16s.nocontrols.rare.spp.f), "data.frame")$PlantSpeciesSampled)
set.seed(50)
permutest(C, permutations = 9999) # spread doesnt vary between forbs species

### 
# Grass species
###
ps.16s.nocontrols.rare.spp.g <- subset_samples(ps.16s.nocontrols.rare.spp, forb == "none")

ps.16s.nc.rare.spp.ord <- ordinate(ps.16s.nocontrols.rare.spp.g, "PCoA", "wunifrac")

ord.plot.spp.g <- plot_ordination(ps.16s.nocontrols.rare.spp.g, ps.16s.nc.rare.spp.ord, color = "PlantSpeciesSampled", shape = "PlantSpeciesSampled") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = PlantSpeciesSampled)) +
  labs(x = "PCoA1 (35.2%)", y = "PCoA2 (24.1%)") +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.margin = margin(c(.1,.1,.1,.1))
  )

ggsave("Figures/Final/Fig-S10b.pdf", ord.plot.spp.g, width = 7, height = 4, units = "in", dpi = 600)


DistWU.Spp <- phyloseq::distance(ps.16s.nocontrols.rare.spp.g, method = "wunifrac", type = "samples")

set.seed(50)
adonis(DistWU.Spp ~ PlantSpeciesSampled, as(sample_data(ps.16s.nocontrols.rare.spp.g), "data.frame"), permutations = 9999) # grass microbiomes are sig different, particularly medusahead

C <- betadisper(DistWU.Spp, as(sample_data(ps.16s.nocontrols.rare.spp.g), "data.frame")$PlantSpeciesSampled)
set.seed(50)
permutest(C, permutations = 9999)

set.seed(50)
pair.spp.g <- adonis.pair(DistWU.Spp, as(sample_data(ps.16s.nocontrols.rare.spp.g), "data.frame")$grass, nper = 9999, corr.method = "BH")
kable(pair.spp.g) # TACA sig diff

#### Ordination (ITS) - Soil ####
ps.ITS.nocontrols.rare.soil <- subset_samples(ps.ITS.nocontrols.rare, SampleType == "Soil")

sample_data(ps.ITS.nocontrols.rare.soil)$SampleSubType <- recode(sample_data(ps.ITS.nocontrols.rare.soil)$SampleSubType, Rhizosphere_soil = "Rhizosphere", Background_soil_sand_mix = "Background soil")

GL_pcoa.soil <- ordinate(
  physeq = ps.ITS.nocontrols.rare.soil, 
  method = "PCoA", 
  distance = "bray")

DistBC <- phyloseq::distance(ps.ITS.nocontrols.rare.soil, method = "bray", type="samples")
set.seed(50)
adonis(DistBC ~ SampleSubType, as(sample_data(ps.ITS.nocontrols.rare.soil), "data.frame"), permutations = 9999) 

C <- betadisper(DistBC, as(sample_data(ps.ITS.nocontrols.rare.soil), "data.frame")$SampleSubType)
set.seed(50)
permutest(C, permutations = 9999) 
plot(C, label = F)

#### Ordination (ITS) - FG ####
ps.ITS.nocontrols.rare <- subset_samples(ps.ITS.nocontrols.rare, !(is.na(Competion)))

sample_data(ps.ITS.nocontrols.rare)$FunGroup <- recode(sample_data(ps.ITS.nocontrols.rare)$FunGroup, grass_x_forb = "Competition")

sample_data(ps.ITS.nocontrols.rare)$FunGroup <- factor(sample_data(ps.ITS.nocontrols.rare)$FunGroup, levels = c("Forb", "Grass", "Competition"))

GL_pcoa <- ordinate(
  physeq = ps.ITS.nocontrols.rare, 
  method = "PCoA", 
  distance = "bray")

DistBC <- phyloseq::distance(ps.ITS.nocontrols.rare, method = "bray", type="samples")
set.seed(50)
adonis(DistBC ~ FunGroup, as(sample_data(ps.ITS.nocontrols.rare), "data.frame"), permutations = 9999) # marginal fungroup diffs

# because all these are marginal, we should not be looking at pairwise comparisons
# 
# sample_data(ps.ITS.nocontrols.rare)$FunGroup <- as.factor(sample_data(ps.ITS.nocontrols.rare)$FunGroup)
# set.seed(50)
# pair.its.com <- adonis.pair(DistBC, as(sample_data(ps.ITS.nocontrols.rare), "data.frame")$FunGroup, nper = 9999, corr.method = "BH")
# kable(pair.its.com)

#### Ordination (ITS) - Spp ####

###
# Forb Species
###

ps.ITS.nocontrols.rare.spp.f <- subset_samples(ps.ITS.nocontrols.rare.spp, grass == "none")

ps.ITS.nc.rare.spp.ord <- ordinate(ps.ITS.nocontrols.rare.spp.f, "PCoA", distance = "bray")

ord.plot.spp.f <- plot_ordination(ps.ITS.nocontrols.rare.spp.f, ps.ITS.nc.rare.spp.ord, color = "PlantSpeciesSampled", shape = "PlantSpeciesSampled") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = PlantSpeciesSampled)) +
  labs(x = "PCoA1 (10.7%)", y = "PCoA2 (8.8%)") +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.margin = margin(c(.1,.1,.1,.1))
  )

ggsave("Figures/Final/Fig-S11a.pdf", ord.plot.spp.f, width = 7, height = 4, units = "in", dpi = 600)

DistBC.spp <- phyloseq::distance(ps.ITS.nocontrols.rare.spp.f, method = "bray", type = "samples")

set.seed(50)
adonis(DistBC.spp ~ PlantSpeciesSampled, as(sample_data(ps.ITS.nocontrols.rare.spp.f), "data.frame"), permutations = 9999)

C <- betadisper(DistBC.spp, as(sample_data(ps.ITS.nocontrols.rare.spp.f), "data.frame")$PlantSpeciesSampled)
set.seed(50)
permutest(C, permutations = 9999) # spread doesnt vary between forbs species

### 
# Grass species
###
ps.ITS.nocontrols.rare.spp.g <- subset_samples(ps.ITS.nocontrols.rare.spp, forb == "none")

ps.ITS.nc.rare.spp.ord <- ordinate(ps.ITS.nocontrols.rare.spp.g, "PCoA", "bray")

ord.plot.spp.g <- plot_ordination(ps.ITS.nocontrols.rare.spp.g, ps.ITS.nc.rare.spp.ord, color = "PlantSpeciesSampled", shape = "PlantSpeciesSampled") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = PlantSpeciesSampled)) +
  labs(x = "PCoA1 (18.3%)", y = "PCoA2 (12.7%)") +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.margin = margin(c(.1,.1,.1,.1))
  )

ggsave("Figures/Final/Fig-S11b.pdf", ord.plot.spp.g, width = 7, height = 4, units = "in", dpi = 600)

DistBC.spp <- phyloseq::distance(ps.ITS.nocontrols.rare.spp.g, method = "bray", type = "samples")

set.seed(50)
adonis(DistBC.spp ~ PlantSpeciesSampled, as(sample_data(ps.ITS.nocontrols.rare.spp.g), "data.frame"), permutations = 9999) 

C <- betadisper(DistWU.Spp, as(sample_data(ps.ITS.nocontrols.rare.spp.g), "data.frame")$PlantSpeciesSampled)
set.seed(50)
permutest(C, permutations = 9999)

#### .---Fig: Ord (16S) - Soil ####

ord.plot.16s.s <- plot_ordination(ps.16s.nocontrols.rare.soil, ps.rare.ord.tr.soil, color = "SampleSubType", shape = "SampleSubType") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = SampleSubType)) +
  labs(x = "PCoA1 (26.8%)", y = "PCoA2 (17.6%)") +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.margin = margin(c(.1,.1,.1,.1))
  ) 

ggsave("Figures/Final/Fig-ord-16s-soil.pdf", ord.plot.16s.s, width = 6, height = 4, units = "in", dpi = 600)


#### .---Fig: Ord (16S) - FG ####
ord.plot <- plot_ordination(ps.16s.nocontrols.rare, ps.rare.ord.tr, color = "FunGroup", shape = "FunGroup") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = FunGroup)) +
  labs(x = "PCoA1 (25.6%)", y = "PCoA2 (18.4%)") +
  theme(
    legend.title = element_blank(),
    legend.position = c(.82, .82),
    legend.text = element_text(size = 10),
    legend.margin = margin(c(.1,.1,.1,.1))
  ) +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3"))

ggsave("Figures/Final/Fig-1.pdf", ord.plot, width = 5, height = 4, units = "in", dpi = 600)


#### .---Fig: Ord (ITS) - Soil ####

ord.plot.ITS.s <- plot_ordination(ps.ITS.nocontrols.rare.soil, GL_pcoa.soil, color = "SampleSubType", shape = "SampleSubType") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = SampleSubType)) +
  labs(x = "PCoA1 (10.6%)", y = "PCoA2 (5.5%)") +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.margin = margin(c(.1,.1,.1,.1))
  ) 

ggsave("Figures/Final/Fig-ord-soil.pdf", ord.plot.ITS.s, width = 6, height = 4, units = "in", dpi = 600)


#### .---Fig: Ord (ITS) - FG ####
ord.plot.ITS <- plot_ordination(ps.ITS.nocontrols.rare, GL_pcoa, color = "FunGroup", shape = "FunGroup") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = FunGroup)) +
  labs(x = "PCoA1 (11.2%)", y = "PCoA2 (5.4%)") +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.margin = margin(c(.1,.1,.1,.1))
  ) +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3"))

ggsave("Figures/Final/Fig-S1.pdf", ord.plot.ITS, width = 6, height = 4, units = "in", dpi = 600)

#### Shannon Diveristy (16S) ####

GL_Alpha <- estimate_richness(ps.16s.nocontrols.rare, measures = "Shannon")

GL_Alpha <- cbind(GL_Alpha, sample_data(ps.16s.nocontrols.rare))

set.seed(50)
kruskal_test(Shannon ~ FunGroup, distribution = approximate(nresample = 9999), data = GL_Alpha)

dunnTest(Shannon ~ FunGroup, data = GL_Alpha, method = "bh") 

#### Shannon Diversity (ITS) ####

GL_Alpha_ITS <- estimate_richness(ps.ITS.nocontrols.rare, measures = "Shannon")
GL_Alpha_ITS <- cbind(GL_Alpha_ITS, sample_data(ps.ITS.nocontrols.rare))

set.seed(50)
kruskal_test(Shannon ~ FunGroup, distribution = approximate(nresample = 9999), data = GL_Alpha_ITS)

#### .---Fig: Div (16S) ####

shannon.sum <- ddply(GL_Alpha, .(FunGroup), summarize, max = max(Shannon))
shannon.sum$cld <- c("a", "b", "ab")

p.rich.16s <- plot_richness(ps.16s.nocontrols.rare, measures = "Shannon", x = "FunGroup", color = "FunGroup") + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  theme_bw(base_size = 15) +
  geom_text(data = shannon.sum, aes(x = FunGroup, label = cld, y = max,
        size = 8,
        vjust = -0.7),
    col = "black") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    strip.text = element_blank()
  ) +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(y = "Shannon Diversity") +
  ylim(NA, 6.8)
  
ggsave("Figures/Final/Fig-2.pdf", p.rich.16s, width = 5, height = 4, units = "in", dpi = 600)

#### .---Fig: Div (ITS) ####
p.rich.ITS <- plot_richness(ps.ITS.nocontrols.rare, measures = "Shannon", x = "FunGroup", color = "FunGroup") + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_blank()
  ) +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(y = "Shannon Diversity")

ggsave("Figures/Final/Fig-S2.pdf", p.rich.ITS, width = 7, height = 4, units = "in", dpi = 600)

#### DESeq (16S) ####
ps.16s.nocontrols.fam <- tax_glom(ps.16s.nocontrols, taxrank = "Family", NArm = FALSE )
ps.16s.nocontrols.nobss <- subset_samples(ps.16s.nocontrols.fam, SampleSubType !="Background_soil_sand_mix")

###
# Grass baseline
###

sample_data(ps.16s.nocontrols.nobss)$FunGroup <- relevel(as.factor(sample_data(ps.16s.nocontrols.nobss)$FunGroup), "Grass")

treat.16s <- phyloseq_to_deseq2(ps.16s.nocontrols.nobss, ~ FunGroup)

dds.16s.treat <- DESeq(treat.16s, test = "Wald", fitType = "parametric")

resultsNames(dds.16s.treat) #gives us comparisons for our contrasts

contrasts <- c("F.G", "FG.G")

contrast.list <- list(F.G = c("FunGroup", "Forb", "Grass"),
                      FG.G = c("FunGroup", "grass_x_forb", "Grass"))


plot.name.list <- list(F.G = "Forb vs. Grass",
                       FG.G = "Competition vs. Grass")

alpha = 0.01
res.list <- list()
plot.list <- list()

FG.F <- c("Clostridiaceae_1", "Veillonellaceae", "Fibrobacteraceae", "Methylophilaceae")
FG.G <- c("Burkholderiaceae", "Rhodocyclaceae")
 
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
  res.alpha.tax$Order = factor(as.character(res.alpha.tax$Order), levels = names(x))
  # Family order
  x = tapply(res.alpha.tax$log2FoldChange, res.alpha.tax$Family, function(x) max(x))
  x = sort(x, TRUE)
  res.alpha.tax$Family = factor(as.character(res.alpha.tax$Family), levels = names(x))
   p <- ggplot(res.alpha.tax, aes(x = Family, y = log2FoldChange)) + 
    geom_point(size = 3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 12),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "none",
      axis.title.x = element_blank(),
        ) +
    ylim(-3.5, 3.5) +
    labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change")))
  
   plot.list[[i]] = p
}

#plot results

plot.list$F.G+ plot.list$FG.G 
plot.FG.G <- plot.list$FG.G 

###
# Extract normalized data
###

rld <- rlog(dds.16s.treat, blind = F) 
norm.rld <- assay(rld)
norm.countg <- data.frame(cbind(ASV = rownames(norm.rld), norm.rld))

res.alpha.tax$ASV <- row.names(res.alpha.tax)

norm.countg <- merge(res.alpha.tax, norm.countg, by = "ASV", all.y = F)
norm.countg <- norm.countg[,c(12,15:141)]
norm.countlg <- pivot_longer(norm.countg, col = -Family, names_to = "Sample", values_to = "norm.counts")

norm.countlg$Sample <- gsub("\\.", "-", norm.countlg$Sample)
norm.countlg$norm.counts <- as.numeric(as.character(norm.countlg$norm.counts))

###
# Forb baseline
###

ps.16s.nocontrols.nobss <- subset_samples(ps.16s.nocontrols.fam, SampleSubType !="Background_soil_sand_mix")

sample_data(ps.16s.nocontrols.nobss)$FunGroup <- relevel(as.factor(sample_data(ps.16s.nocontrols.nobss)$FunGroup), "Forb")

treat.16s <- phyloseq_to_deseq2(ps.16s.nocontrols.nobss, ~ FunGroup)

dds.16s.treat <- DESeq(treat.16s, test = "Wald", fitType = "parametric")

resultsNames(dds.16s.treat) #gives us comparisons for our contrasts

contrasts <- c("G.F", "FG.F")

contrast.list <- list(G.F = c("FunGroup", "Grass", "Forb"),
                      FG.F = c("FunGroup", "grass_x_forb", "Forb"))


plot.name.list <- list(G.F = "Grass vs. Forb",
                       FG.F = "Competition vs. Forb")

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

  
   p <- ggplot(res.alpha.tax[res.alpha.tax$Family %in% fam.new,], aes(x = Family, y = log2FoldChange)) + 
    geom_point(size = 3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 12),
      legend.text = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.key.size = unit(.35, "cm"),
      legend.title = element_text(size = 10),
      legend.position = "none",
      axis.title.x = element_blank()
     ) +
     ylim(-3.5,3.5) +
    labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change")))
  
  plot.list[[i]] = p
}


#plot results

d.plot <- plot.list$G.F + plot.list$FG.F + plot.FG.G

ggsave("Figures/Final/Fig-S3a-talk.jpg", d.plot, width = 11, height = 4, dpi = 600)

###
# Extract normalized data
###

rld <- rlog(dds.16s.treat, blind = F) 
norm.rld <- assay(rld)
norm.count <- data.frame(cbind(ASV = rownames(norm.rld), norm.rld))

df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

norm.count <- norm.count[norm.count$ASV %in% unique(df.res$ASV),]

norm.counts <- merge(norm.count, data.frame(ASV = row.names(df.16s.tax), Family = df.16s.tax[,5]), by = "ASV", all.y = F)
norm.countl <- pivot_longer(norm.counts[,-1], col = -Family, names_to = "Sample", values_to = "norm.counts")

norm.countl$Sample <- gsub("\\.", "-", norm.countl$Sample)
norm.countl$norm.counts <- as.numeric(as.character(norm.countl$norm.counts))

# merge with grass data
norm.count.16s <- rbind(norm.countl, norm.countlg)
norm.count.16s <- unique(norm.count.16s)

#### DESeq (ITS) ####
ps.its.nocontrols.fam <- tax_glom(ps.ITS.nocontrols, taxrank = "Family", NArm = FALSE )

ps.its.nocontrols.nobss <- subset_samples(ps.its.nocontrols.fam, SampleSubType !="Background_soil_sand_mix")

###
# Grass baseline
###

sample_data(ps.its.nocontrols.nobss)$FunGroup <- relevel(as.factor(sample_data(ps.its.nocontrols.nobss)$FunGroup), "Grass")

treat.its <- phyloseq_to_deseq2(ps.its.nocontrols.nobss, ~ FunGroup)

dds.its.treat <- DESeq(treat.its, test = "Wald", fitType = "parametric")

resultsNames(dds.its.treat) #gives us comparisons for our contrasts

contrasts = c("F.G", "FG.G")

contrast.list <- list(F.G = c("FunGroup", "Forb", "Grass"),
                      FG.G = c("FunGroup", "grass_x_forb", "Grass"))


plot.name.list <- list(F.G = "Forb vs. Grass",
                       FG.G = "Competition vs. Grass")

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
 
   p <- ggplot(res.alpha.tax, aes(x = Family, y = log2FoldChange)) + 
    geom_point(size = 3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 12),
      legend.text = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "none",
      legend.key.size = unit(.35, "cm"),
      legend.title = element_text(size = 10),
      axis.title.x = element_blank()
     ) +
    labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change")))
     
  plot.list[[i]] = p
}

#plot results

plot.list$F.G+ plot.list$FG.G 
plot.FG.G <- plot.list$FG.G 

###
# Extract normalized data
###
rld <- rlog(dds.its.treat, blind = F) 
norm.rld <- assay(rld)
norm.countg <- data.frame(cbind(ASV = rownames(norm.rld), norm.rld))

res.alpha.tax$ASV <- row.names(res.alpha.tax)

norm.countg <- merge(res.alpha.tax, norm.countg, by = "ASV", all.y = F)
norm.countg <- norm.countg[,c(12,15:141)]
norm.countlg <- pivot_longer(norm.countg, col = -Family, names_to = "Sample", values_to = "norm.counts")

norm.countlg$Sample <- gsub("\\.", "-", norm.countlg$Sample)
norm.countlg$norm.counts <- as.numeric(as.character(norm.countlg$norm.counts))

###
# Forb baseline
###

sample_data(ps.its.nocontrols.nobss)$FunGroup <- relevel(as.factor(sample_data(ps.its.nocontrols.nobss)$FunGroup), "Forb")

treat.its <- phyloseq_to_deseq2(ps.its.nocontrols.nobss, ~ FunGroup)

dds.its.treat <- DESeq(treat.its, test = "Wald", fitType = "parametric")

resultsNames(dds.its.treat) #gives us comparisons for our contrasts

contrasts <- c("G.F", "FG.F")

contrast.list <- list(G.F = c("FunGroup", "Grass", "Forb"),
                      FG.F = c("FunGroup", "grass_x_forb", "Forb"))


plot.name.list <- list(G.F = "Grass vs. Forb",
                       FG.F = "Competition vs. Forb")

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
  res.alpha.tax$Family <- recode(res.alpha.tax$Family, Unclassified = "Unclassified Sebacinales")
 
  addline_format <- function(x,...){
    gsub('\\s','\n',x)
  }
  
  if(i == "FG.F") {
  p <- ggplot(res.alpha.tax, aes(x = Family, y = log2FoldChange)) + 
    geom_point(size = 3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 12),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "none",
      axis.title.x = element_blank()
     ) +
    labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change"))) +
    scale_x_discrete(breaks = unique(res.alpha.tax$Family), 
    labels = addline_format(c("Unclassified Sebacinales", "Tubeufiaceae")))
  
   plot.list[[i]] = p
  }
  
  if(i == "G.F") {
  p <- ggplot(res.alpha.tax, aes(x = Family, y = log2FoldChange)) + 
    geom_point(size = 3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 12),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "none",
      axis.title.x = element_blank()
     ) +
    labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change")))

   plot.list[[i]] = p
  }
   }

#plot results

des.p.its <- plot.list$G.F + plot.list$FG.F + plot.FG.G

#ggsave("Figures/Final/Fig-S3b.pdf", des.p.its, width = 9, height = 4, dpi = 600)

###
# Extract Normalized abundances
###
rld <- rlog(dds.its.treat, blind = F)
norm.rld <- assay(rld)
norm.count <- data.frame(cbind(ASV = rownames(norm.rld), norm.rld))

res.alpha.tax$ASV <- row.names(res.alpha.tax)

norm.counts <- merge(res.alpha.tax, norm.count, by = "ASV", all.y = F)
norm.counts <- norm.counts[,c(12,15:141)]

norm.countl <- pivot_longer(norm.counts, col = -Family, names_to = "Sample", values_to = "norm.counts")

norm.countl$Sample <- gsub("\\.", "-", norm.countl$Sample)
norm.countl$norm.counts <- as.numeric(as.character(norm.countl$norm.counts))

norm.count.its <- rbind(norm.countl, norm.countlg)

#### .---Fig: DESeq (16S - Main) ####
fam2 <- sort(c("Burkholderiaceae", "Methylophilaceae", "Fibrobacteraceae", "Veillonellaceae", "Clostridiaceae_1", "Rhodocyclaceae"))
fig.norm <- merge(norm.count.16s, mapping[,c(1,35)], by.x = "Sample", by.y = "SampleID_Fix")
fig.norm <- filter(fig.norm, Family %in% fam2)
fig.norm$FunGroup <- recode(fig.norm$FunGroup, grass_x_forb = "Competition")
#fig.norm.sum <- summarySE(fig.norm, groupvars = c("Family", "FunGroup"), measurevar = "norm.counts")
fig.norm.sum <- ddply(fig.norm, .(Family, FunGroup), summarize, norm.counts = max(norm.counts))

fig.norm.sum$cld <- c("a", "b", "a", # burkholderiaceae
                      "a", "b", "b", # Clostridiaceae
                      "a", "b", "b", # fibrobacteraceae
                      "a", "b", "b", # methylophilaceae
                      "a", "b", "a", # rhodocyclaceae
                      "a", "b", "b") # veillonellaceae
plotlist = list()

for(i in fam2) {
  p <- ggplot(fig.norm[fig.norm$Family == i,], aes(x = FunGroup, y = norm.counts, col = FunGroup)) +
    # geom_bar(stat = "identity", position = "dodge", col = "black") +
    # geom_errorbar(aes(ymin = norm.counts - se, ymax = norm.counts + se), width = 0.2, position = position_dodge(1), col = "black") +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.7) +
    labs(title = i, y = "Normalized Abundance") +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 9)
    ) +
    scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
    geom_text(data = fig.norm.sum[fig.norm.sum$Family == i,], mapping =
              aes(x = FunGroup, label = cld, y = norm.counts,
        size = 8,
        #vjust = 0.3),
        vjust = -0.7),
    col = "black") +
     ylim(NA, max(fig.norm[fig.norm$Family == i,]$norm.counts) + 0.055*max(fig.norm[fig.norm$Family == i,]$norm.counts))
  
  plotlist[[i]] <- p
}


p <- grid.arrange(arrangeGrob(grobs = plotlist, ncol = 3, left = textGrob("Normalized Abundance", rot = 90, vjust = 1)))

ggsave("Figures/Final/Fig-3b.jpg", p, width = 8, height = 5, units = "in", dpi = 600)


#### .---Fig: DESeq (16S - Sup) ####
fig.norm <- merge(norm.count.16s, mapping[,c(1,35)], by.x = "Sample", by.y = "SampleID_Fix")
fig.norm <- filter(fig.norm, !Family %in% fam2)
fig.norm$FunGroup <- recode(fig.norm$FunGroup, grass_x_forb = "Competition")
#fig.norm.sum <- summarySE(fig.norm, groupvars = c("Family", "FunGroup"), measurevar = "norm.counts")
fig.norm.sum <- ddply(fig.norm, .(Family, FunGroup), summarize, norm.counts = max(norm.counts))

fam3 <- unique(fig.norm.sum$Family)

fig.norm.sum$cld <- c("a", "ab", "b",
                      "a", "ab", "b",
                      "a", "ab", "b",
                      "a", "b", "ab",
                      "a", "b", "ab",
                      "a", "ab", "b",
                      "a", "ab", "b",
                      "a", "ab", "b",
                      "a", "ab", "b",
                      "a", "b", "ab",
                      "a", "ab", "b",
                      "a", "ab", "b",
                      "a", "ab", "b"
                      )

plotlist = list()

for(i in unique(fig.norm$Family)) {
  p <- ggplot(fig.norm[fig.norm$Family == i,], aes(x = FunGroup, y = norm.counts, col = FunGroup)) +
    # geom_bar(stat = "identity", position = "dodge", col = "black") + 
    # geom_errorbar(aes(ymin = norm.counts - se, ymax = norm.counts + se), width = 0.2, position = position_dodge(1), col = "black") +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.7) +
    labs(title = i, y = "Normalized Abundance") +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 9)
    ) +
    scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
    geom_text(data = fig.norm.sum[fig.norm.sum$Family == i,], mapping =
              aes(x = FunGroup, label = cld, y = norm.counts,
        size = 8,
        #vjust = 0.3),
        vjust = -0.7),
    col = "black") +
   ylim(NA, max(fig.norm[fig.norm$Family == i,]$norm.counts) + 0.1*max(fig.norm[fig.norm$Family == i,]$norm.counts))
  
  plotlist[[i]] <- p
}


p <- grid.arrange(arrangeGrob(grobs = plotlist, ncol = 3, left = textGrob("Normalized Abundance", rot = 90, vjust = 1)))

ggsave("Figures/Final/Fig-S4.pdf", p, width = 7, height = 9, units = "in", dpi = 600)

#### .---Fig: DESeq (ITS) ####
fig.norm <- merge(norm.count.its, mapping[,c(1,35)], by.x = "Sample", by.y = "SampleID_Fix")
fig.norm$FunGroup <- recode(fig.norm$FunGroup, grass_x_forb = "Competition")

#fig.norm.sum <- summarySE(fig.norm, groupvars = c("Family", "FunGroup"), measurevar = "norm.counts")

fig.norm.sum <- ddply(fig.norm, .(Family, FunGroup), summarize, norm.counts = max(norm.counts))

fig.norm.sum$cld <- c("a", "ab", "b",
                      "a", "b", "b",
                      "ab", "a", "b"
                      )

plotlist = list()

for(i in unique(fig.norm$Family)) {
  p <- ggplot(fig.norm[fig.norm$Family == i,], aes(x = FunGroup, y = norm.counts, col = FunGroup, group = FunGroup)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.7) +
    labs(title = i, y = "Normalized Abundance") +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 9)
    ) +
    scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
    geom_text(data = fig.norm.sum[fig.norm.sum$Family == i,], mapping =
              aes(x = FunGroup, label = cld, y = norm.counts,
        size = 8,
        vjust = -0.7),
    col = "black") +
    scale_y_continuous(
  labels = scales::number_format(accuracy = 1,
                                 decimal.mark = ','), limits = c(NA, max(fig.norm[fig.norm$Family == i,]$norm.counts) + 0.1*max(fig.norm[fig.norm$Family == i,]$norm.counts)))
  
  plotlist[[i]] <- p
}


p <- grid.arrange(arrangeGrob(grobs = plotlist, ncol = 3, left = textGrob("Normalized Abundance", rot = 90, vjust = 1)))

ggsave("Figures/Final/Fig-S5.pdf", p, width = 8, height = 2.5, units = "in", dpi = 600)


#### SourceTracker (16S) ####
asv.fam <- data.frame(ASV = row.names(df.16s.tax), Family = df.16s.tax$Family)
asv.fam <- filter(asv.fam, asv.fam$Family %in% c("Methylophilaceae", "Fibrobacteraceae", "Burkholderiaceae", "Veillonellaceae", "Rhodocyclaceae", "Clostridiaceae_1"))

ST_bg_16S <- read.delim("Data/16S/SourceTracker/sink_predictions_Background_soil_contributions.txt")
ST_bg_16S$Source <- "Background"

ST_uk_16S <- read.delim("Data/16S/SourceTracker/sink_predictions_Unknown_contributions.txt")
ST_uk_16S$Source <- "Unknown"

ST_forb_16S <- read.delim("Data/16S/SourceTracker/sink_predictions_Forb_contributions.txt")
sums <- rowSums(ST_forb_16S[,2:ncol(ST_forb_16S)])
ST_forb_16S$Source <- "Forb"
ST_grass_16S <- read.delim("Data/16S/SourceTracker/sink_predictions_Grass_contributions.txt")
ST_grass_16S$Source <- "Grass"

ST_16S <- rbind(ST_uk_16S, ST_bg_16S, ST_forb_16S, ST_grass_16S)

ST_16S <- pivot_longer(ST_16S, col = -c(SampleID, Source), names_to = "ASV", values_to = "prob")
ST_16S2 <- merge(ST_16S, asv.fam, by = "ASV", all.x = F)
test <- ddply(ST_16S2, .(Family), summarize, sum.prob = sum(prob))
ST_16S2 <- ddply(ST_16S2, .(Family, Source), summarize, source.prob = sum(prob))
ST_16S2 <- merge(ST_16S2, test, by = "Family")
ST_16S2$rel.prob <- ST_16S2$source.prob/ST_16S2$sum.prob

####.---Fig: SourceTracker (16S - M) ####
ST_16S2$Family <- factor(ST_16S2$Family, levels = c("Burkholderiaceae", "Clostridiaceae_1", "Fibrobacteraceae", "Methylophilaceae", "Rhodocyclaceae", "Veillonellaceae"))

p.st <- ggplot(ST_16S2, aes(y = rel.prob, x = Family, col = Source, fill = Source)) +
  geom_bar(position="stack", stat="identity", col = "black") + 
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_text(size = 9), 
    legend.text  = element_text(size = 8),
    legend.key.size = unit(.7, "lines"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  ) + 
  scale_fill_manual(values = c("#39568CFF", "magenta4", "#1F968BFF", "#B8DE29FF")) +
  labs(y = "Proportion of ASVs") +
  coord_flip()
  
ggsave("Figures/Final/Fig-4.jpg", p.st, width = 5, height = 2, dpi = 600, units = "in")

#### .---Fig: SourceTracker (16S - S) ####
asv.fam <- data.frame(ASV = row.names(df.16s.tax), Family = df.16s.tax$Family)
#fam4 <- sort(c(as.character(fam3), "Rhodocyclaceae", "Veillonellaceae", "Clostridiaceae_1"))
asv.fam <- filter(asv.fam, asv.fam$Family %in% fam3)

ST_16S3 <- merge(ST_16S, asv.fam, by = "ASV", all.x = F)
test <- ddply(ST_16S3, .(Family), summarize, sum.prob = sum(prob))
ST_16S3 <- ddply(ST_16S3, .(Family, Source), summarize, source.prob = sum(prob))
ST_16S3 <- merge(ST_16S3, test, by = "Family")
ST_16S3$rel.prob <- ST_16S3$source.prob/ST_16S3$sum.prob
#ST_16S3$Family <- factor(ST_16S3$Family, levels = c("Veillonellaceae", "Flavobacteriaceae", "Rhodocyclaceae", "Rubritaleaceae", "Weeksellaceae", "Geodermatophilaceae", "Propionibacteriaceae", "Streptomycetaceae", "Nitrospiraceae", "Azospirillaceae", "Beijerinckiaceae", "Microbacteriaceae", "Steroidobacteraceae", "Acetobacteraceae", "Cyclobacteriaceae"))
  
p.st <- ggplot(ST_16S3, aes(y = rel.prob, x = Family, col = Source, fill = Source)) +
  geom_bar(position="stack", stat="identity", col = "black") + 
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_text(size = 9), 
    legend.text  = element_text(size = 8),
    legend.key.size = unit(.7, "lines"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  ) + 
  scale_fill_manual(values = c("#39568CFF", "magenta4", "#1F968BFF", "#B8DE29FF")) +
  labs(y = "Proportion of ASVs") +
  coord_flip()
  
ggsave("Figures/Final/Fig-S6.jpg", p.st, width = 6, height = 4, dpi = 600, units = "in")

ST_16S <- rbind(ST_16S2, ST_16S3)
ST_16S <- filter(ST_16S, Source != "Background", Source != "Unknown")
ST_16S <- pivot_wider(ST_16S[,c(1,2,5)], names_from = "Source", values_from = "rel.prob")
ST_16S$dom.source <- ifelse(ST_16S$Forb > ST_16S$Grass, "forb", "grass")

#### SourceTracker (ITS) ####
fam.ITS <- c("Tubeufiaceae", "Ceratobasidiaceae", "Unclassified Sebacinales")

asv.fam <- data.frame(ASV = row.names(df.ITS.tax), Order = df.ITS.tax$Order, Family = df.ITS.tax$Family)

levels(asv.fam$Family) <- c(levels(asv.fam$Family), "Unclassified Sebacinales")

asv.fam[asv.fam$Order == "Sebacinales" & asv.fam$Family == "Unclassified",]$Family <- "Unclassified Sebacinales"

asv.fam <- filter(asv.fam, asv.fam$Family %in% fam.ITS)

ST_bg_its <- read.delim("Data/ITS/SourceTracker/sink_predictions_Background_soil_contributions.txt")
ST_bg_its$Source <- "Background"

ST_uk_its <- read.delim("Data/ITS/SourceTracker/sink_predictions_Unknown_contributions.txt")
ST_uk_its$Source <- "Unknown"

ST_forb_its <- read.delim("Data/ITS/SourceTracker/sink_predictions_Forb_contributions.txt")
sums <- rowSums(ST_forb_its[,2:ncol(ST_forb_its)])
ST_forb_its$Source <- "Forb"
ST_grass_its <- read.delim("Data/ITS/SourceTracker/sink_predictions_Grass_contributions.txt")
ST_grass_its$Source <- "Grass"

ST_its <- rbind(ST_uk_its, ST_bg_its, ST_forb_its, ST_grass_its)

ST_its <- pivot_longer(ST_its, col = -c(SampleID, Source), names_to = "ASV", values_to = "prob")
ST_its <- merge(ST_its, asv.fam, by = "ASV", all.x = F)
test <- ddply(ST_its, .(Family), summarize, sum.prob = sum(prob))
ST_its <- ddply(ST_its, .(Family, Source), summarize, source.prob = sum(prob))
ST_its <- merge(ST_its, test, by = "Family")
ST_its$rel.prob <- ST_its$source.prob/ST_its$sum.prob

#### .---Fig: SourceTracker (ITS) ####
p.st <- ggplot(ST_its, aes(y = rel.prob, x = Family, col = Source, fill = Source)) +
  geom_bar(position="stack", stat="identity", col = "black") + 
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.title = element_text(size = 9), 
    legend.text  = element_text(size = 8),
    legend.key.size = unit(.7, "lines"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  ) + 
  scale_fill_manual(values = c("#39568CFF", "magenta4", "#1F968BFF", "#B8DE29FF")) +
  labs(y = "Proportion of ASVs") +
  coord_flip()
  
ggsave("Figures/Final/Fig-S7.pdf", p.st, width = 5, height = 1.5, dpi = 600, units = "in")

#### Log Response Ratio ####
#mapping$FunGroup <- as.factor(mapping$FunGroup)

# Log response ratio
  
map.comp <- mapping %>%
  filter(Competion == "TwoSpecies") %>%
  pivot_longer(c(grass, forb), names_to = "FunGroup2", values_to = "Species")

biomass <- mapping %>%
  filter(Competion == "SingleSpecies") %>%
  filter(TreatRep != "LA3") %>%
  dplyr::group_by(PlantSpeciesSampled) %>%
  dplyr::summarise(avg.wt = mean(Weight.g, na.rm = T)) %>%
  dplyr::rename(Species = PlantSpeciesSampled) %>%
  left_join(map.comp, by = "Species") %>%
  filter(Competion == "TwoSpecies")

biomass$Competitor <- ifelse(biomass$FunGroup2 == "forb", 
                             substr(biomass$Treatment, start = 1, stop = 2),
                             substr(biomass$Treatment, start = 3, stop = 4))

biomass$weight.d <- ifelse(biomass$FunGroup2 == "forb", log(biomass$forb.wt/biomass$avg.wt), log(biomass$grass.wt/biomass$avg.wt))


m.bio <- lmer(LRR ~ FunGroup2 + (1|Treatment) + (1|Treatment:Species), biomass)
plot(fitted(m.bio), resid(m.bio))
qqnorm(resid(m.bio))
qqline(resid(m.bio))
summary(m.bio)

#### .---Fig: LRR ####
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
  scale_fill_manual(values = c("magenta4", "#1F968BFF"))

#ggsave("Figures/Final/Fig-5.pdf", bio.plot, width = 5, height = 4, units = "in", dpi = 600)

#### LRR v NA (16S) ####
lrr.fam.df <- psmelt(ps.16s.nocontrols.fam)

lrr.fam.comp <- lrr.fam.df %>%
  filter(Competion == "TwoSpecies") %>%
  pivot_longer(c(grass, forb), names_to = "FunGroup2", values_to = "Species")

lrr.fam.16s <- lrr.fam.df %>%
  dplyr::filter(Competion == "SingleSpecies") %>%
  dplyr::group_by(PlantSpeciesSampled, Family) %>%
  dplyr::summarise(avg.wt = mean(Weight.g, na.rm = T)) %>%
  dplyr::rename(Species = PlantSpeciesSampled) %>%
  left_join(lrr.fam.comp, by = c("Species", "Family")) %>%
  filter(Competion == "TwoSpecies") 

lrr.fam.16s$weight.d <- ifelse(lrr.fam.16s$FunGroup2 == "forb", log(lrr.fam.16s$forb.wt/lrr.fam.16s$avg.wt), log(lrr.fam.16s$grass.wt/lrr.fam.16s$avg.wt))

lrr.fam.16s <- merge(norm.count.16s, lrr.fam.16s, by = c("Sample", "Family"), all.y = F)

lrr.fam.16s.spp <- lrr.fam.16s

#### LRR v NA (ITS) ####
ITS.fam <- psmelt(ps.its.nocontrols.fam)

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
  
ITS.fam.df3$weight.d <- ifelse(ITS.fam.df3$FunGroup2 == "forb", log(ITS.fam.df3$forb.wt/ITS.fam.df3$avg.wt), log(ITS.fam.df3$grass.wt/ITS.fam.df3$avg.wt))

levels(ITS.fam.df3$Family) <- c(levels(ITS.fam.df3$Family), "Unclassified Sebacinales")

ITS.fam.df3[ITS.fam.df3$Order == "Sebacinales" & ITS.fam.df3$Family == "Unclassified",]$Family <- "Unclassified Sebacinales"

lrr.fam.ITS <- merge(norm.count.its, ITS.fam.df3, by = c("Sample", "Family"), all.y = F)
lrr.fam.ITS.spp <- lrr.fam.ITS

#### Models: LRR v NA (16S) - M ####
fam2 <- c("Burkholderiaceae", "Methylophilaceae", "Fibrobacteraceae", "Veillonellaceae", "Clostridiaceae_1", "Rhodocyclaceae")

Species <- unique(biomass$Species)
group <- c("grass", "forb")
fam.16s <- c(fam2, as.character(fam3))
lrr.16s.m <- expand.grid(Family = fam.16s, FunGroup2 = group, est = NA, se = NA, p = NA, r2c = NA, r2m = NA)

for(i in fam.16s) {
  for(j in group){
  tmp <- lmer(weight.d ~ norm.counts + (1|Species) + (1|Treatment), data = lrr.fam.16s[lrr.fam.16s$Family == i & lrr.fam.16s$FunGroup2 == j,])
 
   lrr.16s.m[lrr.16s.m$Family == i & lrr.16s.m$FunGroup2 == j, "est"] <- summary(tmp)[["coefficients"]][2,1]
  
  lrr.16s.m[lrr.16s.m$Family == i & lrr.16s.m$FunGroup2 == j, "se"] <- summary(tmp)[["coefficients"]][2,2]
    
  lrr.16s.m[lrr.16s.m$Family == i & lrr.16s.m$FunGroup2 == j, "p"] <- summary(tmp)[["coefficients"]][2,5]
       
  lrr.16s.m[lrr.16s.m$Family == i & lrr.16s.m$FunGroup2 == j, "r2c"] <- as.numeric(performance::r2(tmp)$R2_conditional)
  
  lrr.16s.m[lrr.16s.m$Family == i & lrr.16s.m$FunGroup2 == j, "r2m"] <- as.numeric(performance::r2(tmp)$R2_marginal)
      
  }
}

lrr.16s.m$lines <- ifelse(lrr.16s.m$p < 0.05, "dashed", "solid") # this should give me the exact opposite but for some reason it's switched in the graph... 


#### Models: LRR v NA (ITS) ####
fam.ITS <- c("Tubeufiaceae", "Ceratobasidiaceae", "Unclassified Sebacinales")

group <- c("grass", "forb")
lrr.ITS.m <- expand.grid(Family = fam.ITS, FunGroup2 = group, est = NA, se = NA, p = NA)

for(i in fam.ITS) {
  for(j in group){
  tmp <- lmer(weight.d ~ norm.counts + (1|Treatment), data = lrr.fam.ITS[lrr.fam.ITS$Family == i & lrr.fam.ITS$FunGroup2 == j,])
  
  lrr.ITS.m[lrr.ITS.m$Family == i & lrr.ITS.m$FunGroup2 == j, "est"] <- summary(tmp)[["coefficients"]][2,1]
    lrr.ITS.m[lrr.ITS.m$Family == i & lrr.ITS.m$FunGroup2 == j, "se"] <- summary(tmp)[["coefficients"]][2,2]
      lrr.ITS.m[lrr.ITS.m$Family == i & lrr.ITS.m$FunGroup2 == j, "p"] <- summary(tmp)[["coefficients"]][2,5]
  }
}

lrr.ITS.m$lines <- ifelse(lrr.ITS.m$p < 0.05, "dashed", "solid")  # switched
lrr.fam.ITS <- merge(lrr.fam.ITS, lrr.ITS.m, by = c("Family", "FunGroup2"))



#### .---Fig: LRR v NA (16S - M) ####
fam2 <- c("Burkholderiaceae", "Methylophilaceae", "Fibrobacteraceae", "Clostridiaceae_1")

lrr.fam.16s <- merge(lrr.fam.16s, lrr.16s.m, by = c("Family", "FunGroup2"))

plotlist = list()

for(i in fam2) {
  p <- ggplot(lrr.fam.16s[lrr.fam.16s$Family == i,], aes(x = norm.counts, y = weight.d, col = FunGroup2, group = FunGroup2, linetype = FunGroup2)) + 
  geom_smooth(method = "lm", formula = y ~ x, aes(linetype = lines)) +
  geom_point(size = .8) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "none",
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  ) +
  labs(title = i, x = "Normalized Abundance", y = "Log Response Ratio") +
  scale_color_manual(values = c("magenta4", "#1F968BFF"))

  plotlist[[i]] <- p
}

legend <- gtable_filter(ggplotGrob(
  ggplot(lrr.fam.16s[lrr.fam.16s$Family == i,], aes(x = norm.counts, y = weight.d, col = FunGroup2, group = FunGroup2)) +
  geom_point(size = .8) +
  theme(
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("magenta4", "#1F968BFF"))
    ), "guide-box") 

p2 <- grid.arrange(arrangeGrob(grobs = plotlist, 
                               ncol = 2, 
                               left = textGrob("Log Response Ratio", 
                                               rot = 90, 
                                               vjust = 1)),  
                   legend, 
                   widths = unit.c(unit(1, "npc") - legend$width, legend$width), 
                   nrow=1)

ggsave("Figures/Final/Fig-6.pdf", p2, width = 6, height = 5, units = "in", dpi = 600)

#### .---Fig: LRR v NA (16S - S) ####
plotlist = list()
fam4 <- c(as.character(fam3), "Rhodocyclaceae", "Veillonellaceae")

for(i in fam4) {
  p <- ggplot(lrr.fam.16s[lrr.fam.16s$Family == i,], aes(x = norm.counts, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
  geom_smooth(method = "lm", formula = y ~ x, linetype = "dashed") +
  geom_point(size = .8) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "none",
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  ) +
  labs(title = i, x = "Normalized Abundance", y = "Log Response Ratio") +
  scale_color_manual(values = c("magenta4", "#1F968BFF"))
  
  plotlist[[i]] <- p
}

legend <- gtable_filter(ggplotGrob(
  ggplot(lrr.fam.16s[lrr.fam.16s$Family == i,], aes(x = norm.counts, y = weight.d, col = FunGroup2, group = FunGroup2)) +
  geom_point(size = .8) +
  theme(
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("magenta4", "#1F968BFF"))
    ), "guide-box") 

p2 <- grid.arrange(arrangeGrob(grobs = plotlist, 
                               ncol = 3, 
                               left = textGrob("Log Response Ratio", 
                                               rot = 90, 
                                               vjust = 1)),  
                   legend, 
                   widths = unit.c(unit(1, "npc") - legend$width, legend$width), 
                   nrow=1)

ggsave("Figures/Final/Fig-S8.pdf", p2, width = 7, height = 9, units = "in", dpi = 600)


#### .---Fig: LRR v NA (ITS) ####

plotlist = list()

for(i in fam.ITS) {
  p <- ggplot(lrr.fam.ITS[lrr.fam.ITS$Family == i,], aes(x = norm.counts, y = weight.d, col = FunGroup2, group = FunGroup2, linetype = FunGroup2)) + 
  geom_smooth(method = "lm", formula = y ~ x, aes(linetype = lines)) +
  geom_point(size = .8) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "none",
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  ) +
  labs(title = i, x = "Normalized Abundance", y = "Log Response Ratio") +
  scale_color_manual(values = c("magenta4", "#1F968BFF")) +
  scale_x_continuous(
  labels = scales::number_format(accuracy = 1,
                                 decimal.mark = ','))
  plotlist[[i]] <- p
}

legend <- gtable_filter(ggplotGrob(
  ggplot(lrr.fam.ITS[lrr.fam.ITS$Family == i,], aes(x = norm.counts, y = weight.d, col = FunGroup2, group = FunGroup2)) +
  geom_point(size = .8) +
  theme(
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("magenta4", "#1F968BFF"))
    ), "guide-box") 

p2 <- grid.arrange(arrangeGrob(grobs = plotlist, 
                               ncol = 3, 
                               left = textGrob("Log Response Ratio", 
                                               rot = 90, 
                                               vjust = 1,
                                               gp = gpar(fontsize = 10))),  
                   legend, 
                   widths = unit.c(unit(1, "npc") - legend$width, legend$width), 
                   nrow=1)

ggsave("Figures/Final/Fig-S9.pdf", p2, width = 7, height = 2, units = "in", dpi = 600)


#### REVISIONS ####

#### Number of bacterial families ####
bac.fam <- unique(df.16s.tax$Family) # 303 unique families
fun.fam <- unique(df.ITS.tax$Family) # 232 unique families

nrow(df.16s.tax[df.16s.tax$Family == "Unknown_Family",])/nrow(df.16s.tax)
nrow(df.ITS.tax[df.ITS.tax$Family == "Unclassified",])/nrow(df.ITS.tax)

#### Sample numbers ####
GL_Alpha <- filter(GL_Alpha, SampleSubType == "Rhizosphere_soil")

nrow(GL_Alpha[GL_Alpha$FunGroup == "Grass",])
nrow(GL_Alpha[GL_Alpha$FunGroup == "Forb",])
nrow(GL_Alpha[GL_Alpha$FunGroup == "grass_x_forb",])

#### Forb LRR vs Grass biomass ####
a <- ggplot(biomass[biomass$FunGroup2 != "forb",], aes(x = forb.wt, y = weight.d)) +
  theme_classic() +
  geom_smooth(method = "lm", linetype = "dashed") +
  geom_point() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  ) +
  labs(x = "Forb weight (g)", y = "Grass Biomass Log Response Ratio")

m.lb <- lmer(weight.d ~ forb.wt + (1|Species), biomass[biomass$FunGroup2 != "forb",])
plot(fitted(m.lb), resid(m.lb))
qqnorm(resid(m.lb))
qqline(resid(m.lb))
summary(m.lb)

#### Grass LRR vs Forb biomass ####
b <- ggplot(biomass[biomass$FunGroup2 != "grass",], aes(x = grass.wt, y = weight.d)) +
  theme_classic() +
  geom_smooth(method = "lm", linetype = "dashed") +
  geom_point() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  ) +
  labs(x = "Grass weight (g)", y = "Forb Biomass Log Response Ratio")

m.lb2 <- lmer(weight.d ~ grass.wt + (1|Species), biomass[biomass$FunGroup2 != "grass",])
plot(fitted(m.lb2), resid(m.lb2))
qqnorm(resid(m.lb2))
qqline(resid(m.lb2))
summary(m.lb2)

# ggsave("Figures/Final/Revisions/bio-lrr-g.jpg", a, width = 6, height = 5, units = "in", dpi = 600)
# ggsave("Figures/Final/Revisions/bio-lrr-f.jpg", b, width = 6, height = 5, units = "in", dpi = 600)

#### Grass wt vs forb wt ####
m.wt <- lmer(log(forb.wt) ~ log(grass.wt) + (1|Species), biomass[biomass$FunGroup2 != "grass" & biomass$Competion == "TwoSpecies",])
plot(fitted(m.wt), resid(m.wt))
qqnorm(resid(m.wt))
qqline(resid(m.wt))
summary(m.wt)

c <- ggplot(biomass[biomass$FunGroup2 != "grass" & biomass$Competion == "TwoSpecies",], aes(x = grass.wt, y = forb.wt)) +
  geom_smooth(method = "lm", linetype = "dashed") +
  geom_point() +
  labs(x = "Grass Biomass (g)", y = "Forb Biomass (g)") +
  theme_classic() + 
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  ) 

ggsave("Figures/Final/Revisions/bio-corr.jpg", c, width = 6, height = 5, units = "in", dpi = 600)

#### FunGroup Dominance ####
mb3 <- lmer(log(Weight.g) ~ FunGroup + (1|PlantSpeciesSampled), mapping[mapping$Competion == "SingleSpecies" & mapping$PlantSpeciesSampled != "none",])
plot(fitted(mb3), resid(mb3))
qqnorm(resid(mb3))
qqline(resid(mb3))
summary(mb3)

biomass$Weight.g <- ifelse(biomass$FunGroup2 == "grass", biomass$grass.wt, biomass$forb.wt)
mb4 <- lmer(log(Weight.g) ~ FunGroup2 + (1|PlantSpeciesSampled), biomass[biomass$Competion == "TwoSpecies" & biomass$PlantSpeciesSampled != "none",])
mb4.1 <- lmer(log(Weight.g) ~ FunGroup2 + (1|PlantSpeciesSampled) + (1|TreatRep), biomass[biomass$Competion == "TwoSpecies" & biomass$PlantSpeciesSampled != "none",])
anova(mb4, mb4.1)

plot(fitted(mb4), resid(mb4))
qqnorm(resid(mb4))
qqline(resid(mb4))
summary(mb4)

### Graphs
bio.sum <- summarySE(mapping[mapping$Competion == "SingleSpecies" & mapping$PlantSpeciesSampled != "none",], measurevar = "Weight.g", groupvars = "FunGroup", na.rm = T)

a <- ggplot(bio.sum[complete.cases(bio.sum),], aes(y = Weight.g, x = FunGroup, fill = FunGroup)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Weight.g - se, ymax = Weight.g + se, width = 0.2)) +
  xlab("biomass (g)") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 13)#,
    #axis.title.y = element_text(size = 13)
  ) +
  ylim(0,0.19) +
  scale_fill_manual(values = c("magenta4", "#1F968BFF"))

bio.sum2 <- summarySE(biomass[biomass$Competion == "TwoSpecies" & biomass$PlantSpeciesSampled != "none",], measurevar = "Weight.g", groupvars = "FunGroup2", na.rm = T)

b <- ggplot(bio.sum2, aes(y = Weight.g, x = FunGroup2, fill = FunGroup2)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Weight.g - se, ymax = Weight.g + se, width = 0.2)) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 13)
  ) +
  ylim(0,0.19) +
  scale_fill_manual(values = c("magenta4", "#1F968BFF"))

p2 <- grid.arrange(arrangeGrob(grobs = list(a,b), 
                               ncol = 2, 
                               left = textGrob("biomass (g)", 
                                               rot = 90, 
                                               vjust = 1), size = 13), 
                   nrow = 1)

#ggsave("Figures/Final/Fig-talk.jpg", p2, width = 8, height = 4, dpi = 600)

#### Models: P Biomass v NA (16s) ####
bio.16s <- merge(norm.count.16s, mapping[mapping$Competion == "TwoSpecies",], by.x = "Sample", by.y = "SampleID_Fix")

# GRASSES

bio.16s.g <- expand.grid(Family = fam.16s, est = NA, se = NA, p = NA, r2.c = NA, r2.m = NA)

for(i in fam.16s) { 
  tmp <- lmer(log(grass.wt) ~ norm.counts + (1|grass), data = bio.16s[bio.16s$Family == i,])
  bio.16s.g[bio.16s.g$Family == i, "est"] <- summary(tmp)[["coefficients"]][2,1]
  bio.16s.g[bio.16s.g$Family == i, "se"] <- summary(tmp)[["coefficients"]][2,2]
  bio.16s.g[bio.16s.g$Family == i, "p"] <- summary(tmp)[["coefficients"]][2,5]
  bio.16s.g[bio.16s.g$Family == i, "r2.m"] <- r.squaredGLMM(tmp)[1]
  bio.16s.g[bio.16s.g$Family == i, "r2.c"]  <- r.squaredGLMM(tmp)[2]
}

bio.16s.g$lines <- ifelse(bio.16s.g$p < 0.05, "dashed", "solid") 
bio.16s.g <- merge(bio.16s[,-c(4:28,30,32,34)], bio.16s.g, by = c("Family"))
bio.16s.g$FunGroup <- "grass"
colnames(bio.16s.g)[4] <- "Species"
colnames(bio.16s.g)[6] <- "Weight.g"

# Forbs

bio.16s.f <- expand.grid(Family = fam.16s, est = NA, se = NA, p = NA, r2.c = NA, r2.m = NA)

for(i in fam.16s) { 
  tmp <- lmer(log(forb.wt) ~ norm.counts + (1|forb), data = bio.16s[bio.16s$Family == i,])
  bio.16s.f[bio.16s.f$Family == i, "est"] <- summary(tmp)[["coefficients"]][2,1]
  bio.16s.f[bio.16s.f$Family == i, "se"] <- summary(tmp)[["coefficients"]][2,2]
  bio.16s.f[bio.16s.f$Family == i, "p"] <- summary(tmp)[["coefficients"]][2,5]
  bio.16s.f[bio.16s.f$Family == i, "r2.m"] <- r.squaredGLMM(tmp)[1]
  bio.16s.f[bio.16s.f$Family == i, "r2.c"]  <- r.squaredGLMM(tmp)[2]
}

bio.16s.f$lines <- ifelse(bio.16s.f$p < 0.05, "dashed", "solid")

bio.16s.f <- merge(bio.16s[,-c(4:28,33,34)], bio.16s.f, by = c("Family"))
bio.16s.f$FunGroup <- "forb"
bio.16s.f <- bio.16s.f[,-4]
colnames(bio.16s.f)[4] <- "Species"
colnames(bio.16s.f)[6] <- "Weight.g"

bio.16s.fg <- rbind(bio.16s.f, bio.16s.g)

####.---Fig: Bio v NA (16S) - M ####
fam.m <- sort(c("Methylophilaceae", "Burkholderiaceae", "Fibrobacteraceae", "Clostridiaceae_1", "Veillonellaceae", "Rhodocyclaceae"))

plotlist = list()

for(i in fam.m) {
  p <- ggplot(bio.16s.fg[bio.16s.fg$Family == i,], aes(x = norm.counts, y = log(Weight.g), col = FunGroup, group = FunGroup, linetype = FunGroup)) + 
  geom_smooth(data = bio.16s.fg[bio.16s.fg$Family == i & bio.16s.fg$p < 0.10,], method = "lm", formula = y ~ x, aes(linetype = lines)) +
  geom_point(size = .6) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "none",
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.text = element_text(size = 6)
  ) +
    labs(title = i, x = "Normalized Abundance", y = "Log Biomass (g)") +
  scale_color_manual(values = c("magenta4", "#1F968BFF"))

  plotlist[[i]] <- p
}

legend <- gtable_filter(ggplotGrob(
  ggplot(bio.16s.fg[bio.16s.fg$Family == i,], aes(x = norm.counts, y = log(Weight.g), col = FunGroup, group = FunGroup)) +
  geom_point(size = .8) +
  theme(
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("magenta4", "#1F968BFF"))
    ), "guide-box") 

p2 <- grid.arrange(arrangeGrob(grobs = plotlist, 
                               ncol = 3, 
                               left = textGrob("log biomass in pairs (g)", 
                                               rot = 90, 
                                               vjust = 1), 
                               bottom = textGrob("normalized abundance in pairs", hjust = 0.5)),  
                   legend, 
                   widths = unit.c(unit(1, "npc") - legend$width, legend$width), 
                   nrow = 1)


ggsave("Figures/Final/Revisions/Fig-S14.jpg", p2, width = 7, height = 4, units = "in", dpi = 600)

####.---Fig: Bio v NA (16S) - S ####
plotlist = list()

for(i in sort(fam.16s[!fam.16s %in% fam.m])) {
  p <- ggplot(bio.16s.fg[bio.16s.fg$Family == i,], aes(x = norm.counts, y = log(Weight.g), col = FunGroup, group = FunGroup, linetype = FunGroup)) + 
  geom_smooth(data = bio.16s.fg[bio.16s.fg$Family == i & bio.16s.fg$p < 0.10,], method = "lm", formula = y ~ x, aes(linetype = lines)) +
  geom_point(size = .8) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "none",
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  ) +
    labs(title = i, x = "Normalized Abundance", y = "Log Biomass (g)") +
  scale_color_manual(values = c("magenta4", "#1F968BFF"))

  plotlist[[i]] <- p
}

legend <- gtable_filter(ggplotGrob(
  ggplot(bio.16s.fg[bio.16s.fg$Family == i,], aes(x = norm.counts, y = log(Weight.g), col = FunGroup, group = FunGroup)) +
  geom_point(size = .8) +
  theme(
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("magenta4", "#1F968BFF"))
    ), "guide-box") 

p2 <- grid.arrange(arrangeGrob(grobs = plotlist, 
                               ncol = 3, 
                               left = textGrob("log biomass (g)", 
                                               rot = 90, 
                                               vjust = 1)),  
                   legend, 
                   widths = unit.c(unit(1, "npc") - legend$width, legend$width), 
                   nrow=1)

ggsave("Figures/Final/Revisions/Fig-S8.jpg", p2, width = 7, height = 9, units = "in", dpi = 600)

#### Models: P Biomass v NA (ITS) ####
bio.its <- merge(norm.count.its, mapping[mapping$Competion == "TwoSpecies",], by.x = "Sample", by.y = "SampleID_Fix")

# GRASSES
fam.ITS <- c("Tubeufiaceae", "Ceratobasidiaceae", "Unclassified Sebacinales")
bio.its.g <- expand.grid(Family = fam.ITS, est = NA, se = NA, p = NA, r2.c = NA, r2.m = NA)

for(i in fam.ITS) { 
  tmp <- lmer(log(grass.wt) ~ norm.counts + (1|grass), data = bio.its[bio.its$Family == i,])
  bio.its.g[bio.its.g$Family == i, "est"] <- summary(tmp)[["coefficients"]][2,1]
  bio.its.g[bio.its.g$Family == i, "se"] <- summary(tmp)[["coefficients"]][2,2]
  bio.its.g[bio.its.g$Family == i, "p"] <- summary(tmp)[["coefficients"]][2,5]
  bio.its.g[bio.its.g$Family == i, "r2.m"] <- r.squaredGLMM(tmp)[1]
  bio.its.g[bio.its.g$Family == i, "r2.c"]  <- r.squaredGLMM(tmp)[2]
}

bio.its.g$lines <- ifelse(bio.its.g$p < 0.05, "dashed", "solid") 
bio.its.g <- merge(bio.its[,-c(4:28,30,32,34)], bio.its.g, by = c("Family"))
bio.its.g$FunGroup <- "grass"
colnames(bio.its.g)[4] <- "Species"
colnames(bio.its.g)[6] <- "Weight.g"

# Forbs

bio.its.f <- expand.grid(Family = fam.ITS, est = NA, se = NA, p = NA, r2.c = NA, r2.m = NA)

for(i in fam.ITS) { 
  tmp <- lmer(log(forb.wt) ~ norm.counts + (1|forb), data = bio.its[bio.its$Family == i,])
  bio.its.f[bio.its.f$Family == i, "est"] <- summary(tmp)[["coefficients"]][2,1]
  bio.its.f[bio.its.f$Family == i, "se"] <- summary(tmp)[["coefficients"]][2,2]
  bio.its.f[bio.its.f$Family == i, "p"] <- summary(tmp)[["coefficients"]][2,5]
  bio.its.f[bio.its.f$Family == i, "r2.m"] <- r.squaredGLMM(tmp)[1]
  bio.its.f[bio.its.f$Family == i, "r2.c"]  <- r.squaredGLMM(tmp)[2]
}

bio.its.f$lines <- ifelse(bio.its.f$p < 0.05, "dashed", "solid")

bio.its.f <- merge(bio.its[,-c(4:28,33,34)], bio.its.f, by = c("Family"))
bio.its.f$FunGroup <- "forb"
bio.its.f <- bio.its.f[,-4]
colnames(bio.its.f)[4] <- "Species"
colnames(bio.its.f)[6] <- "Weight.g"

bio.its.fg <- rbind(bio.its.f, bio.its.g)

####.---Fig: Bio v NA (ITS) ####

plotlist = list()

for(i in fam.ITS) {
  p <- ggplot(bio.its.fg[bio.its.fg$Family == i,], aes(x = norm.counts, y = log(Weight.g), col = FunGroup, group = FunGroup, linetype = FunGroup)) + 
  geom_smooth(data = bio.its.fg[bio.its.fg$Family == i & bio.its.fg$p < 0.10,], method = "lm", formula = y ~ x, aes(linetype = lines)) +
  geom_point(size = .8) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "none",
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  ) +
  labs(title = i, x = "Normalized Abundance", y = "Log Biomass (g)") +
  scale_color_manual(values = c("magenta4", "#1F968BFF"))

  plotlist[[i]] <- p
}

legend <- gtable_filter(ggplotGrob(
  ggplot(bio.its.fg[bio.its.fg$Family == i,], aes(x = norm.counts, y = log(Weight.g), col = FunGroup, group = FunGroup)) +
  geom_point(size = .8) +
  theme(
    legend.title = element_blank()) 
   +
   scale_color_manual(values = c("magenta4", "#1F968BFF"))
     ), "guide-box") 

p2 <- grid.arrange(arrangeGrob(grobs = plotlist, 
                               ncol = 3, 
                              left = textGrob("Log Biomass (g)", 
                                              rot = 90, 
                                              vjust = 1)), 
                  legend, 
                   widths = unit.c(unit(1, "npc") - legend$width, legend$width),
                   nrow=1)

ggsave("Figures/Final/Revisions/Fig-S9.jpg", p2, width = 7, height = 2, units = "in", dpi = 600)



#### Relative Abundance ####
ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion)))

df_ASV <- psmelt(ps.16s.nocontrols.rare)

grouped_ASV2 <- ddply(df_ASV, .(Family, Sample), summarize, RA = sum(Abundance)/9434)

se <- function(x) sqrt(var(x)/length(x))
avgs_ASV2 <- ddply(grouped_ASV2, .(Family), summarize, Abundance = mean(RA)*100, sd = sd(RA)*100) # RA of all families (over fungroups)

Table3 <- avgs_ASV2[avgs_ASV2$Family %in% fam.16s,]

## ITS
ps.ITS.nocontrols.rare <- subset_samples(ps.ITS.nocontrols.rare, !(is.na(Competion)))
df_ASV.its <- psmelt(ps.ITS.nocontrols.rare)

grouped_ASV.its2 <- ddply(df_ASV.its, .(Order, Family, Sample), summarize, RA = sum(Abundance)/7557)

# se <- function(x) sqrt(var(x)/length(x))
avgs_ASV.its2 <- ddply(grouped_ASV.its2, .(Order, Family), summarize, Abundance = mean(RA)*100, se = se(RA)*100) # RA of all families (over fungroups)

Table3 <- avgs_ASV.its2[avgs_ASV.its2$Family %in% c(f_Ceratobasidiaceae, f_ubeufiaceae, f_),]
