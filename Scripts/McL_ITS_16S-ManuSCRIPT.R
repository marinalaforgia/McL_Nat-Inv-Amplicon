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


# Load 
ps.ITS.nocontrols.rare <- readRDS("Data/ITS/greenhouse_its_itsx_decontam_controlsremoved_rare7557.rds") # rarefied ITS data
ps.ITS.nocontrols.rare <- readRDS("Data/greenhouse_its_itsx_decontam_controlsremoved_rare7557.rds") # rarefied ITS data, Cassie's path


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

#### Ordination (16S) ####
ps.16s.nocontrols.rare <- subset_samples(ps.16s.nocontrols.rare, !(is.na(Competion)))

sample_data(ps.16s.nocontrols.rare)$FunGroup <- recode(sample_data(ps.16s.nocontrols.rare)$FunGroup, grass_x_forb = "Competition")

sample_data(ps.16s.nocontrols.rare)$FunGroup <- factor(sample_data(ps.16s.nocontrols.rare)$FunGroup, levels = c("Forb", "Grass", "Competition"))
  
ps.rare.ord.tr <- ordinate(ps.16s.nocontrols.rare, "PCoA", "wunifrac")

DistWU <- phyloseq::distance(ps.16s.nocontrols.rare, method = "wunifrac", type = "samples")

# adonis compares centroids and is sensitive to differences in spread

adonis(DistWU ~ FunGroup, as(sample_data(ps.16s.nocontrols.rare), "data.frame"), permutations = 9999)

C <- betadisper(DistWU, as(sample_data(ps.16s.nocontrols.rare), "data.frame")$FunGroup)
permutest(C, permutations = 9999) # spread doesnt vary between treatments

pair.16s.com <- adonis.pair(DistWU, as(sample_data(ps.16s.nocontrols.rare), "data.frame")$FunGroup, nper = 9999, corr.method = "BH")
kable(pair.16s.com)

#### Ordination (ITS) ####
ps.ITS.nocontrols.rare <- subset_samples(ps.ITS.nocontrols.rare, !(is.na(Competion)))

sample_data(ps.ITS.nocontrols.rare)$FunGroup <- recode(sample_data(ps.ITS.nocontrols.rare)$FunGroup, grass_x_forb = "Competition")

sample_data(ps.ITS.nocontrols.rare)$FunGroup <- factor(sample_data(ps.ITS.nocontrols.rare)$FunGroup, levels = c("Forb", "Grass", "Competition"))

GL_pcoa <- ordinate(
  physeq = ps.ITS.nocontrols.rare, 
  method = "PCoA", 
  distance = "bray")

DistBC <- phyloseq::distance(ps.ITS.nocontrols.rare, method = "bray", type="samples")

adonis(DistBC ~ FunGroup, as(sample_data(ps.ITS.nocontrols.rare), "data.frame"), permutations = 9999) # marginal fungroup diffs

# because all these are marginal, we should not be looking at pairwise comparisons

sample_data(ps.ITS.nocontrols.rare)$FunGroup <- as.factor(sample_data(ps.ITS.nocontrols.rare)$FunGroup)
pair.its.com <- adonis.pair(DistBC, as(sample_data(ps.ITS.nocontrols.rare), "data.frame")$FunGroup, nper = 9999, corr.method = "BH")
kable(pair.its.com)

#### .---Fig: Ord (16S) ####
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
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3"))

ggsave("Figures/Final/Fig-1.jpg", ord.plot, width = 5, height = 4, units = "in", dpi = 600)


#### .---Fig: Ord (ITS) ####
ord.plot.ITS <- plot_ordination(ps.ITS.nocontrols.rare, GL_pcoa, color = "FunGroup", shape = "FunGroup") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = FunGroup)) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.margin = margin(c(.1,.1,.1,.1))
  ) +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3"))

ggsave("Figures/Final/Fig-S1.jpg", ord.plot.ITS, width = 6, height = 4, units = "in", dpi = 600)

#### Shannon Diveristy (16S) ####

GL_Alpha <- estimate_richness(ps.16s.nocontrols.rare, measures = "Shannon")
GL_Alpha <- cbind(GL_Alpha, sample_data(ps.16s.nocontrols.rare))

kruskal_test(Shannon ~ FunGroup, distribution = approximate(nresample = 9999), data = GL_Alpha)

dunnTest(Shannon ~ FunGroup, data = GL_Alpha, method = "bh") 

#### Shannon Diversity (ITS) ####

GL_Alpha_ITS <- estimate_richness(ps.ITS.nocontrols.rare, measures = "Shannon")
GL_Alpha_ITS <- cbind(GL_Alpha_ITS, sample_data(ps.ITS.nocontrols.rare))

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
  
ggsave("Figures/Final/Fig-2.jpg", p.rich.16s, width = 5, height = 4, units = "in", dpi = 600)

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

ggsave("Figures/Final/Fig-S2.jpg", p.rich.ITS, width = 7, height = 4, units = "in", dpi = 600)

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
  fg.cols <- ifelse(res.alpha.tax$Family %in% FG.G, "red", ifelse(res.alpha.tax$Family %in% FG.F, "blue", "black"))
  
   p <- ggplot(res.alpha.tax, aes(x = Family, y = log2FoldChange)) + 
    geom_point(size = 3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 12, colour = fg.cols),
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
  #fg.cols <- ifelse(res.alpha.tax$Family %in% FG.G, "red", ifelse(res.alpha.tax$Family %in% FG.F, "blue", "black")) 
  #fg.cols <- c("blue", "black", "blue", "black", "black", "blue", "black", "black", "black", "black","black", "black","black", "blue")
  #fg.cols <- c("blue", "blue", "red", "black", "black", "blue", "red", "black", "blue")
  
   p <- ggplot(res.alpha.tax, aes(x = Family, y = log2FoldChange)) + 
    geom_point(size = 3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 12, colour = fg.cols),
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

ggsave("Figures/Final/Fig-S3a.jpg", d.plot, width = 11, height = 4, dpi = 600)

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

# test <- as.data.frame(t(assays(dds.its.treat)[["counts"]]))
# test2 <- assays(dds.its.treat)[["replaceCounts"]])

#plot results

des.p.its <- plot.list$G.F + plot.list$FG.F + plot.FG.G

ggsave("Figures/Final/Fig-S3b.jpg", des.p.its, width = 9, height = 4, dpi = 600)

###
# Extract Normalized abundances
###
rld <- rlog(dds.16s.treat, blind = F)
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
fam2 <- c("Burkholderiaceae", "Methylophilaceae", "Fibrobacteraceae", "Veillonellaceae", "Clostridiaceae_1", "Rhodocyclaceae")
fig.norm <- merge(norm.count.16s, mapping[,c(1,35)], by.x = "Sample", by.y = "SampleID_Fix")
fig.norm <- filter(fig.norm, Family %in% fam2)
fig.norm$FunGroup <- recode(fig.norm$FunGroup, grass_x_forb = "Competition")
#fig.norm.sum <- summarySE(fig.norm, groupvars = c("Family", "FunGroup"), measurevar = "norm.counts")
fig.norm.sum <- ddply(fig.norm, .(Family, FunGroup), summarize, norm.counts = max(norm.counts))

fig.norm.sum$cld <- c("a", "b", "a",
                      "a", "b", "b",
                      "a", "b", "b",
                      "a", "b", "b",
                      "a", "b", "a",
                      "a", "b", "b")
plotlist = list()

for(i in fam2) {
  p <- ggplot(fig.norm[fig.norm$Family == i,], aes(x = FunGroup, y = norm.counts, col = FunGroup, col = FunGroup)) +
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
  p <- ggplot(fig.norm[fig.norm$Family == i,], aes(x = FunGroup, y = norm.counts, col = FunGroup, col = FunGroup)) +
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

ggsave("Figures/Final/Fig-S4.jpg", p, width = 7, height = 9, units = "in", dpi = 600)

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
    ylim(NA, max(fig.norm[fig.norm$Family == i,]$norm.counts) + 0.1*max(fig.norm[fig.norm$Family == i,]$norm.counts))
  
  plotlist[[i]] <- p
}


p <- grid.arrange(arrangeGrob(grobs = plotlist, ncol = 3, left = textGrob("Normalized Abundance", rot = 90, vjust = 1)))

ggsave("Figures/Final/Fig-S5.jpg", p, width = 8, height = 2.5, units = "in", dpi = 600)


#### SourceTracker (16S) ####
asv.fam <- data.frame(ASV = row.names(df.16s.tax), Family = df.16s.tax$Family)
asv.fam <- filter(asv.fam, asv.fam$Family %in% fam2)

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
ST_16S2$rel.prob <- ST_16S2$source.prob/ST_16S$sum.prob

####.---Fig: SourceTracker (16S - M) ####
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
asv.fam <- filter(asv.fam, asv.fam$Family %in% fam3)

ST_16S3 <- merge(ST_16S, asv.fam, by = "ASV", all.x = F)
test <- ddply(ST_16S3, .(Family), summarize, sum.prob = sum(prob))
ST_16S3 <- ddply(ST_16S3, .(Family, Source), summarize, source.prob = sum(prob))
ST_16S3 <- merge(ST_16S3, test, by = "Family")
ST_16S3$rel.prob <- ST_16S3$source.prob/ST_16S3$sum.prob

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
  
ggsave("Figures/Final/Fig-S6.jpg", p.st, width = 6, height = 3, dpi = 600, units = "in")

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
  
ggsave("Figures/Final/Fig-S7.jpg", p.st, width = 5, height = 1.5, dpi = 600, units = "in")

#### Log Response Ratio ####
#mapping$FunGroup <- as.factor(mapping$FunGroup)

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

ggsave("Figures/Final/Fig-5.jpg", bio.plot, width = 5, height = 4, units = "in", dpi = 600)

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

#### Models: LRR v NA (16S) ####

group <- c("grass", "forb")
fam.16s <- c(fam2, fam3)
fam.16s <- as.character(fam.16s)
lrr.16s.m <- expand.grid(Family = fam.16s, FunGroup2 = group, est = NA, se = NA, p = NA)


for(i in fam.16s) {
  for(j in group){
  tmp <- lmer(weight.d ~ norm.counts + (1|Treatment), data = lrr.fam.16s[lrr.fam.16s$Family == i & lrr.fam.16s$FunGroup2 == j,])
  lrr.16s.m[lrr.16s.m$Family == i & lrr.16s.m$FunGroup2 == j, "est"] <- summary(tmp)[["coefficients"]][2,1]
    lrr.16s.m[lrr.16s.m$Family == i & lrr.16s.m$FunGroup2 == j, "se"] <- summary(tmp)[["coefficients"]][2,2]
      lrr.16s.m[lrr.16s.m$Family == i & lrr.16s.m$FunGroup2 == j, "p"] <- summary(tmp)[["coefficients"]][2,5]
  }
}

lrr.16s.m$lines <- ifelse(lrr.16s.m$p < 0.05, "dashed", "solid") # this should give me the exact opposite but for some reason it's switched in the graph... 

#### Models: LRR v NA (ITS) ####

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

#### .---Fig: LRR v NA (16S - M) ####

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
                               ncol = 3, 
                               left = textGrob("Log Response Ratio", 
                                               rot = 90, 
                                               vjust = 1)),  
                   legend, 
                   widths = unit.c(unit(1, "npc") - legend$width, legend$width), 
                   nrow=1)

ggsave("Figures/Final/Fig-6.jpg", p2, width = 7.5, height = 5, units = "in", dpi = 600)

#### .---Fig: LRR v NA (16S - S) ####
plotlist = list()

for(i in fam3) {
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

ggsave("Figures/Final/Fig-S8.jpg", p2, width = 7, height = 9, units = "in", dpi = 600)


#### .---Fig: LRR v NA (ITS) ####
plotlist = list()

for(i in fam.ITS) {
  p <- ggplot(lrr.fam.ITS[lrr.fam.ITS$Family == i,], aes(x = norm.counts, y = weight.d, col = FunGroup2, group = FunGroup2)) + 
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

ggsave("Figures/Final/Fig-S9.jpg", p2, width = 7, height = 2, units = "in", dpi = 600)

