#### WIP deseq2 script ####
# Will merge when complete & running

#DESEQ2 Analysis

library(DESeq2)
library(tidyverse)
library(biobroom)
library(patchwork)
packageVersion("DESeq2") 


### Read in Data ###

## Remove background soil for analysis ##

## 16S ##
ps.16s.nocontrols <- readRDS("Data/") # 16s data
ps.16s.nocontrols <- readRDS("Data/") #16s data, Cassie's path


## ITS ##
ps.ITS.nocontrols <- readRDS("Data/ITS/greenhouse_its_itsx_decontam_controlsremoved.rds") # ITS data
ps.ITS.nocontrols <- readRDS("Data/greenhouse_its_itsx_decontam_controlsremoved.rds") # ITS data, Cassie's path


# Metadata mapping file including biomass and traits (version4)
mapping <- read.csv("/Users/Marina/Documents/UC-Davis/Projects/McL_Nat-Inv-Amplicon/Data/Setup/Grassland-Amplicon-Mapping-File4.csv") 
mapping <- read.csv("Data/Grassland-Amplicon-Mapping-File4.csv") #Cassie's path

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
sample_data(ps.ITS.nocontrols) <- mapping[which(mapping$SampleID_Fix %in% ps.ITS.nocontrols.rare@sam_data$SampleID_Fix),]
sample_data(ps.ITS.nocontrols) # fixed!

sample_data(ps.16s.nocontrols) <- mapping

rm(missing, missing2, df.16s, df.ITS, i)


#### Differential abundance analysis ####

# Write scripts to test for ASVs that differ with:
# Plant species
# Competition
# These should mirror ordination comparisons, etc


# Can be done at ASV level or at higher levels of taxonomic classification (e.g. family, order)
# Starting with only the ASV level

## ITS ##


# Fix taxonomy so that there are no NA's & remove headers (e.g. p__)
# This is important if we want to bind taxonomy to deseq results later

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


## Background soil vs. all other treatments ##
its.treat = phyloseq_to_deseq2(ps.ITS.nocontrols, ~ TreatmentName)

dds.its.treat = DESeq(its.treat, test="Wald", fitType="parametric")

resultsNames(dds.its.treat) #gives us comparisons for our contrasts

contrasts = c("SA.BSS", "ST.BSS",
              "G.BSS", "SAG.BSS", "STG.BSS")

#I think I got them all?
contrast.list <- list(SA.BSS = c("TreatmentName", "Stress_avoiding_forb", "Background_soil_and_sand_mix"),
                      ST.BSS = c("TreatmentName", "Stress_tolerant_forb", "Background_soil_and_sand_mix"),
                      G.BSS = c("TreatmentName", "Invasive_grass", "Background_soil_and_sand_mix"),
                      SAG.BSS = c("TreatmentName", "SA_forb_X_grass", "Background_soil_and_sand_mix"),
                      STG.BSS = c("TreatmentName", "ST_forb_X_grass", "Background_soil_and_sand_mix"))


plot.name.list <- list(SA.BSS = "SA forb vs. Background soil",
                       ST.BSS = "ST forb vs. Background soil",
                       G.BSS = "Invasive grass vs. Background soil",
                       SAG.BSS = "SA forb X grass vs. Background soil",
                       STG.BSS = "ST forb X grass vs. Background soil")

alpha = 0.01
res.list <- vector("list")
plot.list = list()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.its.treat, contrast = contrast.list[[i]], pAdjustMethod = "bonferroni")
  #filter results by p-value
  res.alpha <- res[which(res$padj < alpha), ]
  #Bind taxonomy to results
  res.alpha.tax = cbind(as(res.alpha, "data.frame"), as(tax_table(ps.ITS.nocontrols)[rownames(res.alpha), ], "matrix"))
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
  plot_list[[i]] = p
  }

#plot results

plot_list$SA.BSS + plot_list$ST.BSS + plot_list$G.BSS + plot_list$SAG.BSS + plot_list$STG.BSS

#tidy results into a table to save
df.res <- plyr::ldply(res.list, function(x) x)
names(df.res)[1] <- "Contrast"
names(df.res)[2] <- "ASV"

write.csv(df.res, "its.ddseq.treatment.csv")





# Now choosing a "new" control group to compare to
#sample_data(ps.ITS.nocontrols)$TreatmentName <- relevel(sample_data(ps.ITS.nocontrols)$TreatmentName, "Grass")

# All possible contrasts for TreamentName
# contrasts = c("SA.G", "SA.ST", "SA.SAG", "SA.STG", 
#               "ST.G", "SAG.G", "STG.G", "ST.SAG",
#               "ST.STG", "SAG.STG", "SA.BSS", "ST.BSS",
#               "G.BSS", "SAG.BSS", "STG.BSS")
# 
# 
# contrast.list <- list(SA.G = c("TreatmentName", "Stress_avoiding_forb", "Invasive_grass"),
#                       SA.ST = c("TreatmentName", "Stress_avoiding_forb", "Stress_tolerant_forb"),
#                       SA.SAG = c("TreatmentName", "Stress_avoiding_forb", "SA_forb_X_grass"),
#                       SA.STG = c("TreatmentName", "Stress_avoiding_forb", "ST_forb_X_grass"),
#                       ST.G = c("TreatmentName", "Stress_tolerant_forb", "Invasive_grass"),
#                       SAG.G = c("TreatmentName", "SA_forb_X_grass", "Invasive_grass"),
#                       STG.G = c("TreatmentName", "ST_forb_X_grass", "Invasive_grass"),
#                       ST.SAG = c("TreatmentName", "Stress_tolerant_forb", "SA_forb_X_grass"),
#                       ST.STG = c("TreatmentName", "Stress_tolerant_forb", "ST_forb_X_grass"),
#                       SAG.STG = c("TreatmentName", "SA_forb_X_grass", "ST_forb_X_grass"),
#                       SA.BSS = c("TreatmentName", "Stress_avoiding_forb", "Background_soil_and_sand_mix"),
#                       ST.BSS = c("TreatmentName", "Stress_tolerant_forb", "Background_soil_and_sand_mix"),
#                       G.BSS = c("TreatmentName", "Invasive_grass", "Background_soil_and_sand_mix"),
#                       SAG.BSS = c("TreatmentName", "SA_forb_X_grass", "Background_soil_and_sand_mix"),
#                       STG.BSS = c("TreatmentName", "ST_forb_X_grass", "Background_soil_and_sand_mix"))
