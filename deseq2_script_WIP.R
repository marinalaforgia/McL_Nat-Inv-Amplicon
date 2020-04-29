#DESEQ2 Analysis

library(DESeq2)
library(tidyverse)
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

its.treat = phyloseq_to_deseq2(ps.ITS.nocontrols, ~ TreatmentName)

dds.its.treat = DESeq(its.treat, test="Wald", fitType="parametric")

contrasts = c("SA.G", "SA.ST", "SA.SAG", "SA.STG", 
              "ST.G", "SAG.G", "STG.G", "ST.SAG",
              "ST.STG", "SAG.STG", "SA.BSS", "ST.BSS",
              "G.BSS", "SAG.BSS", "STG.BSS")

#I think I got them all?
contrast.list <- list(SA.G = c("TreatmentName", "Stress_avoiding_forb", "Grass"),
                      SA.ST = c("TreatmentName", "Stress_avoiding_forb", "Stress_tolerant_forb"),
                      SA.SAG = c("TreatmentName", "Stress_avoiding_forb", "SA_forb_X_grass"),
                      SA.STG = c("TreatmentName", "Stress_avoiding_forb", "ST_forb_X_grass"),
                      ST.G = c("TreatmentName", "Stress_tolerant_forb", "Grass"),
                      SAG.G = c("TreatmentName", "SA_forb_X_grass", "Grass"),
                      STG.G = c("TreatmentName", "ST_forb_X_grass", "Grass"),
                      ST.SAG = c("TreatmentName", "Stress_tolerant_forb", "SA_forb_X_grass"),
                      ST.STG = c("TreatmentName", "Stress_tolerant_forb", "ST_forb_X_grass"),
                      SAG.STG = c("TreatmentName", "SA_forb_X_grass", "ST_forb_X_grass"),
                      SA.BSS = c("TreatmentName", "Stress_avoiding_forb", "Background_soil_and_sand_mix"),
                      ST.BSS = c("TreatmentName", "Stress_tolerant_forb", "Background_soil_and_sand_mix"),
                      G.BSS = c("TreatmentName", "Grass", "Background_soil_and_sand_mix"),
                      SAG.BSS = c("TreatmentName", "SA_forb_X_grass", "Background_soil_and_sand_mix"),
                      STG.BSS = c("TreatmentName", "ST_forb_X_grass", "Background_soil_and_sand_mix"))

res.list <- vector("list")

for(i in contrasts) {
  res.list[[paste(r,i,sep = ".")]] <- broom::tidy(results(dds.its.treat, contrast = contrast.list[[i]], pAdjustMethod = "bonferroni"))
}

#tidy results
df <- plyr::ldply(res.list, function(x) x)
names(df)[1] <- "Contrast"
names(df)[2] <- "ASV"
