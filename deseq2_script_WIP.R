#DESEQ2 Analysis

library(DESeq2)

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
mapping <- read.csv("/Data/Grassland-Amplicon-Mapping-File4.csv") #Cassie's path

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
