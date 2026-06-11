library(here)
library(biomaRt)
library(tximport)
library(dplyr)

# reading in AGM data
AGM_meta <- readRDS(here("support", "meta"))
AGM_txi <- readRDS(here("support", "txi_rename"))

# reading in ES cell data
ESbulk_meta <- readRDS(here("support", "meta_ESbulk"))
ESbulk_txi <- readRDS(here("support", "txi_ESbulk"))

# genes list
genes_list <- readRDS(here("support", "genes_list"))

# merging AGM and ES data
merged <- list(abundance = cbind(AGM_txi$abundance, ESbulk_txi$abundance), 
               counts = cbind(AGM_txi$counts, ESbulk_txi$counts),
               length = cbind(AGM_txi$length, ESbulk_txi$length),
               countsFromAbundance = "no")
saveRDS(merged, here("output", "merged"))

# creating new coldata
# fixing a few things to allow rbind
colnames(ESbulk_meta)[1] <- "Sample.name"
ESbulk_meta <- mutate(ESbulk_meta, rep = gsub(".*_", "", Sample.name))

AGM_meta <- mutate(AGM_meta, sort = gsub("_[^_]+$", "", title)) # remove everything after the last underscore
AGM_meta$embryo <- paste0(AGM_meta$embryo, "_embryo")
colnames(AGM_meta)[5] <- "rep"
# taking relevant columns
ES <- select(ESbulk_meta, Sample.name, sort, title, rep)
AGM <- select(AGM_meta, Sample.name, sort, title, rep)
# bind
meta <- rbind(AGM, ES)
  
saveRDS(meta, here("output", "meta"))



