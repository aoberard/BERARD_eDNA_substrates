# Script set-up
library(readxl)
library(dplyr)

# Import raw data
raw_12s_abundance <- readxl::read_excel(here::here("data/raw-data/2024.10.24_12SV5_eDNA_abundance.xlsx"))

raw_16s_abundance <- readxl::read_excel(here::here("data/raw-data/2024.10.24_16Smam_eDNA_abundance.xlsx"))

# Pivot data
pivot_12s <- raw_12s_abundance%>%
  tidyr::pivot_longer(cols =  (which(colnames(raw_12s_abundance) == "observation_sum")+1) :ncol(raw_12s_abundance),
                      names_to = "PCRreplicate",
                      values_to = "reads")

pivot_16s <- raw_16s_abundance%>%
  tidyr::pivot_longer(cols =  (which(colnames(raw_16s_abundance) == "observation_sum")+1) :ncol(raw_16s_abundance),
                      names_to = "PCRreplicate",
                      values_to = "reads")

# Create sample column
pivot_12s <- pivot_12s %>%
  mutate(sample = stringr::str_extract(PCRreplicate, "^[^-]*")) %>%
  mutate(primer = "12SV5")
  
pivot_16s <- pivot_16s %>%
  mutate(sample = stringr::str_extract(PCRreplicate, "^[^-]*")) %>%
  mutate(primer = "16Smamm")

# Pool 12s and 16s by binding rows
pivot_12s <- pivot_12s %>%
  select(intersect(colnames(pivot_12s), colnames(pivot_16s)))

all_edna <- rbind(pivot_12s, pivot_16s )



# Là essayer de combiner au mieux les taxonomy


COMMENT FAIRE? : regrouper des clusters? transformer taxo si deja positif truc a autre truc meme groupe ?



# pour apres creer new colonnes : mutate(in_zoo = stringr::str_detect(sample, "ZOO")
# creer new colonne susbtrates (type) et pk pas type taxon ou autre


# Pooling replicates per sample

rpooled_12s <- pivote_12s %>%
  mutate(positive_replicate = case_when(reads > 0 ~ 1, 
                                        reads == 0 ~ 0))
# For each amorce-samples-cluster sum  sum the replicates ligne  
  
  



