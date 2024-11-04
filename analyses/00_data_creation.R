# Script set-up ----

## Library and functions ----
library(readxl)
library(dplyr)
library(stringr)

#Manually corrected abundance files
raw_12s_abundance <- readxl::read_excel(here::here("data/raw-data/2024.10.30_12SV5_eDNA_abundance_samplecorrected.xlsx"))
raw_16s_abundance <- readxl::read_excel(here::here("data/raw-data/2024.11.04_16Smam_eDNA_abundance_samplecorrected.xlsx"))

#Collection data
odk <- readr::read_csv(here::here("data","raw-data","2024.10.29_odk_eDNA.csv"))
odk_w <- readr::read_csv(here::here("data","raw-data","2024.10.29_odk_eDNA_web.csv"))

#Join collection files
odk <- odk %>% filter(type__prelevement != "ecouvillon_piege")
odk <- left_join(odk, odk_w, by = c("KEY" = "PARENT_KEY"))
odk <- odk %>% select(`info-localite`, id_tube_toile, `ecouvillon_feuille-id_tube_feuille`, `sol_ligne-id_tube_sol_ligne`)
rm(odk_w)


# Global data creation ----

## Pivot longer data ----

#Pivot longer data both file
pivot_12s <- raw_12s_abundance%>%
  tidyr::pivot_longer(cols =  (which(colnames(raw_12s_abundance) == "observation_sum")+1) :ncol(raw_12s_abundance),
                      names_to = "PCRreplicate",
                      values_to = "reads")
pivot_16s <- raw_16s_abundance%>%
  tidyr::pivot_longer(cols =  (which(colnames(raw_16s_abundance) == "observation_sum")+1) :ncol(raw_16s_abundance),
                      names_to = "PCRreplicate",
                      values_to = "reads")

#Create sample column
pivot_12s <- pivot_12s %>%
  mutate(sample = stringr::str_extract(PCRreplicate, "^[^-]*")) %>%
  mutate(primer = "12SV5")
pivot_16s <- pivot_16s %>%
  mutate(sample = stringr::str_extract(PCRreplicate, "^[^-]*")) %>%
  mutate(primer = "16Smamm")


## Bind data from the two primers ----
#Keep shared columns only
pivot_12s <- pivot_12s %>%
  select(intersect(colnames(pivot_12s), colnames(pivot_16s)))

#Binding rows of 12s and 16s
all_edna <- rbind(pivot_12s, pivot_16s )


## Check global data ----

#Number of samples per primers
pivot_12s %>%
  distinct(sample) %>%
  nrow()
pivot_16s %>%
  distinct(sample) %>%
  nrow()

#Identify samples collected absents in abundance files
samples_name <- odk %>% select(`ecouvillon_feuille-id_tube_feuille`,`sol_ligne-id_tube_sol_ligne`,id_tube_toile) %>%
  tidyr::pivot_longer(cols = everything(), values_to = "valeurs") %>%
  pull(valeurs)

setdiff(pivot_12s %>% distinct(sample) %>% pull, samples_name)
setdiff(pivot_16s %>% distinct(sample) %>% pull, samples_name)
setdiff(samples_name, pivot_12s %>% distinct(sample) %>% pull)
setdiff(samples_name, pivot_16s %>% distinct(sample) %>% pull)

#Check differences of final_affiliation
pivot_12s %>%
  distinct(final_affiliation) %>%
  nrow()
pivot_16s %>%
  distinct(final_affiliation) %>%
  nrow()

#Identify samples without correct number of replicates
pivot_12s %>% group_by(sample, final_affiliation) %>%
  filter(n() < 3) %>%
  pull(sample) %>%
  unique()
pivot_12s %>% group_by(sample, final_affiliation) %>%
  filter(n() > 3) %>%
  pull(sample) %>%
  unique()
pivot_16s %>% group_by(sample, final_affiliation) %>%
  filter(n() < 3) %>%
  pull(sample) %>%
  unique()
pivot_16s %>% group_by(sample, final_affiliation) %>%
  filter(n() > 3) %>%
  pull(sample) %>%
  unique()

#Check if taxonomy is similar between final_affiliation from the two primers (should be empty)
all_edna %>%
  group_by(sample, final_affiliation) %>%
  summarize(distinct_classes = n_distinct(Class), # Count distinct values
            distinct_orders = n_distinct(Order),
            distinct_families = n_distinct(Family),
            distinct_genera = n_distinct(Genus),
            distinct_species = n_distinct(Species)) %>%
  filter(distinct_classes > 1 | distinct_orders > 1 | 
           distinct_families > 1 | distinct_genera > 1 | 
           distinct_species > 1)


# ATTENTION ATTENTION Mix affiliation from both primers ----
all_edna <- all_edna %>%
  mutate(case_when( final_affiliation == "Apodemus" ~ "Apodemus_sylvaticus",
                    final_affiliation == "Aegithalos" ~ "Aegithalos_caudatus",
                    final_affiliation == "Passeriformes" ~ "Passeriformes_3",
                    final_affiliation == "Phylloscopus" ~ "Phylloscopus_collybita",
                    final_affiliation == "Sturnus" ~ "Sturnus_vulgaris",
                    final_affiliation == "Bufonidae" ~ "Bufo",
                    final_affiliation %in% c("Columba_livia", "Columba_palumbus") ~ "Columbidae",
                    final_affiliation %in% c("Corvus","Pica_pica") ~ "Corvidae",
                    final_affiliation %in% c("Turdus_philomelos","Turdus_merula") ~ "Turdus"
  ))

# attention a probleme taxonomique aussi ensuite ! (jointure tables + 
# colonne nouvellement créer qui utilise que multiaffiliation pour identifier quel rang est la final affiliation)
# 
# # ILF FAUT PAS FAIRE COMME J'AI COMMENCE A FAIRE, JUSTE remplacer nom affiliation en
# # se disant que va etre pool ensuite lors du pool, car vaut artificiellement augmenter nbr replica
# # parfois, donc faut faire rassemblement affiliation avant decompte replica
# Dans tous les cas le pooling va faire augmentation nbr replica, car reste les mêmes individus
# 
# Ou alors faire pooling amorces après pooling de chaque amorce comme ça on peut
# pb de cluster qui apparaissent que dans un run et pas dans l'autre, l'autre tous negatifs à ça?
# 



# Pooling replicates per sample ----
  
#Create a positive_replicate column to count positive replicates per samples-cluster when pooling
all_edna <- all_edna %>%
  mutate(positive_replicate = case_when(reads > 0 ~ 1, 
                                        reads == 0 ~ 0))

## Pooling for each primiers ----
#Pooling per samples-cluster-primer
edna_ppooled <- all_edna %>%
  group_by(sample, primer, final_affiliation ) %>%
  summarize(across(c(Class, Order, Family, Genus, Species), first),
            sum_positive_replicate = sum(positive_replicate),
            pooled_number = n(),
            sum_reads = sum(reads))

hist(edna_ppooled$sum_positive_replicate)


## Pooling the two primers data ----
#Generate columns to count positive replicates per primers and then pool them per sample -affiliation
edna_gpooled <- edna_ppooled %>%
  mutate( sum_positive_replicate_12s = case_when( primer == "12SV5" ~ sum_positive_replicate, TRUE ~ 0)) %>%
  mutate( sum_positive_replicate_16s = case_when( primer == "16Smamm" ~ sum_positive_replicate, TRUE ~ 0)) %>%
  group_by(sample, final_affiliation) %>%
  summarize(across(c(Class, Order, Family, Genus, Species), first),
            sum_positive_replicate = sum(sum_positive_replicate),
            sum_positive_replicate_12s = sum(sum_positive_replicate_12s),
            sum_positive_replicate_16s = sum(sum_positive_replicate_16s),
            pooled_number = sum(pooled_number),
            pooled_percent_positive = round(100*sum_positive_replicate / pooled_number, 0),
            sum_reads = sum(sum_reads))

#Check pooling
unique(edna_gpooled$pooled_number)  #should equal 6 or 3



# Final variable creation for data filtering ----
#New columns are created in order to make choice about which rows to consider for further analysis

#Variable to identify substrate type
edna_gpooled <- edna_gpooled %>%
  mutate(substrate = case_when(
    stringr::str_detect(sample, "LEAF") ~ "leaf",
    stringr::str_detect(sample, "SOIL") ~ "soil",
    stringr::str_detect(sample, "WEB") ~ "spiderweb",
    .default = "other"
  ))
edna_ppooled <- edna_ppooled %>%
  mutate(substrate = case_when(
    stringr::str_detect(sample, "LEAF") ~ "leaf",
    stringr::str_detect(sample, "SOIL") ~ "soil",
    stringr::str_detect(sample, "WEB") ~ "spiderweb",
    .default = "other"
  ))

#Variable to identify possible domestic taxa
domestic_affiliation_names <- c("Sus_scrofa",
                                "Gallus_gallus",
                                "Bos_taurus",
                                "Capra_hircus",
                                "Ovis_aries",
                                "Canis_lupus",
                                "Equus_caballus",
                                "Meleagris_gallopavo"
                                )

#Phasianus_colchicus ? Numida_meleagris ? Alectoris_rufa? mises avec les domestiques ? (peutetre plus pintade que faisan mais jsp)
# serait surtout numida en vrai parmi eux, mais et encore
################### VERIFIER LISTE MISE DEDANS APRES AVEC MAX
# lui avait pas noter capra hircus, ni cheval ni dinde

edna_gpooled <- edna_gpooled %>%
  mutate(domestic = final_affiliation %in% domestic_affiliation_names)
edna_ppooled <- edna_ppooled %>%
  mutate(domestic = final_affiliation %in% domestic_affiliation_names)

#Variable to identify samples from the ZOO
edna_gpooled <- edna_gpooled %>%
  mutate(in_zoo = stringr::str_detect(sample, "ZOO"))
edna_ppooled <- edna_ppooled %>%
  mutate(in_zoo = stringr::str_detect(sample, "ZOO"))

#Variable to identify level of final_affiliation
edna_gpooled <- edna_gpooled %>%
  mutate(affiliation_level = case_when(
    Species != "Multi-affiliation" ~ "species",
    Genus !="Multi-affiliation" ~ "genus",
    Family !="Multi-affiliation" ~ "family",
    Order !="Multi-affiliation" ~ "order",
    Class !="Multi-affiliation" ~ "class",
  ))

edna_ppooled <- edna_ppooled %>%
  mutate(affiliation_level = case_when(
    Species != "Multi-affiliation" ~ "species",
    Genus !="Multi-affiliation" ~ "genus",
    Family !="Multi-affiliation" ~ "family",
    Order !="Multi-affiliation" ~ "order",
    Class !="Multi-affiliation" ~ "class",
  ))

