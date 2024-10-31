# Script set-up ----

## Library and functions ----
library(readxl)
library(dplyr)
library(stringr)

## Import raw data ▲ attention la data sample replica corrected ajout colonnes vides----

#Manually corrected abundance files
raw_12s_abundance <- readxl::read_excel(here::here("data/raw-data/2024.10.30_12SV5_eDNA_abundance_samplecorrected.xlsx"))
raw_16s_abundance <- readxl::read_excel(here::here("data/raw-data/2024.10.30_16Smam_eDNA_abundance_samplecorrected.xlsx"))

#Collection data
odk <- readr::read_csv(here::here("data","raw-data","2024.10.29_odk_eDNA.csv"))
odk_w <- readr::read_csv(here::here("data","raw-data","2024.10.29_odk_eDNA_web.csv"))

#Join collection files
odk <- odk %>% filter(type__prelevement != "ecouvillon_piege")
odk <- left_join(odk, odk_w, by = c("KEY" = "PARENT_KEY"))
odk <- odk %>% select(`info-localite`, id_tube_toile, `ecouvillon_feuille-id_tube_feuille`, `sol_ligne-id_tube_sol_ligne`)
rm(odk_w)


## Global data creation ----

### Pivot longer data ----

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


### Bind data from the two primers ----
#Keep shared columns only
pivot_12s <- pivot_12s %>%
  select(intersect(colnames(pivot_12s), colnames(pivot_16s)))

#Binding rows of 12s and 16s
all_edna <- rbind(pivot_12s, pivot_16s )


### Check global data ----

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


## Mix affiliation from both primers ----

# all_edna <- all_edna %>%
#   mutate(case_when( final_affiliation == "Apodemus" ~ "Apodemus_sylvaticus",
#                     final_affiliation == "Aegithalos" ~ "Aegithalos_caudatus",
#                     final_affiliation == "Passeriformes" ~ "Passeriformes_3",
#                     final_affiliation == "Phylloscopus" ~ "Phylloscopus_collybita",
#                     final_affiliation == "Sturnus" ~ "Sturnus_vulgaris",
#                     final_affiliation == "Bufonidae" ~ "Bufo"
#   ))
# 
# # Penser que ce serait bien de savoir à quel rang est le final affiliation
# # mais peut-être l'ajouter que a la fin du script avec les autres colonnes pour filtrer
# 
# # ILF FAUT PAS FAIRE COMME J'AI COMMENCE A FAIRE, JUSTE remplacer nom affiliation en 
# # se disant que va etre pool ensuite lors du pool, car vaut artificiellement augmenter nbr replica
# # parfois, donc faut faire rassemblement affiliation avant decompte replica
# 
# 
# AFAIRE (en essayant garder info que issu de pools qqpart ?)
# attention bien chager nom apres si la creation nouveau nom (all_edna utilise apres) sinon garder meme si pas de soucis
# 
# 
# 
# ATTNETION LORS DU POOLING PAR REPLICATS ENSUITE ON FAIT UN POOL DES REPLICAS PAR UN CLUSTER, SAUF QUE CLUSTER ON LE MEME NOM MAIS PAS FORCEMENT
# ATTACHE A LA MEME AFFILIATION SELON LES RUNS, DONC VAUT MIEUX SOIT POOL PAR CLUSTER-amorces-run soit par affiliation (affiliation ce qui sera fait à la fin là
# )



## Pooling replicates per sample ----
  
#Create a positive_replicate column to count positive replicates per samples-cluster when pooling
all_edna <- all_edna %>%
  mutate(positive_replicate = case_when(reads > 0 ~ 1, 
                                        reads == 0 ~ 0))

### By pooling the two primers indiscriminately ----
#Pooling per samples-cluster
edna_gpooled <- all_edna %>%
  group_by(sample, final_affiliation) %>%
  summarize(across(c(Class, Order, Family, Genus, Species), first),
            sum_positive_replicate = sum(positive_replicate),
            pooled_number = n(),
            sum_reads = sum(reads),
            primers = "12S + 16S")


hist(edna_gpooled$sum_positive_replicate)




# Check DATA pooling : que meme nombre que expected -----

unique_sample_count <- n_distinct(all_edna$sample)
unique_affiliation_count <- n_distinct(all_edna$final_affiliation)
expected_rows <- unique_sample_count * unique_affiliation_count
actual_rows <- nrow(edna_gpooled)

# Calculate the expected distinct groups
expected_groups <- expand.grid(sample = unique(all_edna$sample), 
                               final_affiliation = unique(all_edna$final_affiliation))
# Identify missing or collapsed combinations
missing_or_combined <- anti_join(expected_groups, edna_gpooled, by = c("sample", "final_affiliation"))













### By pooling when taking into account primer (alternative) ----
#Pooling per samples-cluster-primer
edna_ppooled <- all_edna %>%
  group_by(sample, primer, final_affiliation, ) %>%
  summarize(across(c(Class, Order, Family, Genus, Species), first),
            sum_positive_replicate = sum(positive_replicate),
            pooled_number = n(),
            sum_reads = sum(reads))

hist(edna_ppooled$sum_positive_replicate)


#(maybe it will be useful at some point ?)

# eclaircir methode analyses

# peut etre commencer par se poser la question de ce qui doit etre decrit (ex nbr moyens de replicas positif par cluster
# et par x et y) - ce qui en somme n'a pas à être analysé pour la suite
  
  



## Final variable addition ----
# ces colonnes doivent permettre de faire les choix pour la suite quant à ce qui doit etre pris en compte dans les analyses
# y'a toujours la sardine, enlever avant ou ? (bien verif que y'avait pas dans les neg et qu'il faut pas l'enlever ou quoi des data)

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
