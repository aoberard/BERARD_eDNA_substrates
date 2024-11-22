# Script set-up ----

## Library and functions ----
library(readxl)
library(dplyr)
library(stringr)

#Manually corrected abundance files
raw_12s_abundance <- readxl::read_excel(here::here("data/raw-data/2024.10.30_12SV5_eDNA_abundance_samplecorrected.xlsx"))
raw_16s_abundance <- readxl::read_excel(here::here("data/raw-data/2024.11.10_16Smam_eDNA_abundance_samplecorrected.xlsx"))

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
  mutate(primer = "16Smam")


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
  summarize(distinct_classes = n_distinct(Class),
            distinct_orders = n_distinct(Order),
            distinct_families = n_distinct(Family),
            distinct_genera = n_distinct(Genus),
            distinct_species = n_distinct(Species)) %>%
  filter(distinct_classes > 1 | distinct_orders > 1 | 
           distinct_families > 1 | distinct_genera > 1 | 
           distinct_species > 1)



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
            sum_reads = sum(reads)) %>%
  ungroup()


## Pooling the two primers data ----

### Mix affiliation from both primers ----

# Conversion of final_affiliation and related taxonomy for shared taxa between primers
edna_ppooled_taxmixed <- edna_ppooled %>%
  mutate(final_affiliation = case_when( 
    final_affiliation == "Apodemus" ~ "Apodemus_sylvaticus",
    final_affiliation == "Aegithalos" ~ "Aegithalos_caudatus",
    final_affiliation == "Phylloscopus" ~ "Phylloscopus_collybita",
    final_affiliation == "Sturnus" ~ "Sturnus_vulgaris",
    final_affiliation == "Bufonidae" ~ "Bufo",
    final_affiliation %in% c("Turdus_philomelos","Turdus_merula") ~ "Turdus",
    final_affiliation %in% c("Corvus","Pica_pica") ~ "Corvidae",
    final_affiliation %in% c("Columba_livia", "Columba_palumbus") ~ "Columbidae",
    final_affiliation == "Passeriformes" ~ "Passeriformes_3",
    TRUE ~ final_affiliation
  )) %>%
  mutate(Species = case_when( 
    final_affiliation == "Apodemus_sylvaticus" ~ "Apodemus_sylvaticus",
    final_affiliation == "Aegithalos_caudatus" ~ "Aegithalos_caudatus",
    final_affiliation == "Phylloscopus_collybita" ~ "Phylloscopus_collybita",
    final_affiliation == "Sturnus_vulgaris" ~ "Sturnus_vulgaris",
    final_affiliation == "Turdus" ~ "Multi-affiliation",
    final_affiliation == "Corvidae" ~ "Multi-affiliation",
    final_affiliation == "Columbidae" ~ "Multi-affiliation",
    TRUE ~ Species
  )) %>%
  mutate(Genus = case_when(
    final_affiliation == "Bufo" ~ "Bufo",
    final_affiliation == "Corvidae" ~ "Multi-affiliation",
    final_affiliation == "Columbidae" ~  "Multi-affiliation",
    TRUE ~ Genus
  ))

#Regroup rows of shared taxa for each primers
edna_ppooled_taxmixed <- edna_ppooled_taxmixed %>%
  group_by(sample, primer, final_affiliation) %>%
  summarise(across(c(Class, Order, Family, Genus, Species), first),
          sum_positive_replicate = max(sum_positive_replicate),      # choice of keeping only highest number of positive replicates
          pooled_number = first(pooled_number),
          sum_reads = sum(sum_reads)) %>%
  ungroup()

#Pool primers and generate columns to count positive replicates per primers
edna_gpooled <- edna_ppooled_taxmixed %>%
  mutate( sum_positive_replicate_12s = case_when( primer == "12SV5" ~ sum_positive_replicate, TRUE ~ 0)) %>%
  mutate( sum_positive_replicate_16s = case_when( primer == "16Smam" ~ sum_positive_replicate, TRUE ~ 0)) %>%
  group_by(sample, final_affiliation) %>%
  summarize(across(c(Class, Order, Family, Genus, Species), first),
            sum_positive_replicate = sum(sum_positive_replicate),
            sum_positive_replicate_12s = sum(sum_positive_replicate_12s),
            sum_positive_replicate_16s = sum(sum_positive_replicate_16s),
            pooled_number = sum(pooled_number),
            pooled_percent_positive = sum_positive_replicate / pooled_number,
            sum_reads = sum(sum_reads))

#Check pooling
unique(edna_gpooled$pooled_number)  #should equal 6 or 3
unique(edna_gpooled$sum_positive_replicate) #should be above 6


# Filter data ----

## Filter choice ▲ ----
filter_zoo <- TRUE
filter_conta <- TRUE
filter_domestic <- TRUE
filter_affiliation_level <- TRUE

## Variable creation for data filtering ----
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
#Phasianus_colchicus, Numida_meleagris and Alectoris_rufa can be bred by humans
#but are not considered to be domestic animals here, as their presence here probably comes from free-ranging individuals.

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
    stringr::str_detect(final_affiliation, Species) ~ "species",
    stringr::str_detect(final_affiliation, Genus) ~ "genus",
    stringr::str_detect(final_affiliation, Family) ~ "family",
    stringr::str_detect(final_affiliation, Order) ~ "order",
    stringr::str_detect(final_affiliation, Class) ~ "class",
    TRUE ~ NA 
  ))

edna_ppooled <- edna_ppooled %>%
  mutate(affiliation_level = case_when(
    stringr::str_detect(final_affiliation, Species) ~ "species",
    stringr::str_detect(final_affiliation, Genus) ~ "genus",
    stringr::str_detect(final_affiliation, Family) ~ "family",
    stringr::str_detect(final_affiliation, Order) ~ "order",
    stringr::str_detect(final_affiliation, Class) ~ "class",
    TRUE ~ NA
  ))

## Apply filters based on selected criteria ----
edna_pfiltered <- edna_ppooled
edna_gfiltered <- edna_gpooled

if (filter_zoo) {
  edna_pfiltered <- edna_pfiltered %>%
    filter(in_zoo == FALSE)
  
  edna_gfiltered <- edna_gfiltered %>%
    filter(in_zoo == FALSE)
}

if (filter_conta) {
  edna_pfiltered <- edna_pfiltered %>%
    filter(!final_affiliation %in% c("Homo_sapiens", "Sardina_pilchardus"))
  
  edna_gfiltered <- edna_gfiltered %>%
    filter(!final_affiliation %in% c("Homo_sapiens", "Sardina_pilchardus"))
}

if (filter_domestic) {
  edna_pfiltered <- edna_pfiltered %>%
    filter(domestic == FALSE)
  
  edna_gfiltered <- edna_gfiltered %>%
    filter(domestic == FALSE)
}

if (filter_affiliation_level) {
  edna_pfiltered <- edna_pfiltered %>%
    filter(affiliation_level %in% c("species", "genus"))
  
  edna_gfiltered <- edna_gfiltered %>%
    filter(affiliation_level %in% c("species", "genus"))
}


# Writing data file ----
edna_gpooled %>%
  write.csv(file = here::here("data","derived-data","edna_gpooled.csv"))

edna_ppooled %>%
  write.csv(file = here::here("data","derived-data","edna_ppooled.csv"))

edna_pfiltered %>%
  write.csv(file = here::here("data","derived-data","edna_pfiltered.csv"))

edna_gfiltered %>%
  write.csv(file = here::here("data","derived-data","edna_gfiltered.csv"))

