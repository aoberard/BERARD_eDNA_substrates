# Script set-up ----

## Library and functions ----
library(purrr)
library(combinat)
library(ggplot2)
library(dplyr)
library(patchwork)


## Graphical parameters ----
palette_substrate <- c(
  "spiderweb" = "#c9c9c9",  
  "leaf" = "#c2d3a6",      
  "soil" = "#e2b19d"       
)

palette_primers <- c("12SV5" = "#FFB3B3AA", "16Smam" = "#A2D1D1AA")

# Ordering elements
class_order <- c("Aves", "Mammalia", "Amphibia", "Lepidosauria")
substrate_order <- c("soil", "leaf", "spiderweb")

# Data quality exploration ----

## Considering reads number ----
#We are seeking for the number of replicates and samples per primer and substrate under possible quality thresholds

#Choose data on which to check quality (global data with small modifications)
data_quality_check <- all_edna  # %>%
  # filter(!stringr::str_detect(sample, "ZOO")) 
  # %>% filter(final_affiliation != "Homo_sapiens")       # possible filtering for human (ubiquitous contamination)

#Calculate possible reads threshold values
threshold_values <- seq(0,
                        data_quality_check %>%
                          group_by(PCRreplicate, primer) %>%
                          summarise(total_reads = sum(reads),
                                    .groups = "drop" ) %>%
                          pull(total_reads) %>%
                          max(),
                        by = 500)

quality_variation <- data_quality_check %>%
    mutate(
      substrate = case_when(
        stringr::str_detect(PCRreplicate, "LEAF") ~ "leaf",
        stringr::str_detect(PCRreplicate, "SOIL") ~ "soil",
        stringr::str_detect(PCRreplicate, "WEB") ~ "spiderweb"
      )
    ) %>%
    group_by(PCRreplicate, primer) %>%
    summarise(
      total_reads = sum(reads),
      substrate = first(substrate),
      sample = first(sample),
      .groups = "drop"
    ) %>%
    ungroup() 

#For each threshold_values, calculate replicate and sample counts under it
quality_variation_replicate <- data.frame()
quality_variation_sample <- data.frame()

for (threshold in threshold_values){
  replicate_data <- quality_variation %>%
    filter(total_reads < threshold) %>%
    group_by(primer, substrate) %>%
    summarise(
      replicate_count = n(),
      .groups = "drop"
    ) %>%
    tidyr::complete(primer, substrate, fill = list(replicate_count = 0)) %>%
    mutate(reads_quality_threshold = threshold)
  
  quality_variation_replicate <- bind_rows(quality_variation_replicate, replicate_data)
  
  sample_data <- quality_variation %>%
    group_by(sample, primer, substrate) %>%
    summarise(
      total_reads = sum(total_reads),
      .groups = "drop"
    ) %>%
    filter(total_reads < threshold) %>% 
    group_by(primer, substrate) %>%
    summarise(
      sample_count = n(),
      .groups = "drop"
    ) %>%
    tidyr::complete(primer, substrate, fill = list(sample_count = 0)) %>%
    mutate(reads_quality_threshold = threshold)
  
  quality_variation_sample <- bind_rows(quality_variation_sample, sample_data)
  
}

quality_variation_final <- left_join(quality_variation_replicate, quality_variation_sample)
quality_variation_final <- quality_variation_final %>%
  mutate(sample_count = if_else(is.na(sample_count), 0, sample_count))


#Plot replicate count against reads threshold
ggplot(quality_variation_final, aes(x = reads_quality_threshold, y = replicate_count, color = substrate)) +
  geom_line(linewidth = 1.1) +  
  facet_wrap(~ primer) +
  scale_color_manual(values = palette_substrate) + 
  labs(
    title = "Variation of replicates quality",
    x = "Reads number threshold",
    y = "Replicates under threshold"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

#Plot sample count against reads threshold
ggplot(quality_variation_final, aes(x = reads_quality_threshold, y = sample_count, color = substrate)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ primer) +
  scale_color_manual(values = palette_substrate) +  
  labs(
    title = "Variation of samples quality",
    x = "Reads number threshold",
    y = "Samples under threshold"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

#Number of replicates and samples per primer and substrate under a chosen quality thresholds 
chosen_quality_threshold <- 1000
quality_variation_final %>%
  filter(reads_quality_threshold == chosen_quality_threshold)


## Considering repeatability ----

### Raw repeatability ----

#For data by primers (edna_pfiltered or edna_ppooled), we plot the number of positive replicates for each substrates and affiliation level 
edna_pfiltered %>%
  group_by(affiliation_level, primer, substrate, sum_positive_replicate) %>%
  summarize(count = n(),
            .groups = 'drop') %>%
  filter(sum_positive_replicate > 0) %>%
  mutate(substrate = factor(substrate, levels = c("soil", "leaf", "spiderweb"))) %>%
  tidyr::complete(affiliation_level, primer, substrate, sum_positive_replicate, fill = list(count = 0)) %>%
  ggplot(aes(x = sum_positive_replicate, y = count, fill = substrate)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(affiliation_level ~ primer) +
  scale_fill_manual(values = palette_substrate) +
  labs(
    title = "Distribution of positive Replicates by affiliation, primer, and substrate",
    x = "Number of positive replicates",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") 

#For data with pooled primers (edna_gfiltered or edna_gpooled), we plot the percent of positives replicates
#Because there are not an equal number of replicates for each affiliation, due to shared and unshared affiliation between primers
edna_gfiltered %>%
  group_by(affiliation_level, substrate, pooled_percent_positive) %>%
  summarize(count = n(),
            .groups = 'drop') %>%
  filter(pooled_percent_positive > 0) %>%
  mutate(substrate = factor(substrate, levels = c("soil", "leaf", "spiderweb"))) %>%
  tidyr::complete(affiliation_level, substrate, pooled_percent_positive, fill = list(count = 0)) %>%
  ggplot(aes(x = pooled_percent_positive, y = count, fill = substrate)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(rows = vars(affiliation_level)) +
  scale_fill_manual(values = palette_substrate) +
  labs(
    title = "Distribution of positive replicates by affiliation and substrate - with pooled primers",
    x = "Percentage of Positive Replicates",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") 


### Normalized ----

# For each primers
edna_pfiltered_normalized <- edna_pfiltered %>%
  group_by(affiliation_level, primer, substrate, sum_positive_replicate) %>%
  summarize(count = n(), .groups = 'drop') %>%
  filter(sum_positive_replicate > 0) %>%
  mutate(substrate = factor(substrate, levels = c("soil", "leaf", "spiderweb"))) %>%
  tidyr::complete(affiliation_level, primer, substrate, sum_positive_replicate, fill = list(count = 0)) 

#Normalize the count by the total per substrate
edna_pfiltered_normalized <- edna_pfiltered_normalized %>%
  group_by(substrate, primer) %>%
  mutate(
    total_count_substrate = sum(count),  
    normalized_count = count / total_count_substrate  
  ) %>%
  ungroup()

#Plot the normalized count
edna_pfiltered_normalized %>%
  ggplot(aes(x = sum_positive_replicate, y = normalized_count, fill = substrate)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(affiliation_level ~ primer) +
  scale_fill_manual(values = palette_substrate) +
  labs(
    title = "Distribution of positive replicates by affiliation, primer, and substrate - with pooled primers",
    x = "Number of positive replicates",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  ylim(0,1)


#For pooled primers

#Total count per substrate
edna_gfiltered_normalized <- edna_gfiltered %>%
  group_by(affiliation_level, substrate, pooled_percent_positive) %>%
  summarize(count = n(), 
            .groups = 'drop') %>%
  filter(pooled_percent_positive > 0) %>%
  mutate(substrate = factor(substrate, levels = c("soil", "leaf", "spiderweb"))) %>%
  tidyr::complete(affiliation_level, substrate, pooled_percent_positive, fill = list(count = 0))

#Normalize the count by the total per substrate
edna_gfiltered_normalized <- edna_gfiltered_normalized %>%
  group_by(substrate) %>%
  mutate(
    total_count_substrate = sum(count),  
    normalized_count = count / total_count_substrate  
  ) %>%
  ungroup()

#Plot the normalized count
edna_gfiltered_normalized %>%
  ggplot(aes(x = pooled_percent_positive, y = normalized_count, fill = substrate)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(rows = vars(affiliation_level)) +
  scale_fill_manual(values = palette_substrate) +
  labs(
    title = "Normalized distribution of positive replicates by affiliation and substrate",
    x = "Percentage of Positive Replicates",
    y = "Normalized Count"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  ylim(0,1)


# Final taxonomic affiliations ----

## Descriptive informations ----

# Affiliations in zoo, out of zoo and unique to zoo
affiliation_in_zoo <- edna_ppooled %>%
  filter(in_zoo == TRUE) %>%
  filter(sum_positive_replicate > 0) %>%
  select(final_affiliation) %>%
  pull() %>%
  unique()

affiliation_out_zoo <- edna_ppooled %>%
  filter(in_zoo == FALSE) %>%
  filter(sum_positive_replicate > 0) %>%
  select(final_affiliation) %>%
  pull() %>%
  unique()

setdiff(affiliation_in_zoo, affiliation_out_zoo)

#Percent of positive samples replicated by more than half of technical replicates - per primer unfiltered data
edna_ppooled %>%
  filter(sum_reads >0) %>%
  group_by(primer) %>%
  summarise(percent_within_repeated = 100 *sum(within_repeated_positive) / n(),
            effectif = n() )

#Percent of positive samples replicated by more than half of technical replicates - per substrate filtered and primer pooled data
edna_gfiltered %>%
  filter(sum_reads >0) %>%
  group_by(substrate) %>%
  summarise(percent_within_repeated = 100 *sum(within_repeated_positive) / n(),
            effectif = n() )

#Total number of distinct affiliation per primers
edna_ppooled %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(primer) %>%
  summarise(distinct_affiliation =  n_distinct(final_affiliation))

edna_pfiltered %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(primer) %>%
  summarise(distinct_affiliation =  n_distinct(final_affiliation))

#Total number of distinct affiliation  per substrate - primers
edna_ppooled %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(substrate, primer) %>%
  summarise(distinct_affiliation =  n_distinct(final_affiliation)) 

#After filters total number of distinct affiliation  per substrate - primers
edna_pfiltered %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(primer, substrate) %>%
  summarise(distinct_affiliation =  n_distinct(final_affiliation))

#After filters and primer pooling total number of distinct affiliation  per substrate
edna_gfiltered %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(substrate) %>%
  summarise(distinct_affiliation =  n_distinct(final_affiliation))

#After filters number of distinct affiliation per taxa Class and per substrates
edna_pfiltered %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(substrate, Class) %>%
  summarise(distinct_affiliation = n_distinct(final_affiliation), .groups = 'drop') %>%
  tidyr::complete(substrate, Class, fill = list(distinct_affiliation = 0))

#Range of distinct affiliation within positive samples per primers
edna_ppooled  %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(sample, primer) %>%
  summarise(distinct_affiliation =  n_distinct(final_affiliation),
            .groups = "drop") %>%
  group_by(primer) %>%
  summarise(max_affiliation = max(distinct_affiliation),
            min_affiliation = min (distinct_affiliation))

# Example of OTU affiliation list
edna_ppooled %>%
  filter(sum_positive_replicate > 0) %>%      
  group_by(primer) %>%                        
  distinct(final_affiliation) %>%             
  ungroup() 

#Number of samples without non-human vertebrate detection
edna_ppooled %>%
  filter(final_affiliation != "Homo_sapiens") %>%      
  group_by(sample, primer) %>%
  summarise(total_sample_reads = sum(sum_reads),
            .groups = "drop") %>%
  filter(total_sample_reads == 0) %>%                 
  group_by(primer) %>%
  summarise(nbr_sample = n())                         

edna_gpooled %>%
  filter(final_affiliation != "Homo_sapiens") %>%      
  group_by(sample) %>%
  summarise(total_sample_reads = sum(sum_reads),
            .groups = "drop") %>%
  filter(total_sample_reads == 0) %>% 
  summarise(nbr_sample = n()) 

#Number of sample with domestic vertebrate detected
edna_gpooled %>%
  filter(final_affiliation %in% domestic_affiliation_names) %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(final_affiliation) %>%
  summarise(positive_sample = n())

#Total number of detection events for each Class
edna_gfiltered %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(Class) %>%
  summarise(detection_event = n())

#Total number of detection events for each primer, substrate
edna_pfiltered %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(primer, substrate) %>%
  summarise(detection_event = n())

#Number of sample per primer substrate after filtering 
edna_pfiltered %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(primer, substrate) %>%
  summarise(positive_sample = n_distinct(sample))

# Means of taxa per substrates
edna_gfiltered %>%
  mutate(positive = if_else(sum_reads >0 , 1, sum_reads)) %>%
  group_by(sample, substrate) %>%
  summarise(identified_taxa = sum(positive),
            .groups = "drop") %>%
  ungroup() %>%
  group_by(substrate) %>%
  summarise(mean = sum(identified_taxa) / n() ,
            .groups = "drop")


## Euleur Plots ----

#Choose data used for euleur plots ▲ 
data_euler_sub <- edna_gfiltered

data_euler_pri <- edna_ppooled_taxmixed 

### For substrates ----

#Identify unique affiliation per substrates
affiliations_by_substrate <- data_euler_sub %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(substrate) %>%
  summarise(
    unique_affiliations = list(unique(final_affiliation))
  ) %>%
  ungroup()

#Create a data frame to store results of combinations
intersections <- tibble(
  substrates = character(),
  intersect_affiliations = list(),
  count_intersect_affiliations = integer()
)

#Loop through combinations in decreasing order of complexity
substrate_combinations <- c(
  list(c("leaf", "soil", "spiderweb")),  
  combn(affiliations_by_substrate$substrate, 2, simplify = FALSE),  #Pairs
  combn(affiliations_by_substrate$substrate, 1, simplify = FALSE)  #Singles
)

for (subset in substrate_combinations) {
  #Get intersecting affiliations for each combination
  subset_data <- affiliations_by_substrate %>%
    filter(substrate %in% subset) %>%
    pull(unique_affiliations)
  
  intersect_affiliations <- reduce(subset_data, intersect)
  
  #Exclude affiliations already assigned to higher-order combinations
  already_assigned <- unlist(intersections$intersect_affiliations)
  intersect_affiliations <- setdiff(intersect_affiliations, already_assigned)
  
  #Keep the results
  intersections <- bind_rows(intersections, tibble(
    substrates = paste(subset, collapse = " & "),
    intersect_affiliations = list(intersect_affiliations),
    count_intersect_affiliations = length(unlist(intersect_affiliations))
  ))
}

# Combine results of single substrates, pairs, and higher-order intersections
substrates_area <- intersections
print(substrates_area)

#Area per known substrates :
sp <- substrates_area$count_intersect_affiliations[substrates_area$substrates == "spiderweb"]
lf <- substrates_area$count_intersect_affiliations[substrates_area$substrates == "leaf"]
so <- substrates_area$count_intersect_affiliations[substrates_area$substrates == "soil"]
splf <- substrates_area$count_intersect_affiliations[substrates_area$substrates == "leaf & spiderweb"]
lfso <- substrates_area$count_intersect_affiliations[substrates_area$substrates == "leaf & soil"]
spso <- substrates_area$count_intersect_affiliations[substrates_area$substrates == "soil & spiderweb"]
splfso <- substrates_area$count_intersect_affiliations[substrates_area$substrates == "leaf & soil & spiderweb"]

#Fit euler plot
euler_sub <- eulerr::euler(c(
  spiderweb = sp,
  leaf = lf,
  soil = so,
  "spiderweb&leaf" = splf,
  "leaf&soil" = lfso,
  "spiderweb&soil" = spso,
  "spiderweb&leaf&soil" = splfso
))

rm(sp, lf, so, splf, lfso, spso, splfso)


#Draw euler plot
plot(euler_sub,
     fills = palette_substrate,
     labels = list(
       labels = c("Spiderweb", "Leaf swabs", "Soil"),    
       col = "gray20",                      
       font = 2,                           
       cex = 1.5                           
     ),
     quantities = list(
       col = "gray20",                       
       cex = 1.2                          
     ),
     edges = list(
       lwd = 0,                            
       col = "#66c2a500"                      
     ),
     
)

### For primers ----
affiliations_by_primer <- data_euler_pri %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(primer) %>%
  summarise(
    unique_affiliations = list(unique(final_affiliation))
  ) %>%
  ungroup()

#Obtain affiliation that are only found in each considered substrate
unique_to_each_primer <- map_df(affiliations_by_primer$primer, function(pr) {
  other_affiliations <- affiliations_by_primer %>%
    filter(primer != pr) %>%
    pull(unique_affiliations) %>%
    unlist()
  
  only_in_pr <- setdiff(unlist(affiliations_by_primer %>% filter(primer == pr) %>% pull(unique_affiliations)), other_affiliations)
  
  tibble(
    primer = pr,
    intersect_affiliations = list(only_in_pr),
    count_intersect_affiliations = length(only_in_pr)
  )
})

#List every possible combination of primers
primer_combinations <- lapply(2:nrow(affiliations_by_primer), function(x) {
  combn(affiliations_by_primer$primer, x, simplify = FALSE)
}) %>% unlist(recursive = FALSE)

#Keep only intersect affiliations between combination
intersections <- map_df(primer_combinations, function(subset) {
  subset_data <- affiliations_by_primer %>%
    filter(primer %in% subset) %>%
    pull(unique_affiliations)
  
  intersect_affiliations <- reduce(subset_data, intersect)
  
  tibble(
    primer = paste(subset, collapse = " & "),
    intersect_affiliations = list(intersect_affiliations),
    count_intersect_affiliations = length(unlist(intersect_affiliations))
  )
})

# Combine results of combination with results solely in one substrate
primers_area <- bind_rows(unique_to_each_primer, intersections)
print(primers_area)

# Fit euler plot
euler_prim <- eulerr::euler(c(
  "12S" = primers_area$count_intersect_affiliations[primers_area$primer == "12SV5"],
  "16S" = primers_area$count_intersect_affiliations[primers_area$primer == "16Smam"],
  "12S&16S" = primers_area$count_intersect_affiliations[primers_area$primer == "12SV5 & 16Smam"]
))

#Draw euler plot

plot(euler_prim,
     fills = palette_primers,
     labels = list(
       labels = c("12SV5", "16Smam"),    
       col = "gray20",                      
       font = 2,                           
       cex = 1.5                           
     ),
     quantities = list(
       col = "gray20",                       
       cex = 1.2                          
     ),
     edges = list(
       lwd = 0,                            
       col = "#66c2a500"                      
     )
)
     

#Drop data used for euleur
rm(data_euler_sub)
rm(data_euler_pri)


## Richness distribution ----

# For only positive samples
edna_gpooled  %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(sample, substrate) %>%
  summarise(distinct_affiliation =  n_distinct(final_affiliation),
            .groups = "drop") %>%
  ggplot(aes(x = distinct_affiliation, y = substrate, fill = substrate, color = substrate)) +
  ggridges::geom_density_ridges(alpha = 0.6, scale = 1.5, rel_min_height = 0.01) +
  scale_fill_manual(values = palette_substrate) +
  scale_color_manual(values = palette_substrate) +
  labs(
    title = "Distribution of distinct affiliations across substrates for positive samples",
    x = "Number of distinct affiliations",
    y = "Substrate",
    fill = "Substrate",
    color = "Substrate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )


## Histogram ----

#Choose data to generate hist ▲
data_hist1 <- edna_pfiltered
data_hist2 <- edna_gfiltered

#Generate data used for histogram  
data_hist1 <- data_hist1 %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(primer, substrate, sample, final_affiliation, Class) %>%
  summarise(
    positive_sample = n(),
    .groups = "drop"
  ) %>%
  group_by(primer, substrate, final_affiliation, Class) %>%
  summarise(
    total_count = sum(positive_sample),
    .groups = "drop"
  ) %>%
  arrange(primer, Class, final_affiliation) %>%
  mutate(Class = factor(Class, levels = class_order))

#Plot for primer 12S
hist12s <- data_hist1 %>%
  filter(primer == "12SV5") %>%
  mutate(final_affiliation = forcats::fct_reorder(final_affiliation, total_count, .desc = TRUE)) %>%
  ggplot(aes(x = final_affiliation, y = total_count, fill = substrate)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  labs(
    title = "Primer: 12SV5",
    x = "Final Affiliation ",
    y = "Count of Positive Samples"
  ) +
  facet_grid(cols = vars(Class), scales = "free_x", space = "free_x") +
  scale_fill_manual(values = palette_substrate) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

#Plot for primer 16S
hist16s <- data_hist1 %>%
  filter(primer == "16Smam") %>%
  mutate(final_affiliation = forcats::fct_reorder(final_affiliation, total_count, .desc = TRUE)) %>%
  ggplot(aes(x = final_affiliation, y = total_count, fill = substrate)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  labs(
    title = "Primer: 16Smam",
    x = "Final Affiliation (Grouped by Class)",
    y = "Count of Positive Samples"
  ) +
  facet_grid(cols = vars(Class), scales = "free_x", space = "free_x") +
  scale_fill_manual(values = palette_substrate) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) 

#Combined histogram
patchwork::wrap_plots(hist12s, hist16s, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  plot_annotation(title = "Histogram of positive samples by primer")

rm(data_hist1)


## Stacked Histogram  ----
data_hist2 <- data_hist2 %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(substrate, sample, final_affiliation, Class) %>%
  summarise(
    positive_sample = n(),
    .groups = "drop"
  ) %>%
  group_by(substrate, final_affiliation, Class) %>%
  summarise(
    total_count = sum(positive_sample),
    .groups = "drop"
  ) %>%
  arrange(Class, final_affiliation) %>%
  mutate(Class = factor(Class, levels = class_order)) 

data_hist_combined <- data_hist2 %>%
  mutate(final_affiliation = forcats::fct_reorder(final_affiliation, total_count, .desc = TRUE))

ggplot(data_hist_combined, aes(x = final_affiliation, y = total_count, fill = substrate)) +
  geom_bar(stat = "identity", position = position_stack()) +
  labs(
    title = "Stacked Histogram of Positive Samples by Substrate",
    x = "Final Affiliation",
    y = "Total Count of Positive Samples"
  ) +
  facet_grid(cols = vars(Class), scales = "free_x", space = "free_x") +
  scale_fill_manual(values = palette_substrate) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) 

rm(data_hist2)


## Rarefaction curves ----

#Choose data to use ▲
data_raref_p <- edna_pfiltered
data_raref_g <- edna_gfiltered

### For each primers ----
#Function to generate inext data
generate_inext_data <- function(primer_type, substrates) {
  mapply(function(sub) {
    df <- data_raref_p %>%
      filter(primer == primer_type, substrate == sub) %>%
      mutate(sum_positive_replicate = if_else(sum_positive_replicate > 0, 1, 0)) %>%
      select(sample, final_affiliation, sum_positive_replicate) %>%
      tidyr::pivot_wider(names_from = sample, values_from = sum_positive_replicate) %>%
      as.data.frame()
    rownames(df) <- df$final_affiliation
    df %>% select(-final_affiliation)
  }, sub = substrates, SIMPLIFY = FALSE)
}

#Application of function 
substrates_list <- data_raref_p$substrate %>% unique() 

data_inext_ppooled <- list(
  "12SV5" = generate_inext_data("12SV5", substrates_list),
  "16Smam" = generate_inext_data("16Smam", substrates_list)
)

#Generate rarefaction curves using INEXT
inext_raw_12s <- iNEXT::iNEXT(data_inext_ppooled[["12SV5"]], q = 0, datatype = "incidence_raw")
inext_raw_16s <- iNEXT::iNEXT(data_inext_ppooled[["16Smam"]], q = 0, datatype = "incidence_raw")

#Plot rarefaction curves 
curv_12s <- iNEXT::ggiNEXT(inext_raw_12s, type = 1) +
  scale_color_manual(values = palette_substrate) +
  scale_fill_manual(values = palette_substrate) +
  ggtitle("Rarefaction Curve - 12SV5") +
  theme_minimal() +
  ylim(c(0,60))

curv_16s <- iNEXT::ggiNEXT(inext_raw_16s, type = 1) +
  scale_color_manual(values = palette_substrate) +
  scale_fill_manual(values = palette_substrate) + 
  ggtitle("Rarefaction Curve - 16Smam") +
  theme_minimal() +
  ylim(c(0,60))

patchwork::wrap_plots(curv_12s, curv_16s, ncol = 2) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

rm(curv_12s)
rm(curv_16s)


### For both primers ----

#Generate inext data from filtered pooled primers
data_inext_gpooled <- lapply(substrates_list, function(sub) {
  df <- data_raref_g %>%
    filter(substrate == sub) %>%
    mutate(sum_positive_replicate = if_else(sum_positive_replicate > 0, 1, 0)) %>%
    select(sample, final_affiliation, sum_positive_replicate) %>%
    tidyr::pivot_wider(names_from = sample, values_from = sum_positive_replicate) %>%
    as.data.frame()
  
  rownames(df) <- df$final_affiliation
  df <- df %>% select(-final_affiliation)
  return(df)
})
names(data_inext_gpooled) <- substrates_list

#Generate rarefaction curves using INEXT
inext_raw_gpooled <- iNEXT::iNEXT(data_inext_gpooled, q = 0, datatype = "incidence_raw")

#Extract plot rarefaction curves data 
curv_gpooled <- iNEXT::ggiNEXT(inext_raw_gpooled, type = 1, se = TRUE)
plot_data <- curv_gpooled$data

#Plot rarefaction curves
plot_data$col <- factor(plot_data$col, levels = substrate_order)

curv_gpooled <- ggplot(data = plot_data, aes(x = x, y = y, colour = col)) +
  geom_ribbon(
    aes(ymin = y.lwr, ymax = y.upr, fill = col),
    alpha = 0.2,
    colour = NA 
  ) +
  geom_line(aes(linetype = lty), size = 2.5) +
  geom_point(
    data = plot_data %>%
      filter(Method == "Extrapolation") %>%
      group_by(col) %>%
      slice_min(order_by = x, n = 1),
    aes(shape = col),
    size = 6
  ) +
  scale_colour_manual(
    values = palette_substrate
  ) +
  scale_fill_manual(
    values = palette_substrate
  ) +
  labs(
    title = NULL,
    x = "Number of samples",
    y = "Species richness"
  ) +
  ylim(c(0, 60)) +  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),  
    axis.text = element_text(color = "gray20", size = 11),
    axis.title = element_text(color = "gray20", size = 12),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
  ) +
  guides(
    lty = "none",  
    fill = "none", 
    colour = guide_legend(override.aes = list(
      linetype = "blank",  
      size = 5  
    ))
  )
curv_gpooled

#Save rarefaction curves
ggsave(filename = here::here("figures","rarefaction_curve.pdf"),
       plot = curv_gpooled, 
       width = 2000, height = 2000, units = "px", device = "pdf", dpi = 300, bg = NULL)


#Remove no more required elements
rm(curv_gpooled)

rm(data_raref_p)
rm(data_raref_g)


# #Simple way to do similar rarefaction curve 
# curv_gpooled <- iNEXT::ggiNEXT(inext_raw_gpooled, type = 1) +
#   scale_color_manual(values = palette_substrate) +
#   scale_fill_manual(values = palette_substrate) +
#   ggtitle("Rarefaction Curve - Pooled 12Sv5 - 16Smam") +
#   ylim(c(0, 60)) +
#   theme_minimal(base_size = 14) +
#   theme(
#     axis.text = element_text(color = "gray20", size = 12),          
#     axis.title = element_text(color = "gray20", size = 14), 
#     legend.position = "bottom"
#   ) 