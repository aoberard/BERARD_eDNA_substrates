# Script set-up ----

## Library and functions ----
library(purrr)
library(combinat)
library(ggplot2)
library(dplyr)
library(patchwork)


## Graphical parameters ----

#Define colors for each substrate
palette_substrate <- c("spiderweb" = "#66c2a5AA", "leaf" = "#fc8d62AA", "soil" = "#FFD700AA")
palette_primers <- c("12SV5" = "#FFB3B3AA", "16Smam" = "#A2D1D1AA")


# Data quality ----

## Considering reads number ----
#We are seeking for the number of replicates and samples per primer and substrate under possible quality thresholds

#Choose data on which to check quality (global data with small modifications)
data_quality_check <- all_edna %>%
  filter(!stringr::str_detect(sample, "ZOO")) 
  # %>% filter(final_affiliation == "Homo_sapiens")       # possible filtering for human (ubiquitous contamination)

#Calculate possible reads threshold values
threshold_values <- seq(0,
                        data_quality_check %>%
                          group_by(PCRreplicate, primer) %>%
                          summarise(total_reads = sum(reads),
                                    .groups = "drop" ) %>%
                          pull(total_reads) %>%
                          max(),
                        by = 100)

#For each threshold_values, calculate replicate and sample counts under it
quality_variation <- map_df(threshold_values, function(threshold) {
  data_quality_check %>%
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
    ungroup() %>%
    filter(total_reads <= threshold) %>%
    group_by(primer, substrate) %>%
    summarise(
      replicate_count = n(),             
      sample_count = n_distinct(sample),
      .groups = "drop"
    ) %>%
    ungroup() %>%
    tidyr::complete(primer, substrate, fill = list(replicate_count = 0, sample_count = 0)) %>%
    mutate(reads_quality_threshold = threshold) 
})

#Plot replicate count against reads threshold
ggplot(quality_variation, aes(x = reads_quality_threshold, y = replicate_count, color = substrate)) +
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
ggplot(quality_variation, aes(x = reads_quality_threshold, y = sample_count, color = substrate)) +
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
quality_variation %>%
  filter(reads_quality_threshold == chosen_quality_threshold)




# Final taxonomic affiliations ----

## Descriptive numbers ----

#Total number of distinct affiliation  per primers
edna_ppooled %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(primer) %>%
  summarise(distinct_affiliation =  n_distinct(final_affiliation))

#After chosen filters total number of distinct affiliation  per primers
edna_pfiltered %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(primer) %>%
  summarise(distinct_affiliation =  n_distinct(final_affiliation))

#Total number of distinct affiliation  per substrate - primers
edna_ppooled %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(substrate, primer) %>%
  summarise(distinct_affiliation =  n_distinct(final_affiliation)) 

# After filters total number of distinct affiliation  per substrate - primers
edna_pfiltered %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(substrate, primer) %>%
  summarise(distinct_affiliation =  n_distinct(final_affiliation))

# Number of distinct affiliation per taxa Class and per substrates
edna_ppooled %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(substrate, Class) %>%
  summarise(distinct_affiliation = n_distinct(final_affiliation), .groups = 'drop') %>%
  tidyr::complete(substrate, Class, fill = list(distinct_affiliation = 0))

# After filters number of distinct affiliation per taxa Class and per substrates
edna_pfiltered %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(substrate, Class) %>%
  summarise(distinct_affiliation = n_distinct(final_affiliation), .groups = 'drop') %>%
  tidyr::complete(substrate, Class, fill = list(distinct_affiliation = 0))


## Euleur Plots ----

#Choose data used for euleur plots ▲ 
data_euler_sub <- edna_gpooled

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

#obtain affiliation that are only found in each considered substrate
unique_to_each_substrate <- map_df(affiliations_by_substrate$substrate, function(sub) {
  other_affiliations <- affiliations_by_substrate %>%
    filter(substrate != sub) %>%
    pull(unique_affiliations) %>%
    unlist()
  
  only_in_sub <- setdiff(unlist(affiliations_by_substrate %>% filter(substrate == sub) %>% pull(unique_affiliations)), other_affiliations)
  
  tibble(
    substrates = sub,
    intersect_affiliations = list(only_in_sub),
    count_intersect_affiliations = length(only_in_sub)
  )
})

#List every possible combination of substrates
substrate_combinations <- lapply(2:nrow(affiliations_by_substrate), function(x) {
  combn(affiliations_by_substrate$substrate, x, simplify = FALSE)
}) %>% unlist(recursive = FALSE)

#Keep only intersect affiliations between combination
intersections <- map_df(substrate_combinations, function(subset) {
  subset_data <- affiliations_by_substrate %>%
    filter(substrate %in% subset) %>%
    pull(unique_affiliations)
  
  intersect_affiliations <- reduce(subset_data, intersect)
  
  tibble(
    substrates = paste(subset, collapse = " & "),
    intersect_affiliations = list(intersect_affiliations),
    count_intersect_affiliations = length(unlist(intersect_affiliations))
  )
})

#Combine results of combination with results solely in one substrate
substrates_area <- bind_rows(unique_to_each_substrate, intersections)
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


## Histogram ----

#Generate data used for histogram ▲ 
data_hist <- edna_pfiltered %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(primer, substrate, sample, final_affiliation) %>%
  summarise(positive_sample = n(), .groups = "drop") %>%
  group_by(primer, substrate, final_affiliation) %>%
  summarise(total_count = sum(positive_sample), .groups = "drop") 

#Plot for primer 12S
hist12s <- data_hist %>%
  filter(primer == "12SV5") %>%
  arrange(primer, desc(total_count)) %>%
  mutate(final_affiliation = factor(final_affiliation, levels = unique(final_affiliation))) %>%
  ggplot(aes(x = final_affiliation, y = total_count, fill = substrate)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  labs(x = "Final Affiliation", y = "Count of Positive Samples") +
  facet_grid(rows = vars(primer)) +
  scale_fill_manual(values = palette_substrate) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "bottom"
  ) +
  ylim(0, 10)

#Plot for primer 16S
hist16s <- data_hist %>%
  filter(primer == "16Smam") %>%
  arrange(primer, desc(total_count)) %>%
  mutate(final_affiliation = factor(final_affiliation, levels = unique(final_affiliation))) %>%
  ggplot(aes(x = final_affiliation, y = total_count, fill = substrate)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  labs(x = "Final Affiliation", y = "Count of Positive Samples") +
  facet_grid(rows = vars(primer)) +
  scale_fill_manual(values = palette_substrate) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "bottom",
  ) +
  ylim(0, 10)

#Combined histogram
patchwork::wrap_plots(hist12s, hist16s, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

rm(hist12s)
rm(hist16s)


## Rarefaction curves ----

#Choose data to use
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

#Plot rarefaction curves 
curv_gpooled <- iNEXT::ggiNEXT(inext_raw_gpooled, type = 1) +
  scale_color_manual(values = palette_substrate) +
  scale_fill_manual(values = palette_substrate) +
  ggtitle("Rarefaction Curve - Pooled 12Sv5 - 16Smam") +
  theme_minimal() +
  ylim(c(0,60))

curv_gpooled
rm(curv_gpooled)

rm(data_raref_p)
rm(data_raref_g)


# Repeatability ----

# reste a clean cette partie !!!!!!!!!!!!!!!!!!!!!!!!!!!----

#Plot number of positive replicates for each substrates and affiliation level 
edna_pfiltered %>%
  group_by(affiliation_level, primer, substrate, sum_positive_replicate) %>%
  summarize(count = n(), .groups = 'drop') %>%
  filter(sum_positive_replicate > 0) %>%
  mutate(substrate = factor(substrate, levels = c("soil", "leaf", "spiderweb"))) %>%
  tidyr::complete(affiliation_level, primer, substrate, sum_positive_replicate, fill = list(count = 0)) %>%
  ggplot(aes(x = sum_positive_replicate, y = count, fill = substrate)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(affiliation_level ~ primer) +
  scale_fill_manual(values = palette_substrate) +
  labs(
    title = "Distribution of Positive Replicates by Affiliation, Primer, and Substrate",
    x = "Number of positive replicates",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")



# pour edn gpooled - il faut utiliser pourcentage de replica positif plutot que nbr de positif car parfois
#3 replicas pools parfois 6 etc vu que parfois taoxn partage ou pas
#########################################

# reflechir, il faut regarder difference de distribution pour chacun, en faisant abstraction de l'effet du nombre ot
# car on sait deja que araigne plus de rpelicat +, mais on veut savoir si difference de tronche (plus souvent de 3 par rapport au total pour chaque substrat p ex , il faut s'affranchir du nombre total)




edna_gfiltered %>%
  mutate(percent_positive_replicate = round(100 * sum_positive_replicate / pooled_number, 1)) %>%
  group_by(affiliation_level, substrate, percent_positive_replicate) %>%
  summarize(count = n(),
            .groups = 'drop') %>%
  filter(percent_positive_replicate > 0) %>%
  mutate(substrate = factor(substrate, levels = c("soil", "leaf", "spiderweb"))) %>%
  tidyr::complete(affiliation_level, substrate, percent_positive_replicate, fill = list(count = 0)) %>%
  ggplot(aes(x = percent_positive_replicate, y = count, fill = substrate)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(rows = vars(affiliation_level)) +
  scale_fill_manual(values = palette_substrate) +
  labs(
    title = "Distribution of Positive Replicates by Affiliation and Substrate",
    x = "Percentage of Positive Replicates",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")







# ESSAYER DE NORMALISER AXE Y, biais lié à nombre 


# Total count per substrate
edna_gfiltered_normalized <- edna_gfiltered %>%
  mutate(percent_positive_replicate = round(100 * sum_positive_replicate / pooled_number, 1)) %>%
  group_by(affiliation_level, substrate, percent_positive_replicate) %>%
  summarize(count = n(), 
            .groups = 'drop') %>%
  filter(percent_positive_replicate > 0) %>%
  mutate(substrate = factor(substrate, levels = c("soil", "leaf", "spiderweb"))) %>%
  tidyr::complete(affiliation_level, substrate, percent_positive_replicate, fill = list(count = 0))

# Normalize the count by the total per substrate
edna_gfiltered_normalized <- edna_gfiltered_normalized %>%
  group_by(substrate) %>%
  mutate(
    total_count_substrate = sum(count),  
    normalized_count = count / total_count_substrate  
  ) %>%
  ungroup()

# Plot the normalized count
edna_gfiltered_normalized %>%
  ggplot(aes(x = percent_positive_replicate, y = normalized_count, fill = substrate)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(rows = vars(affiliation_level)) +
  scale_fill_manual(values = palette_substrate) +
  labs(
    title = "Normalized Distribution of Positive Replicates by Affiliation and Substrate",
    x = "Percentage of Positive Replicates",
    y = "Normalized Count"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  ylim(0,1)






























