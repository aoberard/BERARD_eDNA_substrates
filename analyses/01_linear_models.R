
wild_ppooled <- edna_ppooled %>%
  filter(!final_affiliation %in% c("Homo_sapiens", "Sardina_pilchardus")) %>%
  filter(in_zoo == FALSE)

retirer humain , animaux conta(sardine), animaux domestique, animaux zoo
#### PEUT ETRE ICI OU AVANT GENERER DES DATAFRAME DIFFERENT SELON LE DEGRE DE FILTRE ??
### Peut etre biend 'avoir enregistré precedent dans derived_data avant





# Data exploration ----

library(ggplot2)
# Summarize the data and reorder substrate
summary_data <- edna_ppooled %>%
  filter(!final_affiliation %in% c("Homo_sapiens", "Sardina_pilchardus")) %>%
  filter(in_zoo == FALSE) %>%
  group_by(affiliation_level, primer, substrate, sum_positive_replicate) %>%
  summarize(count = n(), .groups = 'drop') %>%
  mutate(substrate = factor(substrate, levels = c("soil", "leaf", "spiderweb")))

# Plot
ggplot(summary_data, aes(x = sum_positive_replicate, y = count, fill = substrate)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(affiliation_level ~ primer) +
  labs(
    title = "Distribution of Positive Replicates by Affiliation, Primer, and Substrate",
    x = "Sum of Positive Replicates",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Summarize the data and reorder substrate
summary_data <- edna_ppooled %>%
  filter(!final_affiliation %in% c("Homo_sapiens", "Sardina_pilchardus")) %>%
  filter(in_zoo == FALSE) %>%
  group_by(affiliation_level, primer, substrate, sum_positive_replicate) %>%
  summarize(count = n(), .groups = 'drop') %>%
  filter(sum_positive_replicate > 0) %>%
  mutate(substrate = factor(substrate, levels = c("soil", "leaf", "spiderweb")))

# Plot
ggplot(summary_data, aes(x = sum_positive_replicate, y = count, fill = substrate)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(affiliation_level ~ primer) +
  labs(
    title = "Distribution of Positive Replicates by Affiliation, Primer, and Substrate",
    x = "Sum of Positive Replicates",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Summarize the data and reorder substrate
summary_data <- edna_ppooled %>%
  filter(domestic == FALSE) %>%
  filter(in_zoo == FALSE) %>%
  filter(!final_affiliation %in% c("Homo_sapiens", "Sardina_pilchardus")) %>%
  group_by(affiliation_level, primer, substrate, sum_positive_replicate) %>%
  summarize(count = n(), .groups = 'drop') %>%
  filter(sum_positive_replicate > 0) %>%
  mutate(substrate = factor(substrate, levels = c("soil", "leaf", "spiderweb")))

# Plot
ggplot(summary_data, aes(x = sum_positive_replicate, y = count, fill = substrate)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(affiliation_level ~ primer) +
  labs(
    title = "Distribution of Positive Replicates by Affiliation, Primer, and Substrate",
    x = "Sum of Positive Replicates",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


#stats descriptives ----

# Total final_affiliation of wild vertebrates per primers
edna_ppooled %>%
  filter(domestic == FALSE) %>%
  filter(in_zoo == FALSE) %>%
  filter(!final_affiliation %in% c("Homo_sapiens", "Sardina_pilchardus")) %>%
  filter(primer == "12SV5") %>%
  pull(final_affiliation) %>%
  unique() %>%
  length()

edna_ppooled %>%
  filter(domestic == FALSE) %>%
  filter(in_zoo == FALSE) %>%
  filter(!final_affiliation %in% c("Homo_sapiens", "Sardina_pilchardus")) %>%
  filter(primer == "16Smamm") %>%
  pull(final_affiliation) %>%
  unique() %>%
  length()

#nbr affiliation a espece ou genre par substrat
edna_ppooled %>%
  filter(domestic == FALSE) %>%
  filter(in_zoo == FALSE) %>%
  filter(!final_affiliation %in% c("Homo_sapiens", "Sardina_pilchardus")) %>%
  filter(affiliation_level %in% c("genus","species")) %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(primer, substrate) %>%
  summarise( nb_detection =  (length(unique(final_affiliation))))


# Compare mean number of positive clusters per sample, per substrate (per primers)

#Calculate number of positive clusters per sample - substrate - primer
summary_detections_12 <- edna_ppooled %>%
  filter(domestic == FALSE) %>%
  filter(in_zoo == FALSE) %>%
  filter(!final_affiliation %in% c("Homo_sapiens", "Sardina_pilchardus")) %>%
  filter(affiliation_level %in% c("genus","species")) %>%
  filter(primer == "12SV5") %>%
  group_by(primer, substrate, sample) %>%
  summarise(detections = sum(sum_positive_replicate > 0)) 


#Post-hoc comparisons (comparaisons 2 à 2 pour variables à plus de 2 modalités)
attach(summary_detections_12)
TukeyStep1<-aov(detections ~ substrate) #utiliser test de Tukey pour obtenir différence entre groupes (step1)
TukeyHSD(TukeyStep1,'substrate') #utiliser test de Tukey pour obtenir différence entre groupes (step2)
pairwise.wilcox.test(detections, substrate,p.adjust.method = "holm", paired = FALSE) #for p-values
detach(summary_detections_12)


summary_detections_16<- edna_ppooled %>%
  filter(domestic == FALSE) %>%
  filter(in_zoo == FALSE) %>%
  filter(!final_affiliation %in% c("Homo_sapiens", "Sardina_pilchardus")) %>%
  filter(affiliation_level %in% c("genus","species")) %>%
  filter(primer == "16Smamm") %>%
  group_by(primer, substrate, sample) %>%
  summarise(detections = sum(sum_positive_replicate > 0)) 

#Post-hoc comparisons (comparaisons 2 à 2 pour variables à plus de 2 modalités)
attach(summary_detections_16)
TukeyStep1<-aov(detections ~ substrate) #utiliser test de Tukey pour obtenir différence entre groupes (step1)
TukeyHSD(TukeyStep1,'substrate') #utiliser test de Tukey pour obtenir différence entre groupes (step2)
pairwise.wilcox.test(detections, substrate,p.adjust.method = "holm", paired = FALSE) #for p-values
detach(summary_detections_16)


#### ATTENTION  biaisé en partie par le fait que echantillon manquants ... (verif quel data utilisé de base)
##### + moyen de simplifer le sccript beaucoup là


summary_detections <- edna_ppooled %>%
  filter(domestic == FALSE) %>%
  filter(in_zoo == FALSE) %>%
  filter(!final_affiliation %in% c("Homo_sapiens", "Sardina_pilchardus")) %>%
  filter(affiliation_level %in% c("genus", "species")) %>%
  group_by(primer, substrate, sample) %>%
  summarise(detections = sum(sum_positive_replicate > 0), .groups = 'drop')

ggplot(summary_detections, aes(x = substrate, y = detections, fill = substrate)) +
  geom_boxplot() +
  facet_grid(~ primer) +  
  labs(
    title = "Number of Detections per samples",
    x = "Substrate",
    y = "Number of Detections"
  ) +
  theme_minimal() +
  theme(legend.position = "none")





# Mean number of positive replicates per sample - substrates - primer
# juste penser que sera biaiser tant que pas de correction sample replicas sautés ☻


# Dataframe pour regarder nb moyen de replicas positifs, LORSQUE échantillon positif pour un cluster
summary_detections <- edna_ppooled %>%
  filter(domestic == FALSE) %>%
  filter(in_zoo == FALSE) %>%
  filter(!final_affiliation %in% c("Homo_sapiens", "Sardina_pilchardus")) %>%
  filter(affiliation_level %in% c("genus","species")) %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(primer, substrate, sample) 

ggplot(summary_detections, aes(x = substrate, y = sum_positive_replicate, fill = substrate)) +
  geom_boxplot() +
  facet_grid(~ primer) +  
  labs(
    title = "Positive replicate per cluster (only positive cluster)",
    x = "Substrate",
    y = "Positive replicates"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

attach(summary_detections %>% filter(primer == "12SV5"))
TukeyStep1<-aov(sum_positive_replicate ~ substrate) #utiliser test de Tukey pour obtenir différence entre groupes (step1)
TukeyHSD(TukeyStep1,'substrate') #utiliser test de Tukey pour obtenir différence entre groupes (step2)
pairwise.wilcox.test(sum_positive_replicate, substrate,p.adjust.method = "holm", paired = FALSE) #for p-values
detach(summary_detections_16)

attach(summary_detections %>% filter(primer == "16Smamm"))
TukeyStep1<-aov(sum_positive_replicate ~ substrate) #utiliser test de Tukey pour obtenir différence entre groupes (step1)
TukeyHSD(TukeyStep1,'substrate') #utiliser test de Tukey pour obtenir différence entre groupes (step2)
pairwise.wilcox.test(sum_positive_replicate, substrate,p.adjust.method = "holm", paired = FALSE) #for p-values
detach(summary_detections_16)


puis test wilcox








# Diagramme de Venn

# peut-être extraire liste taxon pour chaque substrat, puis faire des listes intersections
# puis determiner le nombre de chaque liste et ensuite 

# Là avec tous taxons (dont au dessus de genre et espece :)


  


# Euleur Plots ----

library(purrr)
library(combinat)

## For substrates ----

#Step1: Identify unique affiliation per substrates
affiliations_by_substrate <- wild_ppooled %>%
  filter(primer == "12SV5") %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(substrate) %>%
  summarise(
    unique_affiliations = list(unique(final_affiliation))
  ) %>%
  ungroup()

#Step2 : obtain affiliation that are only found in each considered substrate
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

#Step3: List every possible combination of substrates
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

# Fit euler plot
fit <- eulerr::euler(c(
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
plot(fit,
     labels = c("Spiderweb", "Leaf", "Soil"),
     quantities = TRUE)

plot(fit,
     fills = c("#66c2a5AA", "#fc8d62AA", "#FFD700AA"),
     labels = list(
       labels = c("Spiderweb", "Leaf", "Soil"),    
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


     

## For primers ----
affiliations_by_primer <- wild_ppooled %>%
  filter(sum_positive_replicate > 0) %>%
  group_by(primer) %>%
  summarise(
    unique_affiliations = list(unique(final_affiliation))
  ) %>%
  ungroup()

#Step2 : obtain affiliation that are only found in each considered substrate
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

#List every possible combination primers
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
fit <- eulerr::euler(c(
  "12S" = primers_area$count_intersect_affiliations[primers_area$primer == "12SV5"],
  "16S" = primers_area$count_intersect_affiliations[primers_area$primer == "16Smamm"],
  "12S&16S" = primers_area$count_intersect_affiliations[primers_area$primer == "12SV5 & 16Smamm"]
))

#Draw euler plot
plot(fit,
     fills = c("#66c2a5AA", "#fc8d62AA"),
     labels = list(
       labels = c("12SV5", "16Smamm"),    
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





# Rarefaction curves ----

# retrive generate abundance data

wild_ppooled %>% tidyr::pivot_wider(names_from = )

wild_ppooled$
library(tidyr)

# Courbe de rarefaction avec accumcomp
library("BiodiversityR")
# a l'aor de bug sur ordi boulot ça donc peutetre test iNext?

#most analysis pipelines require a community matrix (typically having sites as rows, species as columns and abundance values as cell values) and an environmental data set (typically providing numerical and categorical variables for the different sites) as inputs. 
#donc besoin d'un pivot wider probablement


sinon fonction iNEXT permet extrapolation aussi
https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.pdf

library(iNEXT)
iNEXT()

Hill numbers

# courbe rarefaction faisable egalement avec autre indice diversité que richesse, peut-être 
# on integrera a analyse diversité





# MODELE BOURRIN GLM ---- 

#essayer amorces separement ou poolés, à échelle de replica ou de l'echantillon (echelle replica relou potetre si pas colonnes filtres faites debut)?

# pas oublier de tester facteurs organisme effet combiné substrat






#OBJECTIFS RESTANTS : 

#pool 12s - 16s et faire analyse pooles (voir diffrence pris differement)
# analyse replicabilité (notamment avec pool amorce pour max voulait voir - mais peut-être avec pourcentage plutot?)
#analyse glm proba detection
#stats descriptives, dont nbr mauvaise quali par substrat (mauvaise qualite si moins de 500 ou 1000 reads par replicas) - a faire sur data avant pool, pour ensemble des taxons par replicats
#diagramme venn pour affiliations amorces -avant et après pooling affiliations?)


