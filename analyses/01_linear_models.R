


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


mais pareil biaisé en partie par le fait que echantillon manquants ...
moyen de simplifer le sccript beaucoup là


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


puis test wilcox




# MODELE BOURRIN GLM ---- 

# essayer avec et sans avoir rajouter colonnes sample et replicas qui ont sauté 

