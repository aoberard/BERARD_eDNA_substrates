 # A FAIRE APRES DISCU MAX ET NATH ------
 
# modeles 0/1 - proba de detaction
# modele 1,2, 3 - replicabilité
 # faire sur data pooled, peut-être prendre en compte primer en tant que facteur aleatoire
# dans ce cas, pas forcement besoin de prendre données poolés au niveau de l'echantillon, si ?

 # reflechir a si dans ce cas c'est genant le fait d'avoir pool (notamment car certaines affiliations existe pas pour certains primer, pas meme nombre de ligne au total (pas meme nombre de 0 - 1)), mais peut-être si met en facteur aleatoire c'est pg ?
 # le fait est que pool réduit nombre total de taxa (car pool plus filtr ensuite qui rabote le pool)
 
 
# GLM ----
 

 
 
 ## Amorces ensembles ----
 
 
# data glm creation ----
class_to_ignore <- c("Lepidosauria", "Amphibia")
 
d_glm_gpool_detect <- edna_gfiltered %>%
   mutate(detection = if_else(sum_positive_replicate > 0, 1, 0)) %>%
   filter (!Class %in% class_to_ignore)
 
d_glm_gpool_replic <- edna_gfiltered %>%
  filter(sum_positive_replicate > 0) %>%
  filter (!Class %in% class_to_ignore)

 
#modele logit

# detection

m_gpool_detect <- lme4::glmer(data = d_glm_gpool_detect,
                              formula = detection ~ substrate * Class + (1|sample) + (1|final_affiliation),
                              family = binomial(link = "logit"),
                              na.action = "na.fail",
                              control = lme4::glmerControl( optimizer="bobyqa", optCtrl=list(maxfun=2e5) ) )

model_selection <- MuMIn::dredge(m_gpool_detect, rank = "AICc")
model_selection %>% filter(delta <2)

DHARMa::simulateResiduals(m_gpool_detect) %>%
  DHARMa::testResiduals()

summary(m_gpool_detect)














# Repeatability

  hist(d_glm_gpool_replic$sum_positive_replicate)


m_gpool_replic <- lme4::glmer(data = d_glm_gpool_replic,
                              formula = sum_positive_replicate ~ substrate * Class + (1|sample) + (1|final_affiliation) + (1|primer),
                              family = poisson(link="log"),
                              na.action="na.fail",
                              control = lme4::glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5) ) )

model_selection <- MuMIn::dredge(m_gpool_replic, rank = "AICc")
model_selection %>% filter(delta <2)

m_gpool_replic <- lme4::glmer(data = d_glm_gpool_replic,
                              formula = sum_positive_replicate ~ substrate  + (1|sample) + (1|final_affiliation) + (1|primer),
                              family = poisson(link="log"),
                              na.action="na.fail",
                              control = lme4::glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5) ) )

DHARMa::simulateResiduals(m_gpool_replic) %>%
  DHARMa::testResiduals()

summary(m_gpool_replic)
 










m_gpool_replic <-NA
m_gpool_replic <- lme4::glmer(data = d_glm_gpool_replic,
                              formula = (pooled_percent_positive) ~ substrate * Class + (1|sample) ,
                              family = binomial(link = "logit"),
                              weights = pooled_number,
                              na.action = "na.fail",
                              control = lme4::glmerControl( optimizer="bobyqa", optCtrl=list(maxfun=2e5) ) )

model_selection <- MuMIn::dredge(m_gpool_replic, rank = "AICc")
model_selection %>% filter(delta <2)

m_gpool_replic <- lme4::glmer(data = d_glm_gpool_replic,
                              formula = (pooled_percent_positive) ~ substrate + Class + (1|sample) ,
                              family = binomial(link = "logit"),
                              weights = pooled_number,
                              na.action = "na.fail",
                              control = lme4::glmerControl( optimizer="bobyqa", optCtrl=list(maxfun=2e5) ) )

summary(m_gpool_replic)











## Amorces séparement ----
 
 
### 12s d'abord ----
 
# data glm creation ----
d_glm_12s <- edna_pfiltered %>%
   mutate(detection = if_else(sum_positive_replicate > 0, 1, 0)) %>%
   filter(primer == "16Smam") %>%
  filter (!Class %in% c("Lepidosauria", "Amphibia"))

  
 
# MODELE LOGIT
 
m_logit_12s <- lme4::glmer(data = d_glm_12s,
                            formula = detection ~ substrate * Class + (1|sample),
                            family = binomial(link = "logit"),
                            na.action = "na.fail",
                            control = lme4::glmerControl( optimizer="bobyqa", optCtrl=list(maxfun=2e5) ) )
 
model_selection <- MuMIn::dredge(m_logit_12s, rank = "AICc")
model_selection %>% filter(delta <2)

DHARMa::simulateResiduals(m_logit_12s) %>%
  DHARMa::testResiduals()

summary(m_logit_12s)

em <- emmeans::emmeans(m_logit_12s, ~ substrate | Class)
em
plot(em, comparisons = TRUE)
emmeans::contrast(em, "pairwise", adjust = "Tukey")


# MODELE RICHESSE ?

 
#Calculate number of positive clusters per sample - substrate - primer
summary_detections_12 <- edna_pfiltered %>%
  filter(primer == "12SV5") %>%
  group_by(primer, substrate, sample) %>%
  summarise(richness = sum(sum_positive_replicate > 0)) 
 
#

m_rich_12 <- lme4::glmer(
  data = summary_detections_12,
  formula = richness ~ substrate + (1|sample),
  family = poisson(link="log"),
  na.action="na.fail",
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5))
)
 
model_selection <- MuMIn::dredge(m_rich_12, rank = "AICc")
model_selection %>% filter(delta <2)

DHARMa::simulateResiduals(m_rich_12) %>%
  DHARMa::testResiduals()

m_rich_12 |>
  RVAideMemoire::overdisp.glmer()
 
summary(m_rich_12)
 
 
 
 
   
   
   
   
# TEST AU PIF
 
d_glm <- edna_pfiltered %>%
  mutate(detection = if_else(sum_positive_replicate > 0, 1, 0)) %>%
  filter (!Class %in% c("Lepidosauria", "Amphibia"))


m_logit <- lme4::glmer(data = d_glm,
                           formula = detection ~ substrate * Class + primer + (1|sample),
                           family = binomial(link = "logit"),
                           na.action = "na.fail",
                           control = lme4::glmerControl( optimizer="bobyqa", optCtrl=list(maxfun=2e5) ) )

model_selection <- MuMIn::dredge(m_logit, rank = "AICc")
model_selection %>% filter(delta <2)
 

summary(m_logit)

















































# Sens de variation par tests de moyennes  ----

 voir si besoin de changer type de moyenne selon approche glm employée
 
 

# Compare mean number of positive clusters per sample, per substrate (per primers)

#Calculate number of positive clusters per sample - substrate - primer
summary_detections_12 <- edna_pfiltered %>%
  filter(primer == "12SV5") %>%
  group_by(primer, substrate, sample) %>%
  summarise(detections = sum(sum_positive_replicate > 0)) 


#Post-hoc comparisons (comparaisons 2 à 2 pour variables à plus de 2 modalités)
attach(summary_detections_12)
TukeyStep1<-aov(detections ~ substrate) #utiliser test de Tukey pour obtenir différence entre groupes (step1)
TukeyHSD(TukeyStep1,'substrate') #utiliser test de Tukey pour obtenir différence entre groupes (step2)
pairwise.wilcox.test(detections, substrate,p.adjust.method = "holm", paired = FALSE) #for p-values
detach(summary_detections_12)


summary_detections_16<- edna_pfiltered %>%
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


summary_detections <- edna_pfiltered %>%
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
summary_detections <- edna_pfiltered %>%
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
detach(summary_detections %>% filter(primer == "12SV5"))

attach(summary_detections %>% filter(primer == "16Smamm"))
TukeyStep1<-aov(sum_positive_replicate ~ substrate) #utiliser test de Tukey pour obtenir différence entre groupes (step1)
TukeyHSD(TukeyStep1,'substrate') #utiliser test de Tukey pour obtenir différence entre groupes (step2)
pairwise.wilcox.test(sum_positive_replicate, substrate,p.adjust.method = "holm", paired = FALSE) #for p-values
detach(summary_detections %>% filter(primer == "16Smamm"))
