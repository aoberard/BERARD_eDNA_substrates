 # A FAIRE APRES DISCU MAX ET NATH ------
 
# modeles 0/1 - proba de detaction
# modele 1,2, 3 - replicabilité
 # faire sur data pooled, peut-être prendre en compte primer en tant que facteur aleatoire
# dans ce cas, pas forcement besoin de prendre données poolés au niveau de l'echantillon, si ?

 # reflechir a si dans ce cas c'est genant le fait d'avoir pool (notamment car certaines affiliations existe pas pour certains primer, pas meme nombre de ligne au total (pas meme nombre de 0 - 1)), mais peut-être si met en facteur aleatoire c'est pg ?
 # le fait est que pool réduit nombre total de taxa (car pool plus filtr ensuite qui rabote le pool)
 
 
 
 
 
# GLMMs ----
class_to_ignore <- c("Lepidosauria", "Amphibia")
 

## Considering both primers ----

### Model for detection probability ----

#Create data used in the glm
d_glm_detect <- edna_gfiltered %>%
  mutate(detection = if_else(sum_positive_replicate > 0, 1, 0)) %>%
  filter (!Class %in% class_to_ignore)

#Global model specification
m_detect <- lme4::glmer(data = d_glm_detect,
                              formula = detection ~ substrate * Class + (1|sample) + (1|final_affiliation),
                              family = binomial(link = "logit"),
                              na.action = "na.fail",
                              control = lme4::glmerControl( optimizer="bobyqa", optCtrl=list(maxfun=2e5) ) )

#Model selection
model_selection <- MuMIn::dredge(m_detect, rank = "AICc")
model_selection %>% filter(delta <2)

#Best model specification
m_detect <- lme4::glmer(data = d_glm_detect,
                              formula = detection ~ substrate * Class + (1|sample) + (1|final_affiliation),
                              family = binomial(link = "logit"),
                              na.action = "na.fail",
                              control = lme4::glmerControl( optimizer = "bobyqa", optCtrl = list(maxfun=2e5) ) )

#Best model validation
DHARMa::simulateResiduals(m_detect) %>%
  DHARMa::testResiduals()

#Best model look
summary(m_detect)

em <- emmeans::emmeans(m_detect, ~ substrate | Class)
plot(em, comparisons = TRUE)
emmeans::contrast(em, "pairwise", adjust = "Tukey")

ggstats::ggcoef_model(m_detect)


### Model for repeatability ----

#Create data used in the glm
d_glm_replic <- edna_gfiltered %>%
  filter(sum_positive_replicate > 0) %>%
  filter (!Class %in% class_to_ignore)

hist(d_glm_replic$pooled_percent_positive)
d_glm_replic$pooled_percent_positive %>%
  unique()

#Global model specification
m_replic <- lme4::glmer(data = d_glm_replic,
                              formula = (pooled_percent_positive) ~ substrate * Class + (1|sample) ,
                              family = binomial(link = "logit"),
                              weights = pooled_number,
                              na.action = "na.fail",
                              control = lme4::glmerControl( optimizer = "bobyqa", optCtrl = list(maxfun=2e5) ) )

#Model selection
model_selection <- MuMIn::dredge(m_replic, rank = "AICc")
model_selection %>% filter(delta <2)

#Best model specification
m_replic <- lme4::glmer(data = d_glm_replic,
                              formula = (pooled_percent_positive) ~ substrate + Class + (1|sample) ,
                              family = binomial(link = "logit"),
                              weights = pooled_number,
                              na.action = "na.fail",
                              control = lme4::glmerControl( optimizer = "bobyqa", optCtrl = list(maxfun=2e5) ) )

#Best model validation
DHARMa::simulateResiduals(m_replic) %>%
  DHARMa::testResiduals()



########################## pas terriblement ajustement

#beta
m_replic <- glmmTMB::glmmTMB(pooled_percent_positive ~ substrate + Class + (1|sample),
                               family = betabinomial(link = "logit"),
                               weights = pooled_number,
                               data = d_glm_replic)


summary(m_replic)
 
em <- emmeans::emmeans(m_replic, ~ substrate | Class)
plot(em, comparisons = TRUE)
emmeans::contrast(em, "pairwise", adjust = "Tukey")

ggstats::ggcoef_model(m_replic)





m_replic <- lme4::glmer(data = d_glm_replic,
                        formula = (pooled_percent_positive) ~ substrate * Class + (1|sample) ,
                        family = quasibinomial,
                        weights = pooled_number,
                        na.action = "na.fail",
                        control = lme4::glmerControl( optimizer = "bobyqa", optCtrl = list(maxfun=2e5) ) )







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
TukeyStep1<-aov(detections ~ substrate) #(step1)
TukeyHSD(TukeyStep1,'substrate') #(step2)
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
