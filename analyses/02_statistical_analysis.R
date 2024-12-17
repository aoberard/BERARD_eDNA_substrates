# Script set-up ----
 
## Library and functions ----
 library(ggplot2)
 library(dplyr)

## Parameters ----
class_to_ignore <- c("Lepidosauria", "Amphibia")

 
# GLMMs ----

## Model for detection probability ----

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
m_detect_best <- lme4::glmer(data = d_glm_detect,
                              formula = detection ~ substrate * Class + (1|sample) + (1|final_affiliation),
                              family = binomial(link = "logit"),
                              na.action = "na.fail",
                              control = lme4::glmerControl( optimizer = "bobyqa", optCtrl = list(maxfun=2e5) ) )

#Best model validation
DHARMa::simulateResiduals(m_detect_best) %>%
  DHARMa::testResiduals()

#Best model look
summary(m_detect_best)
gtsummary::tbl_regression(m_detect_best)
ggstats::ggcoef_model(m_detect_best)

#LRT test
drop1(m_detect_best,.~.,test="Chisq")

#Post-hoc tests
em_detect <- emmeans::emmeans(m_detect_best, ~ substrate | Class, type = "response")
plot(em_detect, comparisons = TRUE)
emmeans::contrast(em_detect, "pairwise", adjust = "Tukey")

#Plot emmeans prob results
em_summary_detect <- as.data.frame(em_detect)

ggplot(em_summary_detect, aes(x = prob, y = reorder(substrate, prob), color = substrate)) +
  geom_point(size = 4) +
  geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, y = substrate, yend = substrate), linewidth = 1.3) +
  scale_color_manual(values = palette_substrate) + 
  facet_wrap(~ Class, nrow = 2) + 
  theme_minimal(base_size = 14) +
  labs(
    x = "Estimated probability",
    y = "Substrate",
    color = "Substrate"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) 


## Model for within sample repeatability ----
#Create data used in the glm (only positive here)
d_glm_repeat <- edna_gfiltered %>%
  filter(sum_reads > 0) 

#Global model specification
m_repeat <- lme4::glmer(data = d_glm_repeat,
                        formula = within_repeated_positive ~ substrate + (1|sample),
                        family = binomial(link = "logit"),
                        na.action = "na.fail",
                        control = lme4::glmerControl( optimizer="bobyqa", optCtrl=list(maxfun=2e5) ) )

#Model selection
model_selection <- MuMIn::dredge(m_repeat, rank = "AICc")
model_selection %>% filter(delta <2)

#Best model specification
m_repeat_best <- lme4::glmer(data = d_glm_repeat,
                        formula = within_repeated_positive ~ substrate + (1|sample),
                        family = binomial(link = "logit"),
                        na.action = "na.fail",
                        control = lme4::glmerControl( optimizer="bobyqa", optCtrl=list(maxfun=2e5) ) )

#Best model validation
DHARMa::simulateResiduals(m_repeat_best) %>%
  DHARMa::testResiduals()

#Best model look
summary(m_repeat_best)
gtsummary::tbl_regression(m_repeat_best)
ggstats::ggcoef_model(m_repeat_best)

#LRT test
drop1(m_repeat_best,.~.,test="Chisq")

#Post-hoc tests
em_repeat <- emmeans::emmeans(m_repeat_best, ~ substrate , type = "response")
plot(em_repeat, comparisons = TRUE)
emmeans::contrast(em_repeat, "pairwise", adjust = "Tukey")

#Plot emmeans prob results
em_summary_repeat <- as.data.frame(em_repeat)

ggplot(em_summary_repeat, aes(x = prob, y = reorder(substrate, prob), color = substrate)) +
  geom_point(size = 4) +
  geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, y = substrate, yend = substrate), linewidth = 1.3) +
  scale_color_manual(values = palette_substrate) + 
  theme_minimal(base_size = 14) +
  labs(
    x = "Estimated probability",
    y = "Substrate",
    color = "Substrate"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) 











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
