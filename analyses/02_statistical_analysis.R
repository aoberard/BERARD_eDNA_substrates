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
