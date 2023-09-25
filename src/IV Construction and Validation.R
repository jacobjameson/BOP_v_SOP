#=========================================================================
# Purpose: Construct main IV and Validate
# Author: Jacob Jameson 
#=========================================================================

source('src/clean.R')
library(lfe)
library(stargazer)
library(texreg)

#=========================================================================
# IV Construction --------------------------------------------------------
#=========================================================================

# Step 1: leave-out residualize at the ED encounter level
## conditional on shift-level variation, random assignment
## residual from regression represents physician tendency to batch
final$residual_batch <- resid(
  felm(any.batch ~ 0 | dayofweekt + month, data=final))

# Step 2: get batch tendency for each provider
final <- final %>%
  group_by(ED_PROVIDER) %>%
  mutate(Sum_Resid=sum(residual_batch, na.rm=T),
         batch.tendency = (Sum_Resid - residual_batch) / (n() - 1)) %>% 
  ungroup()

#=========================================================================
# Establishing IV Validity -----------------------------------------------
#=========================================================================
## Relevance
## Exclusion
## Monotonicity
#=========================================================================

#=========================================================================
# Relevance --------------------------------------------------------------
## First Stage results
#=========================================================================

# Shift-level FE
first_stage_final <- felm(any.batch ~ batch.tendency | 
                          dayofweekt + month |
                          0| ED_PROVIDER, 
                          data = final)

# Shift-level + complaint FE
first_stage_final_d <- felm(any.batch ~ batch.tendency | 
                            dayofweekt + month + age_groups + 
                            complaint_esi |
                            0|ED_PROVIDER,
                            data = final)

# Shift-level + complaint + individual FE
first_stage_final_d_i <- felm(any.batch ~ batch.tendency | 
                              dayofweekt + month + age_groups + 
                              complaint_esi + GENDER + race |
                              0| ED_PROVIDER , 
                              data = final)
# Latex First-Stage
stargazer(first_stage_final, 
          first_stage_final_d,
          first_stage_final_d_i, title='')

# Print First-Stage
screenreg(list(first_stage_final,
               first_stage_final_d,
               first_stage_final_d_i), 
          include.fstatistic = T)

#=========================================================================
# Reduced-Form Results ---------------------------------------------------
#=========================================================================


#=========================================================================
# Exclusion --------------------------------------------------------------
## Placebo Check
#=========================================================================

final %>%
  group_by(CHIEF_COMPLAINT) %>%
  summarise(mean.batch = mean(any.batch), n =n()) %>%
  arrange(desc(mean.batch))

# The idea here is that we want to show that being a high-batching
# physician is not associated with other traits of a high-quality
# physician. We are going to show that for Urinary Complaints,
# a complaint area where batching is rare, that there are not preferable
# outcomes associated with being a batcher

placebo <- final %>%
  filter(CHIEF_COMPLAINT == 'Urinary Complaints')

#=========================================================================
# LOS
#=========================================================================
placebo1.1 <- felm(ln_ED_LOS ~ batch.tendency + avg_nEDTests | 
                  dayofweekt + month | 0 | ED_PROVIDER, 
                  data = placebo)

placebo1.2 <- felm(ln_ED_LOS ~ batch.tendency + avg_nEDTests | 
                  dayofweekt + month + ESI | 0 | ED_PROVIDER, 
                  data = placebo)

placebo1.3 <- felm(ln_ED_LOS ~ batch.tendency + avg_nEDTests | 
                  dayofweekt + month + ESI + race + GENDER |0| ED_PROVIDER, 
                  data = placebo)

screenreg(list(placebo1.1, placebo1.2, placebo1.3), 
          include.fstatistic = T)

#=========================================================================
# Number of Tests
#=========================================================================

placebo2.1 <- felm(nEDTests ~ batch.tendency + avg_nEDTests | 
                   dayofweekt + month | 0 | ED_PROVIDER, 
                   data = placebo)

placebo2.2 <- felm(nEDTests ~ batch.tendency + avg_nEDTests | 
                   dayofweekt + month + ESI | 0 | ED_PROVIDER, 
                   data = placebo)

placebo2.3 <- felm(nEDTests ~ batch.tendency + avg_nEDTests | 
                   dayofweekt + month + ESI + race + GENDER |0| ED_PROVIDER, 
                   data = placebo)

screenreg(list(placebo2.1, placebo2.2, placebo2.3), 
          include.fstatistic = T)

#=========================================================================
# 72 HR Return
#=========================================================================

placebo3.1 <- felm(RTN_72_HR ~ batch.tendency + avg_nEDTests | 
                   dayofweekt + month | 0 | ED_PROVIDER, 
                   data = placebo)

placebo3.2 <- felm(RTN_72_HR ~ batch.tendency + avg_nEDTests | 
                   dayofweekt + month + ESI | 0 | ED_PROVIDER, 
                   data = placebo)

placebo3.3 <- felm(RTN_72_HR ~ batch.tendency + avg_nEDTests | 
                   dayofweekt + month + ESI + race + GENDER |0| ED_PROVIDER, 
                   data = placebo)

screenreg(list(placebo3.1, placebo3.2, placebo3.3), 
          include.fstatistic = T)






