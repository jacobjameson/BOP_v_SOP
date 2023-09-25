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


