#=========================================================================
# Purpose: Validate IV + Main Results
# Author: Jacob Jameson 
#=========================================================================

library(lfe)
library(stargazer)
library(texreg)
library(xtable)
library(tidyverse)
library(sandwich)
library(lmtest)

data <- read_csv('outputs/data/final.csv') %>%
  filter(ARRIVAL_AGE_DI > 18) %>%
  group_by(complaint_esi) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  group_by(CHIEF_COMPLAINT) %>%
  mutate(batchmean = mean(batched)) %>%
  filter(batchmean > 0.05, imaging == 1)

################################################################################
#-------------------------------------------------------------------------------
# Physician factors to explain batching
#-------------------------------------------------------------------------------
################################################################################

# Model 1 ---- Baseline model --------------------------------------------------
model.1 <- felm(
  batched ~ EXPERIENCE + PROVIDER_SEX + patients_tbs + patients_in_hospital | 
    dayofweekt + month_of_year |0| ED_PROVIDER, 
  data = data)

# Model 2 ---- Add patient characteristics -------------------------------------
model.2 <- felm(
  batched ~ EXPERIENCE + PROVIDER_SEX + patients_tbs + patients_in_hospital |
    dayofweekt + month_of_year + 
    complaint_esi + LAB_PERF  | 0 | ED_PROVIDER, 
  data = data)

stargazer(model.1, model.2, type = "text", single.row = TRUE, digits = 3,
          style = "qje", header = T,
          column.labels = c("Model 1", "Model 2", "Model 3"),
          dep.var.caption = "Dependent Variable: Batched", 
          dep.var.labels.include = FALSE)
      
################################################################################

################################################################################
#-------------------------------------------------------------------------------
# First Stage 
#-------------------------------------------------------------------------------
################################################################################

# Model 1 ---- Baseline model --------------------------------------------------
fs.1 <- felm(
  batched ~ batch.tendency | dayofweekt + month_of_year |0| ED_PROVIDER, 
  data = data
  )

# Model 2 ---- Add patient characteristics -------------------------------------
fs.2 <- felm(
  batched ~ batch.tendency | dayofweekt + month_of_year + 
                             complaint_esi + LAB_PERF | 0 | ED_PROVIDER, 
  data = data
  )

screenreg(list(fs.1, fs.2), include.fstatistic = T)
 
stargazer(fs.1, fs.2, type = "text", single.row = TRUE, digits = 3,
          style = "qje", header = T,
          column.labels = c("Model 1", "Model 2"),
          dep.var.caption = "Dependent Variable: Batched")

################################################################################

################################################################################
#-------------------------------------------------------------------------------
# Reduced Form
#-------------------------------------------------------------------------------
################################################################################

# Model 1a ---- ED LOS ---------------------------------------------------------
edlos.1 <- felm(
  ln_ED_LOS ~ batch.tendency | dayofweekt + month_of_year + dispo + 
                               complaint_esi + LAB_PERF |0| ED_PROVIDER, 
  data = data
  )

# Model 2a ---- TIME TO RESULT -------------------------------------------------
ttr.1 <- felm(
  ln_ttr ~ batch.tendency | dayofweekt + month_of_year + 
                            complaint_esi + LAB_PERF |0| ED_PROVIDER, 
  data = data
  )

# Model 3a ---- NUM IMAGES -----------------------------------------------------
n_images.1 <- felm(
  imgTests ~ batch.tendency | dayofweekt + month_of_year + 
                              complaint_esi + LAB_PERF |0| ED_PROVIDER, 
  data = data
  )

# Model 4a ---- 72HRA ----------------------------------------------------------
hr.1 <- felm(
  RTN_72_HR_ADMIT ~ batch.tendency | dayofweekt + month_of_year + 
                                     complaint_esi + LAB_PERF |0| ED_PROVIDER, 
  data = data
  )

stargazer(edlos.1a, ttr.1a, n_images.1a, hr.1a, type = "text", 
          single.row = TRUE, digits = 3,
          style = "qje", header = T,
          column.labels = c("Model 1", "Model 2", "Model 3", "Model 4"),
          dep.var.caption = "Dependent Variable: Batched")


#-------------------------------------------------------------------------------

# Model 1b ---- ED LOS ----------------------------------------------------------
edlos.1b <- felm(
  ln_ED_LOS ~ batch.tendency + test.inclination | dayofweekt + month_of_year  + 
    complaint_esi + LAB_PERF |0| ED_PROVIDER, 
  data = data
)

# Model 2b ---- TIME TO RESULT --------------------------------------------------
ttr.1b <- felm(
  ln_ttr ~ batch.tendency + test.inclination | dayofweekt + month_of_year + 
    complaint_esi + LAB_PERF |0| ED_PROVIDER, 
  data = data
)

# Model 3b ---- NUM IMAGES ------------------------------------------------------
n_images.1b <- felm(
  imgTests ~ batch.tendency + test.inclination | dayofweekt + month_of_year + 
    complaint_esi + LAB_PERF |0| ED_PROVIDER, 
  data = data
)

# Model 4b ---- 72HRA -----------------------------------------------------------
hr.1b <- felm(
  RTN_72_HR_ADMIT ~ batch.tendency + test.inclination | dayofweekt + month_of_year + 
    complaint_esi + LAB_PERF |0| ED_PROVIDER, 
  data = data
)

stargazer(edlos.1b, ttr.1b, n_images.1b, hr.1b, type = "text", 
          single.row = TRUE, digits = 3,
          style = "qje", header = T,
          column.labels = c("Model 1", "Model 2", "Model 3", "Model 4"),
          dep.var.caption = "Dependent Variable: Batched")

################################################################################

################################################################################
#-------------------------------------------------------------------------------
# IV Results
#-------------------------------------------------------------------------------
################################################################################

# Model 1 ---- ED LOS ----------------------------------------------------------

IV_EDLOS <- felm(
  ln_ED_LOS ~ test.inclination | 
              LAB_PERF + dayofweekt + month_of_year + complaint_esi | 
              (batched ~ batch.tendency + test.inclination) |
              ED_PROVIDER, 
  data = data)

summary(IV_EDLOS)

# do 2SLS by hand
fs <- glm(batched ~ batch.tendency + test.inclination + LAB_PERF + dayofweekt + month_of_year + complaint_esi, data = data)
data$batched_hat <- predict(fs)

ss <- glm(ln_ED_LOS ~ batched_hat + test.inclination + LAB_PERF + dayofweekt + month_of_year + complaint_esi, data = data)

# cluster standard errors
coeftest(ss, vcov = vcovHC(ss, type = "HC1"), cluster = ED_PROVIDER)

batched_hat -0.35100793  0.17888814 -1.9622 0.0497434 *  
  
