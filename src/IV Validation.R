#=========================================================================
# Purpose: Validate IV + Main Results
# Author: Jacob Jameson 
#=========================================================================

library(lfe)
library(stargazer)
library(texreg)
library(xtable)
library(tidyverse)

data <- read_csv('outputs/data/final.csv')

##########################################################################
#=========================================================================
# Establishing IV Validity -----------------------------------------------
#=========================================================================
## Relevance
## Exclusion
## Monotonicity
#=========================================================================
##########################################################################

##########################################################################
#=========================================================================
# Relevance --------------------------------------------------------------
## First Stage results
#=========================================================================
##########################################################################

# Shift-level FE
first_stage_final <- felm(
  any.batch ~ batch.tendency + test.inclination | 
              dayofweekt + month_of_year |0| ED_PROVIDER, 
  data = data)

# Shift-level + complaint FE
first_stage_final_d <- felm(
  any.batch ~ batch.tendency + test.inclination | 
              dayofweekt + month_of_year + age_groups + 
              complaint_esi | 0 | ED_PROVIDER,
  data = data)

# Shift-level + complaint + individual FE
first_stage_final_d_i <- felm(
  any.batch ~ batch.tendency + test.inclination | 
              dayofweekt + month_of_year + age_groups + 
              complaint_esi + GENDER + race | 0 | ED_PROVIDER,
  data = data)

# Save the results to a .txt file
sink("outputs/tables/First Stage.txt")

screenreg(list(first_stage_final,
               first_stage_final_d,
               first_stage_final_d_i), 
          include.fstatistic = T)

stargazer(list(first_stage_final,
               first_stage_final_d,
               first_stage_final_d_i),
          type = "text", header = FALSE, 
          title = "First Stage", style = 'QJE')

sink()

##########################################################################
#=========================================================================
# Reduced-Form Results ---------------------------------------------------
#=========================================================================
##########################################################################

# Model 1 does not control testing inclination
model1.LOS <- felm(
  ln_ED_LOS ~ batch.tendency + patients_in_hospital | 
              dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data = data)

model1.ntest <- felm(
  nEDTests ~ batch.tendency + patients_in_hospital | 
             dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data = data)

model1.72 <- felm(
  RTN_72_HR ~ batch.tendency  + patients_in_hospital | 
              dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data = data)

# Model 2 controls for testing inclination
model2.LOS <- felm(
  ln_ED_LOS ~ batch.tendency + test.inclination + patients_in_hospital | 
              dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data = data)

model2.ntest <- felm(
  nEDTests ~ batch.tendency + test.inclination + patients_in_hospital | 
             dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data = data)

model2.72 <- felm(
  RTN_72_HR ~ batch.tendency + test.inclination + patients_in_hospital| 
              dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data = data)

#=========================================================================
# All coefficients are scaled by the difference in tendency/inclination
# between the ninetieth and tenth percentile physicians for 
# interpretability. 
#=========================================================================

percentile_10.b <- quantile(data$batch.tendency, probs = 0.10)
percentile_90.b <- quantile(data$batch.tendency, probs = 0.90)

percentile_10.t <- quantile(data$test.inclination, probs = 0.10)
percentile_90.t <- quantile(data$test.inclination, probs = 0.90)

coeffecient.scale.b <- percentile_90.b - percentile_10.b
coeffecient.scale.t <- percentile_90.t - percentile_10.t

# Save the results to a .txt file
sink("outputs/tables/Reduced Form.txt")

stargazer(list(model1.LOS, 
               model1.ntest,
               model1.72), type = "text", 
          header = FALSE, title = "Reduced Form", style = 'QJE')

stargazer(list(model2.LOS, 
               model2.ntest,
               model2.72), type = "text", 
          header = FALSE, title = "Reduced Form", style = 'QJE')



# Function for scaling the coefficients in regressions for 
# interpretability
scale_coefficients <- function(model, scale_factors) {

  scaled_model <- model
  
  scaled_model$coefficients[1, "Estimate"] <-
    model$coefficients[1, "Estimate"] * scale_factors[1]
  
  scaled_model$coefficients[1, "Cluster s.e."] <- 
    model$coefficients[1, "Cluster s.e."] * scale_factors[1]
  
  if (length(scale_factors) == 2) {
    
    scaled_model$coefficients[2, "Estimate"] <- 
      model$coefficients[2, "Estimate"] * scale_factors[2]
    
    scaled_model$coefficients[2, "Cluster s.e."] <-
      model$coefficients[2, "Cluster s.e."] * scale_factors[2]
  }
  
  return(scaled_model)
}

scale_factors.1 <- c(coeffecient.scale.b, coeffecient.scale.t)
scale_factors.2 <- c(coeffecient.scale.b)

# Scaled regression results
summary(model1.LOS)$call
scale_coefficients(
  summary(model1.LOS), 
  scale_factors.1)$coefficients[,c('Estimate', 'Cluster s.e.')]

summary(model1.ntest)$call
scale_coefficients(
  summary(model1.ntest), 
  scale_factors.1)$coefficients[,c('Estimate', 'Cluster s.e.')]

summary(model1.72)$call
scale_coefficients(
  summary(model1.72), 
  scale_factors.1)$coefficients[,c('Estimate', 'Cluster s.e.')]

summary(model2.LOS)$call
scale_coefficients(
  summary(model2.LOS), 
  scale_factors.2)$coefficients[,c('Estimate', 'Cluster s.e.')]

summary(model2.ntest)$call
scale_coefficients(
  summary(model2.ntest), 
  scale_factors.2)$coefficients[,c('Estimate', 'Cluster s.e.')]

summary(model2.72)$call
scale_coefficients(
  summary(model2.72), 
  scale_factors.2)$coefficients[,c('Estimate', 'Cluster s.e.')]

sink()

##########################################################################
#=========================================================================
# Exclusion --------------------------------------------------------------
## Placebo Check
#=========================================================================
##########################################################################

data <- read_csv('outputs/data/all_clean.csv')

placebo.complaints <- data %>%
  group_by(CHIEF_COMPLAINT) %>%
  summarise(mean.batch = mean(any.batch), n=n()) %>%
  filter(mean.batch <= 0.1)

placebo.complaints <- placebo.complaints$CHIEF_COMPLAINT

placebo <- data %>%
  filter(CHIEF_COMPLAINT %in% placebo.complaints)

# Save the results to a .txt file
sink("outputs/tables/Placebo Check.txt")

#=========================================================================
# Time controls only
#=========================================================================

placebo1.LOS <- felm(
  ln_ED_LOS ~ batch.tendency + test.inclination | 
              dayofweekt + month_of_year |0| ED_PROVIDER, 
  data = placebo)

placebo1.ntest <- felm(
  nEDTests ~ batch.tendency + test.inclination | 
              dayofweekt + month_of_year |0| ED_PROVIDER, 
  data = placebo)

placebo1.72 <- felm(
  RTN_72_HR ~ batch.tendency + test.inclination | 
              dayofweekt + month_of_year |0| ED_PROVIDER, 
  data = placebo)


stargazer(list(placebo1.LOS,
               placebo1.ntest,
               placebo1.72),
          type = "text", header = FALSE, 
          title = "Placebo Check-- no controls", style = 'QJE')

#=========================================================================
# Time controls + complaint severity
#=========================================================================

placebo2.LOS <- felm(
  ln_ED_LOS ~ batch.tendency + test.inclination | 
              dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data = placebo)

placebo2.ntest <- felm(
  nEDTests ~ batch.tendency + test.inclination | 
             dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data = placebo)

placebo2.72 <- felm(
  RTN_72_HR ~ batch.tendency + test.inclination | 
              dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data = placebo)

stargazer(list(placebo2.LOS,
               placebo2.ntest,
               placebo2.72),
          type = "text", header = FALSE, 
          title = "Placebo Check-- Complaint", style = 'QJE')

#=========================================================================
# Time controls + complaint + person level + volums
#=========================================================================

placebo3.LOS <- felm(
  ln_ED_LOS ~ batch.tendency + test.inclination + 
              patients_in_hospital | 
              dayofweekt + month_of_year + complaint_esi + 
              race + GENDER |0| ED_PROVIDER, 
  data = placebo)

placebo3.ntest <- felm(
  nEDTests ~ batch.tendency + test.inclination + 
             patients_in_hospital | 
             dayofweekt + month_of_year + complaint_esi + 
             race + GENDER |0| ED_PROVIDER, 
  data = placebo)

placebo3.72 <- felm(
  RTN_72_HR ~ batch.tendency + test.inclination + 
              patients_in_hospital | 
              dayofweekt + month_of_year + complaint_esi + 
              race + GENDER |0| ED_PROVIDER, 
  data = placebo)

stargazer(list(placebo3.LOS,
               placebo3.ntest,
               placebo3.72),
          type = "text", header = FALSE, 
          title = "Placebo Check-- all controls", style = 'QJE')


sink()

##########################################################################
#=========================================================================
# Monotonicity -----------------------------------------------------------
## first stage should be weakly positive for all subsamples 
## (Dobbie, Goldin, and Yang 2018)

## instrument constructed by leaving out a particular subsample
## has predictive power over that same left-out subsample 
## (Bhuller et al. 2020).
#=========================================================================
##########################################################################

data <- read_csv('outputs/data/final.csv')

complaints <- unique(data$CHIEF_COMPLAINT)

model_results <- list()
# Iterate through each complaint and run the regression

for(complaint in complaints) {
  
  data_subset <- data %>% 
    filter(CHIEF_COMPLAINT == complaint)
  
  model <- felm(
    any.batch ~ batch.tendency | 
                dayofweekt + month_of_year + ESI | 0 | ED_PROVIDER,
    data = data_subset
  )
  model_results[[complaint]] <- model
}

# Save the results to a .txt file
sink("outputs/tables/Monotonicity.txt")

print("first stage should be weakly positive for 
      all subsamples (Dobbie, Goldin, and Yang 2018)")

stargazer(model_results[1:5], type = "text", 
          title = "Regression Results by Complaint",
          column.labels = complaints[1:5],
          style = 'QJE')

stargazer(model_results[5:10], type = "text", 
          title = "Regression Results by Complaint",
          column.labels = complaints[5:10],
          style = 'QJE')

stargazer(model_results[10:15], type = "text", 
          title = "Regression Results by Complaint",
          column.labels = complaints[10:15],
          style = 'QJE')


model_results <- list()

for(complaint in complaints) {
  
  data_subset <- data %>% 
    filter(CHIEF_COMPLAINT == complaint)
  
  data_subset$residual_batch <- resid(
    felm(any.batch ~ 0 | dayofweekt + month_of_year, data=data_subset))
  
  data_subset <- data_subset %>%
    group_by(ED_PROVIDER) %>%
    mutate(Sum_Resid=sum(residual_batch, na.rm=T),
           batch.tendency = (Sum_Resid - residual_batch) / (n() - 1))
  
  model <- felm(
    any.batch ~ batch.tendency | 
      dayofweekt + month_of_year + ESI | 0 | ED_PROVIDER,
    data = data_subset
  )
  model_results[[complaint]] <- model
}

print("instrument constructed by leaving out a particular subsample
      has predictive power over that same left-out subsample 
      (Bhuller et al. 2020).")


stargazer(model_results[1:5], type = "text", 
          title = "Regression Results by Complaint",
          column.labels = complaints[1:5],
          style = 'QJE')

stargazer(model_results[5:10], type = "text", 
          title = "Regression Results by Complaint",
          column.labels = complaints[5:10],
          style = 'QJE')

stargazer(model_results[10:15], type = "text", 
          title = "Regression Results by Complaint",
          column.labels = complaints[10:15],
          style = 'QJE')

sink()

##########################################################################
#=========================================================================
# IV Results ------------------------------------------------------------
#=========================================================================
##########################################################################

# Main Results
model.IV.LOS <- felm(
  ln_ED_LOS ~ test.inclination | dayofweekt + month_of_year + complaint_esi | 
              (any.batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

model.IV.ntest <- felm(
  nEDTests ~ test.inclination | dayofweekt + month_of_year + complaint_esi | 
            (any.batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)
       
model.IV.72 <- felm(
  RTN_72_HR ~ test.inclination | dayofweekt + month_of_year + complaint_esi | 
             (any.batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)


# Save the results to a .txt file
sink("outputs/tables/2SLS Results.txt")

stargazer(list(model.IV.LOS,
               model.IV.ntest,
               model.IV.72), type = "text", 
          covariate.labels = c('Testing Inclination', 'Batch'),
          dep.var.labels = c('Log LOS', 'Number of Tests', '72 Hour Return',
                             'Admission'),
          header = F, 
          title = "2SLS Results: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')

stargazer(list(model.IV.LOS,
               model.IV.ntest,
               model.IV.72), type = "latex", 
          covariate.labels = c('Testing Inclination', 'Batch'),
          dep.var.labels = c('Log LOS', 'Number of Tests', '72 Hour Return'),
          header = F, 
          title = "2SLS Results: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')

sink()


#=========================================================================
# Exploring certain subgroups
#=========================================================================

# Function to map specific complaints to general categories
map_complaints <- function(complaint) {
  if (complaint %in% c('Abdominal Complaints', 'Gastrointestinal Issues')) {
    return('Gastrointestinal/Abdominal Issues')
  } else if (complaint %in% c('Chest Pain', 'Cardiac Arrhythmias')) {
    return('Cardiac/Chest-Related Issues')
  } else if (complaint %in% c('Shortness of Breath', 'Upper Respiratory Symptoms')) {
    return('Respiratory-Related Issues')
  } else if (complaint %in% c('Neurological Issue', 'Dizziness/Lightheadedness/Syncope')) {
    return('Neurological/Syncope Issues')
  } else if (complaint %in% c('Extremity Complaints', 'Back or Flank Pain', 
                              'Falls, Motor Vehicle Crashes, Assaults, and Trauma')) {
    return('Musculoskeletal/Extremity Issues')
  } else {
    return('General/Other Symptoms')
  }
}

table(data$CHIEF_COMPLAINT)
# Apply the function to create a new column for grouped complaints
data$GROUPED_COMPLAINT <- sapply(data$CHIEF_COMPLAINT, map_complaints)

# Check the new grouping
table(data$GROUPED_COMPLAINT)

# Unique chief complaints
unique_complaints <- unique(data$GROUPED_COMPLAINT)
 data %>%
   filter(GROUPED_COMPLAINT == 'Musculoskeletal/Extremity Issues') %>%
   group_by(batch, sequenced, any.batch) %>%
   summarize(n = n()) %>% view()


# Save the results to a .txt file
sink("outputs/tables/2SLS Results Complaints.txt")

# Loop over each unique complaint
for(complaint in unique_complaints) {
  
  subset <- filter(data, GROUPED_COMPLAINT == complaint)
  
  # Run your regression models on this subset
  model.IV.LOS <- felm(ln_ED_LOS ~ test.inclination | dayofweekt + month_of_year + ESI | 
                         (lab_image_batch ~ batch.tendency + test.inclination)| ED_PROVIDER, data = subset)
  
  model.IV.ntest <- felm(nEDTests ~ test.inclination | dayofweekt + month_of_year + ESI | 
                           (lab_image_batch ~ batch.tendency + test.inclination)| ED_PROVIDER, data = subset)
  
  model.IV.72 <- felm(RTN_72_HR ~ test.inclination | dayofweekt + month_of_year + ESI | 
                        (lab_image_batch ~ batch.tendency + test.inclination)| ED_PROVIDER, data = subset)

  # Generate and print the stargazer table for the current complaint
  stargazer(list(model.IV.LOS, model.IV.ntest, model.IV.72), 
            type = "text", 
            covariate.labels = c('Testing Inclination', 'Batch'),
            dep.var.labels = c('Log LOS', 'Number of Tests', '72 Hour Return'),
            header = F, 
            title = paste("2SLS Results for", complaint, ": Length of Stay, Number of Tests, and 72-Hour"), 
            style = 'QJE')
}

sink()

#=========================================================================
# Specific Types of Test Batching: Lab+Image Batch
#=========================================================================

# Save the results to a .txt file
sink("outputs/tables/2SLS Results Specific Batch Types.txt")

model.IV.LOS <- felm(
  ln_ED_LOS ~ test.inclination | dayofweekt + month_of_year + complaint_esi | 
    (lab_image_batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

model.IV.ntest <- felm(
  nEDTests ~ test.inclination | dayofweekt + month_of_year + complaint_esi | 
    (lab_image_batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

model.IV.72 <- felm(
  RTN_72_HR ~ test.inclination | dayofweekt + month_of_year + complaint_esi | 
    (lab_image_batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)


stargazer(list(model.IV.LOS,
               model.IV.ntest,
               model.IV.72), type = "text", 
          covariate.labels = c('Testing Inclination', 'Batch'),
          dep.var.labels = c('Log LOS', 'Number of Tests', '72 Hour Return'),
          header = F, 
          title = "2SLS Results for Lab+Image Batch: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')


stargazer(list(model.IV.LOS,
               model.IV.ntest,
               model.IV.72), type = "latex", 
          covariate.labels = c('Testing Inclination', 'Batch'),
          dep.var.labels = c('Log LOS', 'Number of Tests', '72 Hour Return'),
          header = F, 
          title = "2SLS Results for Lab+Image Batch: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')

#=========================================================================
# Specific Types of Test Batching: Image+Image Batch
#=========================================================================

model.IV.LOS <- felm(
  ln_ED_LOS ~ test.inclination | dayofweekt + month_of_year + complaint_esi | 
    (image_image_batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

model.IV.ntest <- felm(
  nEDTests ~ test.inclination | dayofweekt + month_of_year + complaint_esi | 
    (image_image_batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

model.IV.72 <- felm(
  RTN_72_HR ~ test.inclination | dayofweekt + month_of_year + complaint_esi | 
    (image_image_batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)


stargazer(list(model.IV.LOS,
               model.IV.ntest,
               model.IV.72), type = "text", 
          covariate.labels = c('Testing Inclination', 'Batch'),
          dep.var.labels = c('Log LOS', 'Number of Tests', '72 Hour Return'),
          header = F, 
          title = "2SLS Results for Image+Image Batch: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')



stargazer(list(model.IV.LOS,
               model.IV.ntest,
               model.IV.72), type = "latex", 
          covariate.labels = c('Testing Inclination', 'Batch'),
          dep.var.labels = c('Log LOS', 'Number of Tests', '72 Hour Return'),
          header = F, 
          title = "2SLS Results for Image+Image Batch: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')



model.IV.LOS <- felm(
  ln_ED_LOS ~ test.inclination | dayofweekt + month_of_year + complaint_esi | 
    (lab_image_batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

model.IV.ntest <- felm(
  nEDTests ~ test.inclination | dayofweekt + month_of_year + complaint_esi | 
    (lab_image_batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

model.IV.72 <- felm(
  RTN_72_HR ~ test.inclination | dayofweekt + month_of_year + complaint_esi | 
    (lab_image_batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

stargazer(list(model.IV.LOS,
               model.IV.ntest,
               model.IV.72), type = "text", 
          covariate.labels = c('Testing Inclination', 'Batch'),
          dep.var.labels = c('Log LOS', 'Number of Tests', '72 Hour Return'),
          header = F, 
          title = "2SLS Results for Lab+Image Batch: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')



stargazer(list(model.IV.LOS,
               model.IV.ntest,
               model.IV.72), type = "latex", 
          covariate.labels = c('Testing Inclination', 'Batch'),
          dep.var.labels = c('Log LOS', 'Number of Tests', '72 Hour Return'),
          header = F, 
          title = "2SLS Results for Lab+Image Batch: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')

sink()
