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

data <- data %>%
  mutate(
    Time_arrival_to_triage = TRIAGE_COMPLETED_REL - ARRIVAL_DTTM_REL,
    Time_triage_to_firstcontact = FIRST_CONTACT_DTTM_REL - TRIAGE_COMPLETED_REL,
    Time_to_Result_Lab = LAB_RESULT_DTTM_REL - LAB_ORDER_REL,
    Time_to_Result_XRay = PLAIN_XRAY_RESULT_DTTM_REL - XR_ORDER_REL,
    Time_to_Result_Ultrasound = US_RESULT_DTTM_REL - US_ORDER_DTTM_REL,
    Time_to_Result_CTcon = CT_WITH_CONTR_RESULT_DTTM_REL - CON_CT_ORDER_REL,
    Time_to_Result_CTnon = CT_WITHOUT_CONTR_RESULT_DTTM_REL - NON_CON_CT_ORDER_REL
  )

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
# Save the results to a .txt file
sink("outputs/tables/First Stage.txt")

# Shift-level FE
first_stage_final <- felm(
  batched ~ batch.tendency | 
              dayofweekt + month_of_year |0| ED_PROVIDER, 
  data = data)

# Shift-level + complaint FE
first_stage_final_d <- felm(
  batched ~ batch.tendency  | 
              dayofweekt + month_of_year + 
              complaint_esi | 0 | ED_PROVIDER, 
  data = data)

# Shift-level + complaint + LAB
first_stage_final_d_i <- felm(
  batched ~ batch.tendency  | 
              dayofweekt + month_of_year + 
              complaint_esi + LAB_PERF | 0 |  ED_PROVIDER, 
  data = data)


screenreg(list(first_stage_final,
               first_stage_final_d,
               first_stage_final_d_i), 
          include.fstatistic = T)

stargazer(list(first_stage_final,
               first_stage_final_d,
               first_stage_final_d_i),
          type = "text", header = FALSE, 
          title = "First Stage", style = 'QJE')

# Shift-level FE
first_stage_final <- glm(
  batched ~ batch.tendency + dayofweekt + month_of_year, 
  family = binomial(link = "probit"),
  data = data)

# Shift-level + complaint FE
first_stage_final_d <- glm(
  batched ~ batch.tendency + dayofweekt + month_of_year + complaint_esi,
  family = binomial(link = "probit"),
  data = data)

# Shift-level + complaint + LAB
first_stage_final_d_i <- glm(
  batched ~ batch.tendency + dayofweekt + month_of_year + 
            complaint_esi + LAB_PERF,
  family = binomial(link = "probit"),
  data = data)

vcov1 <- vcovCL(first_stage_final, cluster = data$ED_PROVIDER)
vcov2 <- vcovCL(first_stage_final_d, cluster = data$ED_PROVIDER)
vcov3 <- vcovCL(first_stage_final_d_i, cluster = data$ED_PROVIDER)

robust_summary1 <- coeftest(first_stage_final, vcov = vcov1)
robust_summary2 <- coeftest(first_stage_final_d, vcov = vcov2)
robust_summary3 <- coeftest(first_stage_final_d_i, vcov = vcov3)


stargazer(first_stage_final, first_stage_final_d, first_stage_final_d_i,
          type = "text", header = FALSE, 
          omit = c('LAB_PERF', "dayofweekt", "month_of_year",
                   "complaint_esi", "GENDER", "race", "Constant"),
          se = list(sqrt(diag(vcov1)), sqrt(diag(vcov2)), sqrt(diag(vcov3))),
          title = "First Stage", style = 'QJE')


sink()

##########################################################################
#=========================================================================
# Reduced-Form Results ---------------------------------------------------
#=========================================================================
##########################################################################

# Save the results to a .txt file
sink("outputs/tables/Reduced Form.txt")

# Model 1 does not control testing inclination --------------------------------
model1.LOS <- felm(
  ln_ED_LOS ~ batch.tendency  | LAB_PERF +imaging +
              dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data)

model1.ntest <- felm(
  imgTests ~ batch.tendency  | LAB_PERF + imaging +
             dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data)

model1.72 <- felm(
  RTN_72_HR_ADMIT ~ batch.tendency | LAB_PERF + imaging +
              dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data)


stargazer(list(model1.72, 
               model1.LOS,
               model1.ntest), type = "text", 
          header = FALSE, title = "Reduced Form", style = 'QJE')

# Model 2 controls for testing inclination ------------------------------------
model2.LOS <- felm(
  ln_ED_LOS ~ batch.tendency + test.inclination  | LAB_PERF  + 
    dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data)

model2.ntest <-  felm(
  imgTests ~ batch.tendency  + test.inclination | LAB_PERF  + 
              dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER , 
  data)


model2.72 <- felm(
  RTN_72_HR_ADMIT ~ batch.tendency + test.inclination | LAB_PERF  + 
              dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data)


stargazer(list(model2.72,
               model2.LOS, 
               model2.ntest), type = "text", 
          header = FALSE, title = "Reduced Form", style = 'QJE')


# Mediation Analysis of Reduced Form -------------------------------------------

mediator_model <- lm(
  imgTests ~ batch.tendency +
              dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
  data = data)

outcome_model.1 <- lm(
  ln_ED_LOS ~ batch.tendency + imgTests +
              dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
  data = data)

med.out <- mediation::mediate(
  mediator_model, outcome_model.1, treat = "batch.tendency", 
  mediator = "imgTests", cluster = data$ED_PROVIDER)

print('Mediation Analysis of Reduced Form')
print('Outcome: ln_ED_LOS')
outcome_model.1$call
print('Mediator: imgTests')
mediator_model$call

summary(med.out)

outcome_model.2 <- lm(
  RTN_72_HR_ADMIT ~ batch.tendency + imgTests +
    dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
  data = data)

print('Mediation Analysis of Reduced Form')
print('Outcome: RTN_72_HR_ADMIT')
outcome_model.1$call
print('Mediator: imgTests')
mediator_model$call

med.out <- mediation::mediate(
  mediator_model, outcome_model.2, treat = "batch.tendency", 
  mediator = "imgTests", cluster = data$ED_PROVIDER)

summary(med.out)

# Mediation Analysis of Reduced Form with testing inc --------------------------

mediator_model <- lm(
  imgTests ~ batch.tendency + test.inclination +
    dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
  data = data)

outcome_model.1 <- lm(
  ln_ED_LOS ~ batch.tendency + imgTests + test.inclination +
    dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
  data = data)

med.out <- mediation::mediate(
  mediator_model, outcome_model.1, treat = "batch.tendency", 
  mediator = "imgTests", cluster = data$ED_PROVIDER)

print('Mediation Analysis of Reduced Form')
print('Outcome: ln_ED_LOS')
outcome_model.1$call
print('Mediator: imgTests')
mediator_model$call

summary(med.out)

outcome_model.2 <- lm(
  RTN_72_HR_ADMIT ~ batch.tendency + imgTests + test.inclination +
    dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
  data = data)

print('Mediation Analysis of Reduced Form')
print('Outcome: RTN_72_HR_ADMIT')
outcome_model.1$call
print('Mediator: imgTests')
mediator_model$call

med.out <- mediation::mediate(
  mediator_model, outcome_model.2, treat = "batch.tendency", 
  mediator = "imgTests", cluster = data$ED_PROVIDER)

summary(med.out)

sink()

##########################################################################
#=========================================================================
# IV Results ------------------------------------------------------------
#=========================================================================
##########################################################################

# Save the results to a .txt file
sink("outputs/tables/IV Results.txt")

# Main Results -------------------------------------------------------------

# first stage
fs <- lm(batched ~ batch.tendency + test.inclination + 
         dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
         data = data)

data$p_batched <- predict(fs)

# second stage

ss <- lm(imgTests ~ p_batched + test.inclination + imaging + 
         dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
         data = data)

# poisson regression
ss <- glm(imgTests ~ p_batched +  imaging + test.inclination +
    dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
    data = data, family = quasipoisson)

coeftest(ss, vcov = vcovCL(ss, cluster = data$ED_PROVIDER))

# negative binomial regression
library(MASS)
ss <- glm.nb(imgTests ~ p_batched + imaging + 
    dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
    data = data)

# cluster SE on ED_PROVIDER
library(lmtest)
library(sandwich)

coeftest(ss, vcov = vcovCL(ss, cluster = data$ED_PROVIDER))



#data <- data %>% filter(ED_LOS > total_testing_time)

time <- felm(
  log(total_testing_time) ~ test.inclination | LAB_PERF + dayofweekt + month_of_year + complaint_esi | 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER , 
  data = data)

summary(time)

LOS <- felm(
  ln_ED_LOS ~ test.inclination | imaging + LAB_PERF + dayofweekt + month_of_year + complaint_esi | 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER , 
  data = data)

ntest <- felm(
  imgTests ~ test.inclination |  imaging + LAB_PERF + dayofweekt + complaint_esi + month_of_year| 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

rtn <- felm(
  RTN_72_HR_ADMIT~ test.inclination | imaging + LAB_PERF + dayofweekt + complaint_esi + month_of_year| 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

stargazer(LOS, ntest, rtn, time,type = "text", 
          header = FALSE, title = "IV Results", style = 'QJE')


ctcon <- felm(
  CON_CT_PERF ~ test.inclination | XR_PERF + NON_CON_CT_PERF + US_PERF + LAB_PERF + dayofweekt + complaint_esi + month_of_year| 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

ctnon <- felm(
  NON_CON_CT_PERF ~ test.inclination | XR_PERF + CON_CT_PERF + US_PERF + LAB_PERF + dayofweekt + complaint_esi + month_of_year| 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

xray <- felm(
  XR_PERF ~ test.inclination |  CON_CT_PERF + NON_CON_CT_PERF + US_PERF + LAB_PERF + dayofweekt + complaint_esi + month_of_year| 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

us <- felm(
  US_PERF ~ test.inclination | CON_CT_PERF + NON_CON_CT_PERF + XR_PERF + LAB_PERF + dayofweekt + complaint_esi + month_of_year| 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

stargazer(ctcon, ctnon, xray, us, type = "text", 
          header = FALSE, title = "IV Results", style = 'QJE')


ctcon <- felm(
  Time_to_Result_CTcon ~ test.inclination | LAB_PERF + dayofweekt + complaint_esi + month_of_year| 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = subset(data, CON_CT_PERF == 1))

ctnon <- felm(
  Time_to_Result_CTnon ~ test.inclination |  LAB_PERF + dayofweekt + complaint_esi + month_of_year| 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = subset(data, NON_CON_CT_PERF == 1))

xray <- felm(
  Time_to_Result_XRay ~ test.inclination | LAB_PERF + dayofweekt + complaint_esi + month_of_year| 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = subset(data, XR_PERF == 1))

us <- felm(
  Time_to_Result_Ultrasound ~ test.inclination | LAB_PERF + dayofweekt + complaint_esi + month_of_year| 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = subset(data, US_PERF == 1))


stargazer(ctcon, ctnon, xray, us, type = "text", 
          header = FALSE, title = "IV Results", style = 'QJE')




first_stage <- lm(batched ~ batch.tendency + test.inclination  + dayofweekt + 
                    month_of_year + LAB_PERF + complaint_esi,  data)

data$p_batched <- predict(first_stage, newdata = data)

# Interaction model to explore the role of ultrasound in batching on ED_LOS
los_model <- lm(ln_ED_LOS ~ p_batched  + US_PERF + CON_CT_PERF + NON_CON_CT_PERF +
                  XR_PERF + LAB_PERF + test.inclination + complaint_esi +
                  dayofweekt + month_of_year, data = data)

summary(los_model)




LOS <- felm(
  ln_ED_LOS ~ test.inclination | LAB_PERF + dayofweekt + month_of_year + complaint_esi | 
    (batched ~ batch.tendency)| ED_PROVIDER, 
  data = data)

ntest <- felm(
  imgTests ~ test.inclination | imaging + LAB_PERF +  dayofweekt + complaint_esi + month_of_year | 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER + complaint_esi, 
  data = data)

summary(ntest)

rtn <- felm(
  RTN_72_HR_ADMIT ~ test.inclination | LAB_PERF + dayofweekt + complaint_esi + month_of_year | 
    (batched ~ batch.tendency) | ED_PROVIDER, 
  data = data)

stargazer(LOS, ntest, rtn, time, type = "text", 
          header = FALSE, title = "IV Results", style = 'QJE')

# manual calculation of the IV results
IV <- function(data){
  
  # first stage
  first_stage_lm <- lm(batched ~ batch.tendency + dayofweekt + test.inclination +
                              month_of_year + complaint_esi + LAB_PERF,  data)
  
  first_stage_probit <- glm(batched ~ batch.tendency + dayofweekt + test.inclination +
                         month_of_year + complaint_esi + LAB_PERF,  data, family = binomial(link = "probit"))
  
  data$p_batched_lm <- predict(first_stage_lm, newdata = data, type = "response")
  data$p_batched_probit <- predict(first_stage_probit, newdata = data, type = "response")
  
  # second stage
  second_stage.1 <- lm(ln_ED_LOS ~ p_batched_lm + LAB_PERF + dayofweekt + 
                         test.inclination +  month_of_year + complaint_esi, data) 
  
  second_stage.2 <- lm(imgTests ~ p_batched_lm + LAB_PERF + dayofweekt + 
                         test.inclination + month_of_year + complaint_esi, data) 
  
  second_stage.3 <- lm(RTN_72_HR_ADMIT ~ p_batched_lm + LAB_PERF + dayofweekt + 
                         test.inclination + month_of_year + complaint_esi, data)
  
  n_clusters <- length(unique(data$ED_PROVIDER))

  # cluster robust standard errors of second stage
  vcov1 <- vcovCL(second_stage.1, cluster = data$ED_PROVIDER)
  robust_summary1 <- coeftest(second_stage.1, vcov = vcov1, type = "HC1", 
                              df = n_clusters - 1)
  
  vcov2 <- vcovCL(second_stage.2, cluster = data$ED_PROVIDER)
  robust_summary3 <- coeftest(second_stage.2, vcov = vcov2, type = "HC1",
                              df = n_clusters - 1)
  
  vcov3 <- vcovCL(second_stage.3, cluster = data$ED_PROVIDER)
  robust_summary3 <- coeftest(second_stage.3, vcov = vcov3, type = "HC1", 
                              df = n_clusters - 1)
  
  stargazer(second_stage.1, second_stage.2, second_stage.3,
            type = "text", header = FALSE, 
            omit = c('LAB_PERF', "dayofweekt", "month_of_year",
                     "complaint_esi", "GENDER", "race", "Constant"),
            se = list(sqrt(diag(vcov1)), sqrt(diag(vcov2)), sqrt(diag(vcov3))),
            title = "Manual 2SLS", style = 'QJE')
  
  
  # second stage
  second_stage.1 <- lm(ln_ED_LOS ~ p_batched_probit + LAB_PERF + dayofweekt + 
                         test.inclination +  month_of_year + complaint_esi, data) 
  second_stage.2 <- lm(imgTests ~ p_batched_probit + LAB_PERF + dayofweekt + 
                         test.inclination + month_of_year + complaint_esi, data) 
  second_stage.3 <- lm(RTN_72_HR_ADMIT ~ p_batched_probit + LAB_PERF + dayofweekt + 
                         test.inclination + month_of_year + complaint_esi, data) 
  
  n_clusters <- length(unique(data$ED_PROVIDER))
  
  # cluster robust standard errors of second stage
  vcov1 <- vcovCL(second_stage.1, cluster = data$ED_PROVIDER)
  robust_summary1 <- coeftest(second_stage.1, vcov = vcov1, type = "HC1", 
                              df = n_clusters - 1)
  
  vcov2 <- vcovCL(second_stage.2, cluster = data$ED_PROVIDER)
  robust_summary3 <- coeftest(second_stage.2, vcov = vcov2, type = "HC1",
                              df = n_clusters - 1)
  
  vcov3 <- vcovCL(second_stage.3, cluster = data$ED_PROVIDER)
  robust_summary3 <- coeftest(second_stage.3, vcov = vcov3, type = "HC1", 
                              df = n_clusters - 1)
  
  stargazer(second_stage.1, second_stage.2, second_stage.3,
            type = "text", header = FALSE, 
            omit = c('LAB_PERF', "dayofweekt", "month_of_year",
                     "complaint_esi", "GENDER", "race", "Constant"),
            se = list(sqrt(diag(vcov1)), sqrt(diag(vcov2)), sqrt(diag(vcov3))),
            title = "Manual 2SLS", style = 'QJE')
}

IV(data = data)

complaints <- data %>%
  group_by(CHIEF_COMPLAINT) %>%
  summarise(n = n()) %>%
  filter(n > 300) %>%
  pull(CHIEF_COMPLAINT) 
complaints

summary(felm(
  imgTests ~ test.inclination | imaging + ESI + LAB_PERF + dayofweekt + month_of_year| 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = subset(data, CHIEF_COMPLAINT == "Dizziness/Lightheadedness/Syncope")))

summary(felm(
  ln_ED_LOS ~ test.inclination | imaging + ESI + LAB_PERF + dayofweekt + month_of_year| 
    (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = subset(data, CHIEF_COMPLAINT == "Dizziness/Lightheadedness/Syncope")))


for (complaint in complaints){
  
  summary(felm(
    imgTests ~ test.inclination | imaging + LAB_PERF + dayofweekt + month_of_year| 
      (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
    data = filter(data, complaint_esi == complaint)))
  
}


data$ESI <- as.factor(data$ESI)


# mediation analysis function 

  
  dat <- data
  dat$total_testing_time_log <- log(dat$total_testing_time)
  # first stage
  first_stage <- lm(batched ~ batch.tendency + test.inclination  + dayofweekt + 
                    month_of_year + LAB_PERF + complaint_esi,  dat)
  
  dat$p_batched <- predict(first_stage, newdata = dat)

  mediator_model <- lm(
    imgTests ~ p_batched + test.inclination  + dayofweekt + month_of_year +
               LAB_PERF + complaint_esi,
    data = dat)
  
  
  outcome_model.1 <- lm(
    total_testing_time_log ~ p_batched + imgTests + test.inclination  +
                dayofweekt + month_of_year  + complaint_esi + 
                LAB_PERF , 
    data = dat)
  
  outcome_model.2 <- lm(
    ln_ED_LOS ~ p_batched + imgTests + test.inclination  +
      dayofweekt + month_of_year  + complaint_esi +  LAB_PERF, 
    data = dat)
  
  med.out1 <- mediation::mediate(
    mediator_model, outcome_model.1, treat = "p_batched", 
    mediator = "imgTests",  cluster = dat$ED_PROVIDER)
  

  print(summary(felm(
    imgTests ~ test.inclination | complaint_esi  + imaging + LAB_PERF + dayofweekt + month_of_year| 
      (batched ~ batch.tendency + test.inclination)| ED_PROVIDER, 
    data = dat)))
  
  print('Mediation Analysis of Reduced Form')
  print('Outcome: ln_ED_LOS')
  outcome_model.1$call
  print('Mediator: imgTests')
  mediator_model$call
  
  print(summary(med.out1))
  
  print('Mediation Analysis of Reduced Form')
  print('Outcome: RTN_72_HR_ADMIT')
  outcome_model.1$call
  print('Mediator: imgTests')
  mediator_model$call
  
  med.out2 <- mediation::mediate(
    mediator_model, outcome_model.2, treat = "p_batched", 
    mediator = "imgTests", cluster = dat$ED_PROVIDER)
  
print(summary(med.out2))



  
sink()

library(survival)

data

##########################################################################
#=========================================================================
# Complaint Analysis  ----------------------------------------------------
#=========================================================================
##########################################################################

# for each unique complaint, estimate the regression model
complaints <- unique(data$CHIEF_COMPLAINT)

for (complaint in complaints){

  data$complaint <- ifelse(data$CHIEF_COMPLAINT == complaint, 1, 0)
  
  # first stage
  first_stage <- lm(batched ~ batch.tendency + test.inclination + dayofweekt + 
                    month_of_year + as.factor(ESI) + LAB_PERF,  
                    subset(data, complaint == 1))
  
  data$p_batched <- predict(first_stage, newdata = data, type = "response")
  
  # second stage
  second_stage.1 <- lm(ln_ED_LOS ~ p_batched + LAB_PERF + dayofweekt + 
                         test.inclination + month_of_year + as.factor(ESI) + LAB_PERF,  
                       subset(data, complaint == 1))
  
  second_stage.2 <- lm(imgTests ~ p_batched + LAB_PERF + dayofweekt + 
                         test.inclination +month_of_year + as.factor(ESI) + LAB_PERF,  
                       subset(data, complaint == 1))
  
  second_stage.3 <- lm(RTN_72_HR_ADMIT ~ p_batched + LAB_PERF + dayofweekt + 
                       test.inclination + month_of_year + as.factor(ESI) + LAB_PERF,  
                       subset(data, complaint == 1))
  
  n_clusters <- length(unique(subset(data, complaint == 1)$ED_PROVIDER))

  # cluster robust standard errors of second stage
  vcov1 <- vcovCL(second_stage.1, cluster = data$ED_PROVIDER)
  
  robust_summary1 <- coeftest(second_stage.1, cluster = data$ED_PROVIDER, type = "HC1", 
                              df = n_clusters - 1)
  
  vcov2 <- vcovCL(second_stage.2, cluster = data$ED_PROVIDER)
  
  robust_summary3 <- coeftest(second_stage.2, vcov = vcov2, type = "HC1",
                              df = n_clusters - 1)
  
  vcov3 <- vcovCL(second_stage.3, cluster = data$ED_PROVIDER)
  
  robust_summary3 <- coeftest(second_stage.3, vcov = vcov3, type = "HC1", 
                              df = n_clusters - 1)
  
  stargazer(second_stage.1, second_stage.2, second_stage.3,
            type = complaint, header = FALSE, 
            omit = c('LAB_PERF', "dayofweekt", "month_of_year", "ESI"),
            title = complaint, style = 'QJE')
  
}



comps <- c("Shortness of Breath", "Fatigue and Weakness", 'Fevers, Sweats or Chills',
           "Falls, Motor Vehicle Crashes, Assaults, and Trauma",
           'Chest Pain', 'Neurological Issue', 'Abdominal Complaints')


comps <- c("Shortness of Breath")

for (complaint in comps){
  
  data$complaint <- ifelse(data$CHIEF_COMPLAINT == complaint, 1, 0)
  
  subd <- subset(data, complaint == 1)
  
  subd$residual_batch <- resid(
    felm(batched ~ 0 | dayofweekt + month_of_year + ESI + LAB_PERF |0|ED_PROVIDER, data=subd)
  )
  subd$residual_ntests <- resid(
    felm(imgTests ~ 0 | dayofweekt + month_of_year + ESI + LAB_PERF |0|ED_PROVIDER, data=subd))
  
  # Step 2: get batch tendency for each provider
  subd <- subd %>%
    group_by(ED_PROVIDER) %>%
    mutate(Sum_Resid=sum(residual_batch, na.rm=T),
           batch.tendency = (Sum_Resid - residual_batch) / (n() - 1),
           
           Sum_Resid_ntests=sum(residual_ntests, na.rm=T),
           test.inclination = (Sum_Resid_ntests - residual_ntests) / (n() - 1)) %>%
    ungroup()
  
  LOS <- felm(
    ln_ED_LOS ~ test.inclination | LAB_PERF + dayofweekt + ESI + complaint_esi | 
      (batched ~ batch.tendency)| ED_PROVIDER, subd)
  
  ntest <- felm(
    imgTests ~ test.inclination | LAB_PERF +  dayofweekt + ESI + month_of_year | 
      (batched ~ batch.tendency)| ED_PROVIDER, subd)
  
  rtn <- felm(
    RTN_72_HR_ADMIT ~ test.inclination | LAB_PERF + dayofweekt + ESI + month_of_year | 
      (batched ~ batch.tendency) | ED_PROVIDER, subd)
  
  stargazer(LOS, ntest, rtn, type = "text", 
            header = FALSE, title = complaint, style = 'QJE')
}



significant <- c('Abdominal Complaints', 'Neurological Issue', 'Falls, Motor Vehicle Crashes, Assaults, and Trauma',
                 'Fevers, Sweats or Chills', 'Neurological Issue', 'Shortness of Breath')





for (complaint in comps){
  
  data$complaint <- ifelse(data$CHIEF_COMPLAINT == complaint, 1, 0)
  
  subd <- subset(data, complaint == 1)
  
  subd$residual_batch <- resid(
    felm(batched ~ 0 | dayofweekt + month_of_year + ESI + LAB_PERF |0|ED_PROVIDER, data=subd)
  )
  subd$residual_ntests <- resid(
    felm(imgTests ~ 0 | dayofweekt + month_of_year + ESI + LAB_PERF |0|ED_PROVIDER, data=subd))
  
  # Step 2: get batch tendency for each provider
  subd <- subd %>%
    group_by(ED_PROVIDER) %>%
    mutate(Sum_Resid=sum(residual_batch, na.rm=T),
           batch.tendency = (Sum_Resid - residual_batch) / (n() - 1),
           
           Sum_Resid_ntests=sum(residual_ntests, na.rm=T),
           test.inclination = (Sum_Resid_ntests - residual_ntests) / (n() - 1)) %>%
    ungroup()
  
  
  fs <- lm(batched ~ batch.tendency + test.inclination + dayofweekt + 
             month_of_year +  as.factor(ESI) + LAB_PERF, data=subd)
  
  subd$pbatch <- predict(fs, newdata = subd, type = "response")
  
  # second stage
  second_stage.1 <- lm(ln_ED_LOS ~ pbatch + LAB_PERF + dayofweekt + 
                         test.inclination + month_of_year + as.factor(ESI) + LAB_PERF,  
                       subd)
  
  second_stage.2 <- lm(imgTests ~ pbatch + LAB_PERF + dayofweekt + 
                         test.inclination +month_of_year + as.factor(ESI) + LAB_PERF,  
                       subd)
  
  second_stage.3 <- lm(RTN_72_HR_ADMIT ~ pbatch + LAB_PERF + dayofweekt + 
                         test.inclination + month_of_year + as.factor(ESI) + LAB_PERF,  
                       subd)
  
  
  n_clusters <- length(subd$ED_PROVIDER)
  
  # cluster robust standard errors of second stage
  vcov1 <- vcovCL(second_stage.1, cluster = subd$ED_PROVIDER)
  
  robust_summary1 <- coeftest(second_stage.1, cluster = subd$ED_PROVIDER, type = "HC1", 
                              df = n_clusters - 1)
  
  vcov2 <- vcovCL(second_stage.2, cluster = subd$ED_PROVIDER)
  
  robust_summary3 <- coeftest(second_stage.2, vcov = vcov2, type = "HC1",
                              df = n_clusters - 1)
  
  vcov3 <- vcovCL(second_stage.3, cluster = subd$ED_PROVIDER)
  
  robust_summary3 <- coeftest(second_stage.3, vcov = vcov3, type = "HC1", 
                              df = n_clusters - 1)
  
  stargazer(second_stage.1, second_stage.2, second_stage.3,
            type = 'text', header = FALSE, title = complaint,
            
            omit = c('LAB_PERF', "dayofweekt", "month_of_year", "ESI"),
            style = 'QJE')

}






################################################################################

# not include imaging --- 

data$ESI <- as.factor(data$ESI)

first_stage <- lm(batched ~ batch.tendency + test.inclination + dayofweekt + 
                    month_of_year + complaint_esi + LAB_PERF, data=data)

data$pbatch <- predict(first_stage, newdata = data, type = "response")

second_stage.1 <- lm(ln_ED_LOS ~ pbatch + LAB_PERF + dayofweekt + test.inclination +
                      month_of_year + complaint_esi, data=data)

second_stage.2 <- lm(imgTests ~ pbatch + LAB_PERF + dayofweekt + test.inclination +
                      month_of_year + complaint_esi , data=data)

second_stage.3 <- lm(RTN_72_HR_ADMIT ~ pbatch + LAB_PERF + dayofweekt + test.inclination +
                      month_of_year + complaint_esi, data=data)

vcov1 <- vcovCL(second_stage.1, cluster = data$ED_PROVIDER)
vcov2 <- vcovCL(second_stage.2, cluster = data$ED_PROVIDER)
vcov3 <- vcovCL(second_stage.3, cluster = data$ED_PROVIDER)

stargazer(second_stage.1, second_stage.2, second_stage.3,
          type = 'latex', header = FALSE, title = 'All Complaints',
          se = list(sqrt(diag(vcov1)), sqrt(diag(vcov2)), sqrt(diag(vcov3))),
          omit = c('LAB_PERF', "dayofweekt", "month_of_year", "complaint_esi",
                   'test.inclination'),
          style = 'QJE')

# ---- include imaging 

ndata <- subset(data, imaging == 1)

first_stage <- lm(batched ~ batch.tendency + test.inclination + dayofweekt + 
                    month_of_year + complaint_esi + LAB_PERF, data=ndata)

ndata$pbatch <- predict(first_stage, newdata = ndata, type = "response")

second_stage.1 <- lm(ln_ED_LOS ~ pbatch + LAB_PERF + dayofweekt + test.inclination +
                       month_of_year + complaint_esi + imaging, data=ndata)

second_stage.2 <- lm(imgTests ~ pbatch + LAB_PERF + dayofweekt + test.inclination +
                       month_of_year + complaint_esi + imaging , data=ndata)

second_stage.3 <- lm(RTN_72_HR_ADMIT ~ pbatch + LAB_PERF + dayofweekt + test.inclination +
                       month_of_year + complaint_esi + imaging, data=ndata)

vcov1 <- vcovCL(second_stage.1, cluster = ndata$ED_PROVIDER)
vcov2 <- vcovCL(second_stage.2, cluster = ndata$ED_PROVIDER)
vcov3 <- vcovCL(second_stage.3, cluster = ndata$ED_PROVIDER)

stargazer(second_stage.1, second_stage.2, second_stage.3,
          type = 'latex', header = FALSE, title = 'All Complaints',
          se = list(sqrt(diag(vcov1)), sqrt(diag(vcov2)), sqrt(diag(vcov3))),
          omit = c('LAB_PERF', "dayofweekt", "month_of_year", "complaint_esi",
                   'test.inclination', 'imaging'),
          style = 'QJE')

mediator_model <- lm(
  imgTests ~ pbatch + test.inclination +
  dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
  data = ndata)

outcome_model.1 <- lm(
  ln_ED_LOS ~ pbatch + imgTests + test.inclination +
    dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
  data = ndata)

med.out <- mediation::mediate(
  mediator_model, outcome_model.1, treat = "pbatch", 
  mediator = "imgTests", cluster = ndata$ED_PROVIDER)

summary(med.out)


#-----

mediator_model <- lm(
  imgTests ~ pbatch + test.inclination +
    dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
  data = data)

outcome_model.1 <- lm(
  ln_ED_LOS ~ pbatch + imgTests + test.inclination +
    dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
  data = data)

med.out <- mediation::mediate(
  mediator_model, outcome_model.1, treat = "pbatch", 
  mediator = "imgTests", cluster = data$ED_PROVIDER)

summary(med.out)

outcome_model.2 <- lm(
  RTN_72_HR_ADMIT ~ pbatch + imgTests + test.inclination +
  dayofweekt + month_of_year + complaint_esi + LAB_PERF, 
  data = data)

med.out2 <- mediation::mediate(
  mediator_model, outcome_model.2, treat = "pbatch", 
  mediator = "imgTests", cluster = data$ED_PROVIDER)

summary(med.out2)


################################################################################










library(sandwich)
library(lmtest)
library(stargazer)

run_2sls_clustered <- function(data, 
                               iv, 
                               treatment, 
                               outcomes, 
                               controls, 
                               cluster_var,
                               relevant_vars) {
  
  # First stage
  fs_formula <- as.formula(paste(treatment, "~", iv, "+", paste(controls, collapse = " + ")))
  fs <- lm(fs_formula, data = data)
  
  # Predicted values
  data$predicted_treatment <- predict(fs, newdata = data, type = "response")
  
  # Second stage for each outcome
  second_stages <- list()
  for (outcome in outcomes) {
    ss_formula <- as.formula(paste(outcome, "~ predicted_treatment +", paste(controls, collapse = " + ")))
    ss <- lm(ss_formula, data = data)
    
    # Calculate clustered standard errors
    clustered_se <- vcovCL(ss, cluster = data[[cluster_var]])
    
    # Get the coefficients with clustered standard errors
    clustered_coef <- coeftest(ss, vcov = clustered_se)
    
    second_stages[[outcome]] <- clustered_coef
  }
  
  # Prepare models for stargazer
  stargazer_models <- lapply(second_stages, function(x) {
    model <- x[relevant_vars, ]
    class(model) <- c("coeftest", "matrix")
    return(model)
  })
  
  # Generate stargazer output
  stargazer_output <- capture.output(
    stargazer(stargazer_models, 
              type = "text", 
              column.labels = outcomes,
              model.names = FALSE,
              covariate.labels = relevant_vars)
  )
  
  # Return results
  return(list(first_stage = fs, 
              second_stages = second_stages, 
              stargazer_output = stargazer_output))
}


stargazer(second_stages, 
          type = "text", 
          column.labels = outcomes,
          model.names = FALSE,
          covariate.labels = relevant_vars)

# Define your variables
data <- subd
iv <- "batch.tendency"
treatment <- "batched"
outcomes <- c("ln_ED_LOS", "imgTests", "RTN_72_HR_ADMIT")
controls <- c("test.inclination", "dayofweekt", "month_of_year", "ESI", "LAB_PERF")
cluster_var <- "ED_PROVIDER"
relevant_vars <- "predicted_treatment"

# Run the 2SLS with clustered SE
results <- run_2sls_clustered(data, iv, treatment, outcomes, controls, cluster_var, relevant_vars)

# Print the stargazer output
cat(results$stargazer_output, sep = "\n")

# If you want to access individual models:
# First stage model
print(summary(results$first_stage))

# Second stage models
print(results$second_stages$ln_ED_LOS)
print(results$second_stages$imgTests)
print(results$second_stages$RTN_72_HR_ADMIT)






























final$residual_batch <- resid(
  felm(batched ~ 0 | dayofweekt + month_of_year + complaint_esi + LAB_PERF |0|ED_PROVIDER, data=final)
)

final$residual_ntests <- resid(
  felm(imaging ~ 0 | dayofweekt + month_of_year + complaint_esi + LAB_PERF |0|ED_PROVIDER, data=final))

# Step 2: get batch tendency for each provider
final <- final %>%
  group_by(ED_PROVIDER) %>%
  mutate(Sum_Resid=sum(residual_batch, na.rm=T),
         batch.tendency = (Sum_Resid - residual_batch) / (n() - 1),
         
         Sum_Resid_ntests=sum(residual_ntests, na.rm=T),
         test.inclination = (Sum_Resid_ntests - residual_ntests) / (n() - 1)) %>%
  ungroup()









##########################################################################
#=========================================================================
# Exclusion --------------------------------------------------------------
## Placebo Check
#=========================================================================
##########################################################################

data <- read_csv('outputs/data/all_clean.csv') 

complaints <- data %>%
  group_by(CHIEF_COMPLAINT) %>%
  summarise(mbatch = mean(batched),
            mimg = mean(imaging)) %>%
  filter(mimg < 0.25 & mbatch < 0.05) %>%
  pull(CHIEF_COMPLAINT) 

datap <- data %>%
  filter(CHIEF_COMPLAINT %in% complaints)
  

# Placebo Check 
model1.LOS <- felm(
  ln_ED_LOS ~ batch.tendency + test.inclination   | LAB_PERF + imaging +
    dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data)

model1.ntest <- felm(
  imgTests ~ batch.tendency + test.inclination  | LAB_PERF + imaging +
    dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data)

model1.72 <- felm(
  RTN_72_HR_ADMIT ~ batch.tendency + test.inclination | LAB_PERF + imaging +
    dayofweekt + month_of_year + complaint_esi |0| ED_PROVIDER, 
  data)


stargazer(list(model1.72, 
               model1.LOS,
               model1.ntest), type = "text", 
          header = FALSE, title = "Reduced Form", style = 'QJE')







placebo.complaints <- data %>%
  filter(ED_LOS > 0) %>%
  mutate(img = ifelse(imgTests > 0, 1, 0)) %>%
  group_by(CHIEF_COMPLAINT) %>%
  summarise(mean.img = mean(img), n=n()) %>%
  filter(mean.img < 0.25) 

placebo.complaints <- placebo.complaints$CHIEF_COMPLAINT
placebo.complaints
placebo <- data %>%
  filter(CHIEF_COMPLAINT %in% placebo.complaints)

# Save the results to a .txt file
sink("outputs/tables/Placebo Check.txt")

#=========================================================================
# Time controls only
#=========================================================================

placebo1.LOS <- felm(
  ln_ED_LOS ~ batch.tendency + test.inclination |  LAB_PERF + 
              dayofweekt + month_of_year |0| ED_PROVIDER, 
  data = placebo)

placebo1.ntest <- felm(
  imgTests ~ batch.tendency + test.inclination | LAB_PERF + 
              dayofweekt + month_of_year |0| ED_PROVIDER, 
  data = placebo)

placebo1.72 <- felm(
  RTN_72_HR ~ batch.tendency + test.inclination |  LAB_PERF + 
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

summary(felm(
  ln_ED_LOS ~ 0 | LAB_PERF + dayofweekt + month_of_year + complaint_esi | 
    (batched ~ batch.tendency)| ED_PROVIDER, 
  data = data))














library(AER)
library(lmtest)
library(sandwich)
library(MASS)  # For probit regression


iv_probit_linear_regression <- function(data, outcome_var, instrument_var, control_vars, cluster_var) {
  
  # First stage: Probit regression
  first_stage_formula <- as.formula(paste("batched ~", instrument_var, "+", paste(control_vars, collapse = " + ")))
  first_stage <- glm(first_stage_formula, family = binomial(link = "probit"), data = data)
  
  # Generate predicted probabilities
  data$predicted_batched_prob <- predict(first_stage, newdata = data, type = "response")
  
  # Second stage: Linear regression
  second_stage_formula <- as.formula(paste(outcome_var, "~ predicted_batched_prob +", paste(control_vars, collapse = " + ")))
  second_stage <- lm(second_stage_formula, data = data)
  
  # Calculate clustered standard errors
  clustered_se <- vcovCL(second_stage, cluster = data[[cluster_var]], type = "HC1")
  
  # Get coefficient, SE, and p-value for predicted_batched_prob
  coef_batched <- coef(second_stage)["predicted_batched_prob"]
  se_batched <- sqrt(diag(clustered_se))["predicted_batched_prob"]
  t_value <- coef_batched / se_batched
  p_value <- 2 * pt(abs(t_value), df = nrow(data) - length(coef(second_stage)), lower.tail = FALSE)
  
  # Get R-squared
  r_squared <- summary(second_stage)$r.squared
  
  # Return results
  return(list(
    coefficient = coef_batched,
    std_error = se_batched,
    p_value = p_value,
    r_squared = r_squared
  ))
}

ln_ED_LOS_result <- iv_probit_linear_regression(
  data = data,
  outcome_var = "ln_ED_LOS",
  instrument_var = "batch.tendency",
  control_vars = c("LAB_PERF", "dayofweekt", "complaint_esi", "month_of_year", 'test.inclination'),
  cluster_var = "ED_PROVIDER"
)
ln_ED_LOS_result

nEDTests_result <- iv_probit_linear_regression(
  data = data,
  outcome_var = "nEDTests",
  instrument_var = "batch.tendency",
  control_vars = c("LAB_PERF", "dayofweekt", "complaint_esi", "month_of_year", 'test.inclination'),
  cluster_var = "ED_PROVIDER"
)
nEDTests_result

RTN_72_HR_ADMIT_result <- iv_probit_linear_regression(
  data = data,
  outcome_var = "RTN_72_HR_ADMIT",
  instrument_var = "batch.tendency",
  control_vars = c("LAB_PERF", "dayofweekt", "complaint_esi", "month_of_year", 'test.inclination'),
  cluster_var = "ED_PROVIDER"
)



library(boot)
library(sandwich)
library(lmtest)

iv_probit_mediation_clustered <- function(data, outcome_var, mediator_var, instrument_var, control_vars, cluster_var, n_bootstrap = 1000) {
  
  # Helper function for bootstrapping
  boot_function <- function(data, indices) {
    # Create bootstrap sample respecting clusters
    cluster_ids <- unique(data[[cluster_var]])
    sampled_clusters <- sample(cluster_ids, replace = TRUE)
    boot_data <- do.call(rbind, lapply(sampled_clusters, function(id) {
      data[data[[cluster_var]] == id, ]
    }))
    
    # First stage: Probit regression
    first_stage_formula <- as.formula(paste("batched ~", instrument_var, "+", paste(control_vars, collapse = " + ")))
    first_stage <- glm(first_stage_formula, family = binomial(link = "probit"), data = boot_data)
    
    # Generate predicted probabilities
    boot_data$predicted_batched_prob <- predict(first_stage, newdata = boot_data, type = "response")
    
    # Mediator model
    mediator_formula <- as.formula(paste(mediator_var, "~ predicted_batched_prob +", paste(control_vars, collapse = " + ")))
    mediator_model <- lm(mediator_formula, data = boot_data)
    
    # Outcome model
    outcome_formula <- as.formula(paste(outcome_var, "~", mediator_var, "+ predicted_batched_prob +", paste(control_vars, collapse = " + ")))
    outcome_model <- lm(outcome_formula, data = boot_data)
    
    # Calculate clustered standard errors
    clustered_se_mediator <- vcovCL(mediator_model, cluster = boot_data[[cluster_var]], type = "HC1")
    clustered_se_outcome <- vcovCL(outcome_model, cluster = boot_data[[cluster_var]], type = "HC1")
    
    # Calculate effects
    a <- coef(mediator_model)["predicted_batched_prob"]
    b <- coef(outcome_model)[mediator_var]
    c_prime <- coef(outcome_model)["predicted_batched_prob"]
    
    # Calculate standard errors
    se_a <- sqrt(diag(clustered_se_mediator))["predicted_batched_prob"]
    se_b <- sqrt(diag(clustered_se_outcome))[mediator_var]
    se_c_prime <- sqrt(diag(clustered_se_outcome))["predicted_batched_prob"]
    
    print('Round complete')
    # Return effects and standard errors
    return(c(a = a, b = b, c_prime = c_prime, ab = a*b, se_a = se_a, se_b = se_b, se_c_prime = se_c_prime))
  }
  
  # Perform bootstrapping
  boot_results <- boot(data = data, statistic = boot_function, R = n_bootstrap)
  
  # Calculate confidence intervals
  ci_ab <- boot.ci(boot_results, type = "perc", index = 4)  # Indirect effect
  ci_c_prime <- boot.ci(boot_results, type = "perc", index = 3)  # Direct effect
  
  # Calculate z-scores and p-values
  z_indirect <- mean(boot_results$t[,4]) / mean(boot_results$t[,5] * boot_results$t[,6])  # Using delta method approximation
  z_direct <- mean(boot_results$t[,3]) / mean(boot_results$t[,7])
  p_indirect <- 2 * (1 - pnorm(abs(z_indirect)))
  p_direct <- 2 * (1 - pnorm(abs(z_direct)))
  
  # Prepare results
  results <- list(
    indirect_effect = mean(boot_results$t[,4]),
    direct_effect = mean(boot_results$t[,3]),
    total_effect = mean(boot_results$t[,3]) + mean(boot_results$t[,4]),
    ci_indirect = ci_ab$percent[4:5],
    ci_direct = ci_c_prime$percent[4:5],
    z_indirect = z_indirect,
    z_direct = z_direct,
    p_indirect = p_indirect,
    p_direct = p_direct
  )
  
  return(results)
}


mediation_results <- iv_probit_mediation_clustered(
  data = data,
  outcome_var = "ln_ED_LOS",
  mediator_var = "nEDTests",
  instrument_var = "batch.tendency",
  control_vars = c("LAB_PERF", "dayofweekt", "complaint_esi", "month_of_year", 'test.inclination'),
  cluster_var = "ED_PROVIDER",
  n_bootstrap = 10
)

print(mediation_results)




mediation_results <- iv_probit_mediation_clustered(
  data = data,
  outcome_var = "RTN_72_HR_ADMIT",
  mediator_var = "nEDTests",
  instrument_var = "batch.tendency",
  control_vars = c("LAB_PERF", "dayofweekt", "complaint_esi", "month_of_year", 'test.inclination'),
  cluster_var = "ED_PROVIDER",
  n_bootstrap = 10
)

print(mediation_results)












instrument_var = "batch.tendency"
control_vars = c("LAB_PERF", "dayofweekt", "complaint_esi", "month_of_year", "test.inclination")
# First stage: Probit regression
first_stage_formula <- as.formula(paste("batched ~", instrument_var, "+", paste(control_vars, collapse = " + ")))
first_stage <- glm(first_stage_formula, family = binomial(link = "probit"), data = data)
predict(first_stage, newdata = data, type = "response")


# Generate predicted probabilities
data$predicted_batched_prob <- predict(first_stage, newdata = data, type = "response")



summary(felm(
  ln_ED_LOS ~ 0 | LAB_PERF + dayofweekt + month_of_year + complaint_esi | 
    (batched ~ batch.tendency)| ED_PROVIDER, 
  data = data))



first_stage <- felm(batched ~ batch.tendency | LAB_PERF + dayofweekt + month_of_year + complaint_esi | 0 | ED_PROVIDER, data = data)
summary(first_stage)

reduced_form <- felm(ln_ED_LOS ~ batch.tendency | LAB_PERF + dayofweekt + month_of_year + complaint_esi | 0 | ED_PROVIDER, data = data)
summary(reduced_form)

iv_no_fe <- felm(ln_ED_LOS ~ LAB_PERF + dayofweekt + month_of_year + complaint_esi | 0 | (batched ~ batch.tendency) | ED_PROVIDER, data = data)
summary(iv_no_fe)

library(AER)
iv_model <- ivreg(ln_ED_LOS ~ batched + LAB_PERF + dayofweekt + month_of_year + complaint_esi | batch.tendency + LAB_PERF + dayofweekt + month_of_year + complaint_esi, data = data)
summary(iv_model, diagnostics = TRUE)


library(sandwich)
library(lmtest)

coeftest(iv_model, vcov = vcovHC(iv_model, type = "HC1"))

# Main Results
model.IV.LOS1 <- felm(
  ln_ED_LOS ~ 0 | LAB_PERF + dayofweekt + month_of_year + complaint_esi | 
              (batched ~ batch.tendency)| ED_PROVIDER, 
  data = data)

summary(model.IV.LOS1)

model.IV.ntest1 <- felm(
  imgTests ~ test.inclination | LAB_PERF +  dayofweekt + complaint_esi + month_of_year | 
    (batched ~ batch.tendency)| ED_PROVIDER, 
  data = data)

plot(data$test.inclination, data$batch.tendency)


model.IV.721 <- felm(
  RTN_72_HR_ADMIT ~ 0 | LAB_PERF + dayofweekt + complaint_esi + month_of_year | 
    (batched ~ batch.tendency)| ED_PROVIDER, 
  data = data)


library(AER)

iv_model <- ivreg(ln_ED_LOS ~ batched + LAB_PERF + dayofweekt + complaint_esi + 
                    month_of_year | 
                    batch.tendency + LAB_PERF + dayofweekt + complaint_esi + 
                    month_of_year, 
                  data = data)

# Summary of the IV regression
summary(iv_model)

# Weak instruments test
summary(iv_model, diagnostics = TRUE)











mod <- lm(batched ~ batch.tendency + LAB_PERF + dayofweekt + complaint_esi + 
     month_of_year, data = data)

data$batch.p <- predict(mod, newdata = data, type = "response")

mod.2 <- lm(ln_ED_LOS ~ batch.p + LAB_PERF + dayofweekt + complaint_esi + 
     month_of_year, data = data)

stargazer(mod, mod.2, type = "text", title = "IV Results", column.labels = c("LOS", "nTest"), style = 'QJE',
          omit = c('month_of_year', 'dayofweekt', 'complaint_esi', 'LAB_PERF'))

stargazer(model.IV.LOS1, type = "text", title = "IV Results", column.labels = c("LOS"), style = 'QJE')


summary(mod.2)


stargazer(model.IV.LOS1, model.IV.ntest1, model.IV.721, 
          type = "text", 
          title = "IV Results",
          column.labels = c("LOS", "nTest", "RTN"),
          style = 'QJE')



model.IV.LOS2 <- felm(
  ln_ED_LOS ~ 0 | LAB_PERF + dayofweekt + month_of_year + complaint_esi | 
    (batched ~ batch.tendency.2)| ED_PROVIDER, 
  data = data)

model.IV.ntest2 <- felm(
  nEDTests ~ 0  | LAB_PERF +  dayofweekt + complaint_esi + month_of_year | 
    (batched ~ batch.tendency.2)| ED_PROVIDER, 
  data = data)

model.IV.722 <- felm(
  RTN_72_HR_ADMIT ~ 0 | LAB_PERF + dayofweekt + complaint_esi + month_of_year | 
    (batched ~ batch.tendency.2)| ED_PROVIDER, 
  data = data)




stargazer(model.IV.LOS2, model.IV.ntest2, model.IV.722, 
          type = "text", 
          title = "IV Results",
          column.labels = c("LOS", "nTest", "RTN"),
          style = 'QJE')


# do 2SLS by hand
batch <- glm(batched ~ batch.tendency + test.inclination + LAB_PERF + dayofweekt + 
               month_of_year + complaint_esi, data = data, family = 'binomial')

data$batched.hat1 <-  predict(batch, newdata = data, type = "response")

batch <- glm(batched ~ batch.tendency.2 + LAB_PERF + dayofweekt + 
               month_of_year + complaint_esi, data = data, family = 'binomial')

data$batched.hat2 <-  predict(batch, newdata = data, type = "response")

return1 <- felm(RTN_72_HR_ADMIT ~ batched.hat + test.inclination | LAB_PERF + dayofweekt + 
                month_of_year + complaint_esi|0| ED_PROVIDER, data = data)

ntest1 <- felm(nEDTests ~ batched.hat + test.inclination | LAB_PERF + dayofweekt + 
                 month_of_year + complaint_esi|0| ED_PROVIDER, data = data)

los1 <- felm(ln_ED_LOS ~ batched.hat + test.inclination | LAB_PERF + dayofweekt + 
                 month_of_year + complaint_esi|0| ED_PROVIDER, data = data)



return2 <- felm(RTN_72_HR_ADMIT ~ batched.hat2 | LAB_PERF + dayofweekt + 
                  month_of_year + complaint_esi|0| ED_PROVIDER, data = data)

ntest2 <- felm(nEDTests ~ batched.hat2 | LAB_PERF + dayofweekt + 
                 month_of_year + complaint_esi|0| ED_PROVIDER, data = data)

los2 <- felm(ln_ED_LOS ~ batched.hat2 | LAB_PERF + dayofweekt + 
               month_of_year + complaint_esi|0| ED_PROVIDER, data = data)


stargazer(return1, ntest1, los1, type = "text", 
          title = "2SLS Results: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')


stargazer(return2, ntest2, los2, type = "text", 
          title = "2SLS Results: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')






model.los <- felm(
  ln_ED_LOS ~ test.inclination + LAB_PERF | complaint_esi + month_of_year | 
    (batched ~ batch.tendency)| ED_PROVIDER, 
  data = data)

model.72 <- felm(
  RTN_72_HR_ADMIT ~ test.inclination + LAB_PERF | complaint_esi | 
    (batched ~ batch.tendency)| ED_PROVIDER, 
  data = data)


summary(model.72)

summary(model.IV.ntest)

library(AER)
iv2 = ivreg(RTN_72_HR ~ batched + LAB_PERF + dayofweekt | 
              batch.tendency + test.inclination + LAB_PERF + dayofweekt, data = data)

summary(iv2, vcov = sandwich, diagnostics = TRUE)

summary(iv2)


# do 2SLS by hand
batch <- lm(batched ~ batch.tendency + LAB_PERF + dayofweekt + month_of_year + complaint_esi, data = data)
data$batched.hat <- predict(batch)
ntests <- lm(ln_ED_LOS ~ batched.hat + LAB_PERF + dayofweekt + month_of_year + complaint_esi, data = data)
summary(ntests)




summary(felm(
  ln_ED_LOS ~ LAB_PERF + test.inclination | dayofweekt + month_of_year + complaint_esi | 
   (batched ~ batch.tendency + LAB_PERF + test.inclination) | ED_PROVIDER , data = data
))

library(AER)
summary(ivreg(ln_ED_LOS ~ batched + LAB_PERF + dayofweekt + month_of_year + complaint_esi |
                  batch.tendency + LAB_PERF + dayofweekt + month_of_year + complaint_esi,
      data = data))


summary(model.IV.LOS)
summary(model.IV.72)

stargazer(list(model.IV.72,
               model.IV.LOS,
               model.IV.ntest), type = "text", 
          title = "2SLS Results: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')

# first stage
fs <- lm(
  batched ~ batch.tendency + test.inclination + LAB_PERF + dayofweekt + month_of_year + complaint_esi, 
  data = data)

data$predicted_batching <- predict(fs)

outcome_model <- lm(ln_ED_LOS ~ predicted_batching + nEDTests + LAB_PERF + test.inclination +
                    dayofweekt + month_of_year + complaint_esi, data = data)

summary(outcome_model)

summary(model.IV.LOS)

mediator_model <- lm(nEDTests ~ predicted_batching + test.inclination + LAB_PERF +
                     dayofweekt + month_of_year + complaint_esi, data = data)


med.out.2 <- mediation::mediate(
  mediator_model, outcome_model, treat = "predicted_batching", 
  mediator = "nEDTests", data = data[-c(4664, 24155,45554), ], n.boot = 100)

summary(med.out.1)
summary(med.out.2)

# Save the results to a .txt file
sink("outputs/tables/2SLS Results.txt")

stargazer(list(model.IV.72,
               model.IV.LOS,
               model.IV.ntest), type = "text", 
          header = F, 
          title = "2SLS Results: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')

stargazer(list(model.IV.LOS,
               model.IV.ntest), type = "latex", 
          covariate.labels = c('Testing Inclination', 'Batch'),
          dep.var.labels = c('Log LOS', 'Number of Tests', '72 Hour Return'),
          header = F, 
          title = "2SLS Results: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')

med.out.1$call
summary(med.out.1)

med.out.2$call
summary(med.out.2)


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
   group_by(batch, sequenced, batched) %>%
   summarize(n = n()) %>% view()


# Save the results to a .txt file
sink("outputs/tables/2SLS Results Complaints.txt")

subset <- filter(data, GROUPED_COMPLAINT == unique_complaints[3])
# Loop over each unique complaint
i = 0
for(complaint in unique_complaints) {
  
  i <- i + 1
  subset <- filter(data, GROUPED_COMPLAINT == unique_complaints[i])
  # Run your regression models on this subset
  model.IV.LOS <- felm(ln_ED_LOS ~ 0 | LAB_PERF + dayofweekt + month_of_year + ESI | 
                         (batched ~ batch.tendency )| ED_PROVIDER, data = subset)
  
  model.IV.ntest <- felm(nEDTests ~ 0 | LAB_PERF + dayofweekt + month_of_year + ESI | 
                           (batched ~ batch.tendency )| ED_PROVIDER, data = subset)
  
  model.IV.72 <- felm(RTN_72_HR ~ 0 | LAB_PERF + dayofweekt + month_of_year + ESI | 
                        (batched ~ batch.tendency )| ED_PROVIDER, data = subset)

  # Generate and print the stargazer table for the current complaint
  stargazer(list(model.IV.LOS, model.IV.ntest, model.IV.72), 
            type = "text", 
            dep.var.labels = c('Log LOS', 'Number of Tests', '72 Hour Return'),
            header = F, 
            title = paste("2SLS Results for", complaint, ": Length of Stay, Number of Tests, and 72-Hour"), 
            style = 'QJE')
}

sink()



m.2 <- felm(nEDTests ~ test.inclination | dayofweekt + month_of_year + CHIEF_COMPLAINT | 
               (any.batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
             data = subset(data, ESI >= 1 & ESI < 3))

m.3 <- felm(nEDTests ~ test.inclination | dayofweekt + month_of_year + CHIEF_COMPLAINT | 
               (any.batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
             data = subset(data, ESI >= 3))

star1 <- stargazer(m.2, m.3, type = "text", 
                   covariate.labels = c('Testing Inclination', 'Batch'),
                   dep.var.labels = 'Number of Tests',
                   header = F, 
                   title = "2SLS Results: Number of Tests by ESI", 
                   style = 'QJE')

m.2 <- felm(ln_ED_LOS ~ test.inclination | dayofweekt + month_of_year + CHIEF_COMPLAINT | 
              (any.batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
            data = subset(data, ESI >= 1 & ESI < 3))

m.3 <- felm(ln_ED_LOS ~ test.inclination | dayofweekt + month_of_year + CHIEF_COMPLAINT | 
              (any.batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
            data = subset(data, ESI >= 3))

star1 <- stargazer(m.2, m.3, type = "text", 
                   covariate.labels = c('Testing Inclination', 'Batch'),
                   dep.var.labels = 'LOS',
                   header = F, 
                   title = "2SLS Results: Number of Tests by ESI", 
                   style = 'QJE')





# Main Results
summary(felm(
  ln_ED_LOS ~ test.inclination + any.batch*as.factor(ESI)| dayofweekt + month_of_year + CHIEF_COMPLAINT | 
    (any.batch*as.factor(ESI) ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data))

felm(
  nEDTests ~ test.inclination + any.batch | dayofweekt + month_of_year + complaint_esi | 
    (any.batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)

 felm(
  RTN_72_HR ~ test.inclination | dayofweekt + month_of_year + complaint_esi | 
    (any.batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data)




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

#=========================================================================
# Robustness, Change definition of batching
#=========================================================================

data$any.batch[!is.na(data$sequenced)] <- 0

data$residual_batch <- resid(
  felm(any.batch ~ 0 | dayofweekt + month_of_year + complaint_esi, data=final))


data <- data %>%
  group_by(ED_PROVIDER) %>%
  mutate(Sum_Resid=sum(residual_batch, na.rm=T),
         batch.tendency = (Sum_Resid - residual_batch) / (n() - 1)) %>%
  ungroup()


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


stargazer(list(model.IV.LOS,
               model.IV.ntest,
               model.IV.72), type = "latex", 
          covariate.labels = c('Batch'),
          dep.var.labels = c('Log LOS', 'Number of Tests', '72 Hour Return',
                             'Admission'),
          header = F, 
          title = "2SLS Results: Length of Stay, Number of Tests, and 72-Hour", 
          style = 'QJE')

#=========================================================================
# Mediation Analysis
data <- read_csv('../outputs/data/final.csv')


model.first.stage <- felm(
  nEDTests ~ 0 + test.inclination | dayofweekt + month_of_year + complaint_esi | 
    (any.batch ~ batch.tendency + test.inclination) | ED_PROVIDER, 
  data = data
)

# Extract the fitted values from the first stage
data$Predicted_nEDTests <- model.first.stage$fitted.values

# Second Stage: Use the predicted 'Number of Tests' to predict 'LOS'
model.second.stage <- felm(
  ln_ED_LOS ~ Predicted_nEDTests + 0 | dayofweekt + month_of_year + complaint_esi | (any.batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data
)

# Estimate the direct effect of 'Batch' on 'LOS', controlling for the mediator 'Number of Tests'
model.direct.effect <- felm(
  ln_ED_LOS ~ any.batch + nEDTests + 0 | 
    dayofweekt + month_of_year + complaint_esi | 
    (any.batch ~ batch.tendency + test.inclination)| ED_PROVIDER, 
  data = data
)

# Output the summary of the models
summary(model.first.stage)      # Outputs the results of the first stage regression
summary(model.second.stage)     # Outputs the results of the second stage regression
summary(model.direct.effect)    # Outputs the results for the direct effect estimation

# Calculating the indirect effect:
# Coefficient for 'any.batch' from the first stage
coef_first_stage <- coef(model.first.stage)[["`any.batch(fit)`"]]

# Coefficient for 'Predicted_nEDTests' from the second stage
coef_second_stage <- coef(model.second.stage)[["Predicted_nEDTests"]]

# Indirect effect is the product of the two coefficients
indirect_effect <- coef_first_stage * coef_second_stage

# Output the indirect effect
print(paste("Indirect effect of 'Batch' on 'LOS' through 'Number of Tests':", indirect_effect))
#=========================================================================

library(lfe)

# Step 1: Model the effect of Batch on Number of Tests (the mediator)
# We use 'felm' to account for potential high-dimensional fixed effects
model_batch_to_tests <- felm(
  nEDTests ~ any.batch | dayofweekt + month_of_year + complaint_esi | 0 | ED_PROVIDER, 
  data = data
)

# Step 2: Model the effect of Number of Tests on LOS, controlling for Batch
# Again, we include fixed effects for time and complaint_esi
model_tests_to_los <- felm(
  ln_ED_LOS ~ nEDTests + any.batch | dayofweekt + month_of_year + complaint_esi | 0 | ED_PROVIDER, 
  data = data
)

# Extract the fitted values from the first model to use in the second model
data$predicted_num_tests <- model_batch_to_tests$fitted.values

# Step 3: Compute the total causal effect
# Average Causal Mediation Effect (ACME) is the product of the coefficients from
# the effect of Batch on Number of Tests and the effect of Number of Tests on LOS
acme <- coef(model_batch_to_tests)[[1]] * coef(model_tests_to_los)[[1]]

# Average Direct Effect (ADE) is the coefficient of Batch on LOS directly
ade <- coef(model_tests_to_los)[[2]]
# The total causal effect is the sum of ACME and ADE
total_causal_effect <- acme + ade

# Output the total causal effect
print(paste("Average Causal Mediation Effect (ACME):", acme))
print(paste("Average Direct Effect (ADE):", ade))
print(paste("Total Causal Effect of Batch on LOS:", total_causal_effect))


#=========================================================================
library(mediation)

med.fit <- lm(nEDTests ~ any.batch + test.inclination + 
              complaint_esi, 
              data=data) # Mediator model

out.fit <- lm(ln_ED_LOS ~ nEDTests + any.batch + test.inclination + 
              complaint_esi, 
              data=data) # Outcome model

med.out <- mediate(med.fit, out.fit, treat="any.batch", mediator="nEDTests")

summary(med.out)




# Load necessary library
library(stats)

# First Stage: Use lm to predict any.batch based on batch.tendency and other controls
first_stage <- lm(any.batch ~ batch.tendency + test.inclination  + complaint_esi, data = data)
data$predicted_any_batch <- predict(first_stage, newdata = data)

# Mediator Model: Model nEDTests as a function of predicted any.batch and other controls
mediator_model <- lm(nEDTests ~ predicted_any_batch + test.inclination  + complaint_esi, data = data)

# Outcome Model: Model ln_ED_LOS as a function of nEDTests (mediator) and predicted any.batch
outcome_model <- lm(ln_ED_LOS ~ nEDTests + predicted_any_batch + complaint_esi, data = data)


med.out <- mediate(mediator_model, outcome_model, 
                   treat="predicted_any_batch", 
                   mediator="nEDTests", sims=1000)

summary(med.out)




dat <- data
# first stage
first_stage <- lm(batched ~ batch.tendency + test.inclination + imaging + dayofweekt + 
                    month_of_year + LAB_PERF + complaint_esi,  dat)

dat$p_batched <- predict(first_stage, newdata = dat)

mod <- glm(imgTests ~ p_batched:complaint_esi + test.inclination + imaging + dayofweekt + 
      month_of_year + LAB_PERF, data = dat)

#cluster SE
coeftest(mod, vcov = vcovHC(mod, type = "HC1", cluster = dat$ED_PROVIDER))












