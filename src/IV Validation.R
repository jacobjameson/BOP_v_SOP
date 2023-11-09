#=========================================================================
# Purpose: Validate IV
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
first_stage_final <- felm(any.batch ~ batch.tendency + test.inclination | 
                          dayofweekt + month_of_year |
                          0| ED_PROVIDER, 
                          data = data)

# Shift-level + complaint FE
first_stage_final_d <- felm(any.batch ~ batch.tendency + test.inclination | 
                            dayofweekt + month_of_year + age_groups + 
                            complaint_esi |
                            0|ED_PROVIDER,
                            data = data)

# Shift-level + complaint + individual FE
first_stage_final_d_i <- felm(any.batch ~ batch.tendency + test.inclination | 
                              dayofweekt + month_of_year + age_groups + 
                              complaint_esi + GENDER + race |
                              0| ED_PROVIDER , 
                              data = data)

# Save the results to a .txt file
sink("outputs/tables/First Stage.txt")

screenreg(list(first_stage_final,
               first_stage_final_d,
               first_stage_final_d_i), 
          include.fstatistic = T)

stargazer(list(first_stage_final,first_stage_final_d,
               first_stage_final_d_i), type = "text", 
          header = FALSE, title = "First Stage", style = 'QJE')

sink()


##########################################################################
#=========================================================================
# Mediation Analysis -----------------------------------------------------
#=========================================================================
##########################################################################
library(mediation)

data$ESI <- as.factor(data$ESI)

reduced.form.LOS.M <- 
  lm(avg_nEDTests ~ 
         batch.tendency + dayofweekt + month_of_year + 
         age_groups + CHIEF_COMPLAINT + ESI, data = data
  )

reduced.form.LOS.Y <- 
  lm(ln_ED_LOS ~ 
       batch.tendency + avg_nEDTests + dayofweekt + 
       month_of_year + age_groups +  CHIEF_COMPLAINT + ESI, data = data
  )

mediation <- mediate(reduced.form.LOS.M, reduced.form.LOS.Y,
                     treat='batch.tendency', mediator='avg_nEDTests',
                     boot=TRUE, sims=1)

summary(mediation)

##########################################################################
#=========================================================================
# Reduced-Form Results ---------------------------------------------------
#=========================================================================
##########################################################################

reduced.form.LOS.0 <- 
  felm(ln_ED_LOS ~ 
         batch.tendency + patients_in_hospital | 
         dayofweekt + month_of_year + complaint_esi  |
         0|ED_PROVIDER, 
       data = data
  )

reduced.form.ntest.0 <- 
  felm(nEDTests ~ 
         batch.tendency + patients_in_hospital | 
         dayofweekt + month_of_year + complaint_esi |
         0|ED_PROVIDER, 
       data = data
  )

reduced.form.72.0 <- 
  felm(RTN_72_HR ~ 
         batch.tendency  + patients_in_hospital| 
         dayofweekt + month_of_year + complaint_esi  |
         0|ED_PROVIDER, 
       data = data
  )

reduced.form.LOS <- 
  felm(ln_ED_LOS ~ 
         batch.tendency + test.inclination + patients_in_hospital| 
         dayofweekt + month_of_year + complaint_esi |
         0|ED_PROVIDER, 
       data = data
   )

reduced.form.ntest <- 
  felm(nEDTests ~ 
         batch.tendency + test.inclination + patients_in_hospital | 
         dayofweekt + month_of_year + complaint_esi |
         0|ED_PROVIDER, 
       data = data
  )

reduced.form.72 <- 
  felm(RTN_72_HR ~ 
         batch.tendency + test.inclination  + patients_in_hospital| 
         dayofweekt + month_of_year + complaint_esi |
         0|ED_PROVIDER, 
       data = data
  )


data %>%
  group_by(ED_PROVIDER) %>%
  summarize(test.inclination = mean(test.inclination), 
            batch.tendency = mean(batch.tendency)) %>%
ggplot()  +
  geom_point(aes(y=test.inclination, x = batch.tendency),size=5, stroke=1) +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  geom_smooth(aes(y=test.inclination,x = batch.tendency), method = "lm",se = T) +
  theme_bw() +
  theme(plot.background=element_rect(fill='white'),
        panel.border = element_blank(),
        axis.text.y  = element_text(size=20, color='black'), 
        axis.text.x  = element_text(size=20),
        plot.margin = unit(c(0.5, 0.2, 0.2, 0.2), "cm"),
        panel.grid.major=element_line(color='grey85',size=0.3),
        legend.position = 'none',
        axis.title.y =  element_text(color = 'black',size = 20),
        axis.title.x = element_text(color = 'black',size = 20),
        strip.text.x = element_text(color = 'black', size = 20, face = "bold"),
        plot.title = element_text(color = "black", size = 30, 
                                  face = "bold", margin = margin(0,0,30,0), hjust = 0),
        plot.subtitle = element_text(color = "black", size = 14, 
                                     margin = margin(0,0,30,0),hjust = 0),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16, face = 'bold'),
        legend.background = element_rect(fill = "grey96", color = NA)) +
  labs(x='\nBatch-Ordering Tendency',
       y='Test Ordering Inclination\n')


# All coefficients are scaled by the 
# difference in batch tendency between the ninetieth and tenth lenient physicians  
# for interpretability. 
percentile_10.b <- quantile(data$batch.tendency, probs = 0.10)
percentile_90.b <- quantile(data$batch.tendency, probs = 0.90)

percentile_10.t <- quantile(data$test.inclination, probs = 0.10)
percentile_90.t <- quantile(data$test.inclination, probs = 0.90)

coeffecient.scale.b <- percentile_90.b - percentile_10.b
coeffecient.scale.t <- percentile_90.t - percentile_10.t


# Save the results to a .txt file
sink("outputs/tables/Reduced Form.txt")

stargazer(list(reduced.form.LOS.0, 
               reduced.form.ntest.0,
               reduced.form.72.0), type = "text", 
          header = FALSE, title = "Reduced Form", style = 'QJE')

stargazer(list(reduced.form.LOS, 
               reduced.form.ntest,
               reduced.form.72), type = "text", 
          header = FALSE, title = "Reduced Form", style = 'QJE')


# Scale your models

scale_factors.1 <- c(coeffecient.scale.b, coeffecient.scale.t)
scale_factors.2 <- c(coeffecient.scale.b)

scale_coefficients <- function(model, scale_factors) {

  scaled_model <- model
  
  scaled_model$coefficients[1, "Estimate"] <- model$coefficients[1, "Estimate"] * scale_factors[1]
  scaled_model$coefficients[1, "Cluster s.e."] <- model$coefficients[1, "Cluster s.e."] * scale_factors[1]
  
  if (length(scale_factors) == 2) {
    
    scaled_model$coefficients[2, "Estimate"] <- model$coefficients[2, "Estimate"] * scale_factors[2]
    scaled_model$coefficients[2, "Cluster s.e."] <- model$coefficients[2, "Cluster s.e."] * scale_factors[2]
  }
  
  return(scaled_model)
}

scale_coefficients(summary(reduced.form.LOS.0), scale_factors.1)
scale_coefficients(summary(reduced.form.ntest.0), scale_factors.1)
scale_coefficients(summary(reduced.form.72.0), scale_factors.1)

scale_coefficients(summary(reduced.form.LOS), scale_factors.2)
scale_coefficients(summary(reduced.form.ntest), scale_factors.2)
scale_coefficients(summary(reduced.form.72), scale_factors.2)

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
  summarise(mean.batch = mean(any.batch), n =n()) %>%
  filter(mean.batch <= 0.1) 

placebo.complaints <- placebo.complaints$CHIEF_COMPLAINT

placebo <- data %>%
  filter(CHIEF_COMPLAINT %in% placebo.complaints)

# The idea here is that we want to show that being a high-batching
# physician is not associated with other traits of a high-quality
# physician. We are going to show that for Urinary Complaints,
# a complaint area where batching is rare, that there are not preferable
# outcomes associated with being a batcher

#=========================================================================
# LOS
#=========================================================================

# Save the results to a .txt file
sink("outputs/tables/Placebo Check.txt")

placebo1.1 <- felm(ln_ED_LOS ~ batch.tendency + avg_nEDTests | 
                  dayofweekt + month_of_year | 0 | ED_PROVIDER, 
                  data = placebo)

placebo1.2 <- felm(ln_ED_LOS ~ batch.tendency + avg_nEDTests | 
                  dayofweekt + month_of_year + complaint_esi | 0 | ED_PROVIDER, 
                  data = placebo)

placebo1.3 <- felm(ln_ED_LOS ~ batch.tendency + avg_nEDTests | 
                  dayofweekt + month_of_year + complaint_esi + race + GENDER |0| ED_PROVIDER, 
                  data = placebo)

stargazer(list(placebo1.1, placebo1.2, placebo1.3), type = "text", 
          header = FALSE, title = "Placebo Check", style = 'QJE')

#=========================================================================
# Number of Tests
#=========================================================================

placebo2.1 <- felm(nEDTests ~ batch.tendency + avg_nEDTests | 
                   dayofweekt + month_of_year | 0 | ED_PROVIDER, 
                   data = placebo)

placebo2.2 <- felm(nEDTests ~ batch.tendency + avg_nEDTests | 
                   dayofweekt + month_of_year + complaint_esi | 0 | ED_PROVIDER, 
                   data = placebo)

placebo2.3 <- felm(nEDTests ~ batch.tendency + avg_nEDTests | 
                   dayofweekt + month_of_year + complaint_esi + race + GENDER |0| ED_PROVIDER, 
                   data = placebo)


stargazer(list(placebo2.1, placebo2.2, placebo2.3), type = "text", 
          header = FALSE, title = "Placebo Check", style = 'QJE')

#=========================================================================
# 72 HR Return
#=========================================================================

placebo3.1 <- felm(RTN_72_HR ~ batch.tendency + avg_nEDTests | 
                   dayofweekt + month_of_year | 0 | ED_PROVIDER, 
                   data = placebo)

placebo3.2 <- felm(RTN_72_HR ~ batch.tendency + avg_nEDTests | 
                   dayofweekt + month_of_year + complaint_esi | 0 | ED_PROVIDER, 
                   data = placebo)

placebo3.3 <- felm(RTN_72_HR ~ batch.tendency + avg_nEDTests | 
                   dayofweekt + month_of_year + complaint_esi + race + GENDER |0| ED_PROVIDER, 
                   data = placebo)


stargazer(list(placebo3.1, placebo3.2, placebo3.3), type = "text", 
          header = FALSE, title = "Placebo Check", style = 'QJE')

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












# Sub Group Analysis 1
## Use tendency constructed from all complaints
subgroup_analysis <- function(dependent_var, caption_text){
  
  unique_complaints <- unique(final$CHIEF_COMPLAINT)
  
  results_df <- data.frame(
    CHIEF_COMPLAINT = character(),
    Coefficient = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(complaint in unique_complaints){
    
    sub_data <- final[final$CHIEF_COMPLAINT == complaint, ]
    result <- summary(
      felm(as.formula(
        paste(dependent_var,"~ avg_nEDTests | 
                            dayofweekt + month  | 
                            (any.batch ~ batch.tendency + avg_nEDTests) | 0")),
                           data = sub_data))
    
    coeff <- result$coef["`any.batch(fit)`", "Estimate"]
    p_val <- result$coef["`any.batch(fit)`", "Pr(>|t|)"]
    
    results_df <- rbind(results_df, 
                        data.frame(CHIEF_COMPLAINT = complaint, 
                                   Coefficient = coeff, P_value = p_val))
  }
  
  # Formatting
  results_df$Coefficient <- round(results_df$Coefficient, 3)
  results_df$P_value <- round(results_df$P_value, 3)
  results_df$Significance <- ifelse(results_df$P_value < 0.01, "***", 
                                    ifelse(results_df$P_value < 0.05, "**", 
                                           ifelse(results_df$P_value < 0.1, "*", "")))
  results_df$Coefficient_Starred <- paste0(results_df$Coefficient, 
                                           results_df$Significance)
  results_df <- results_df[, c("CHIEF_COMPLAINT", "Coefficient_Starred", "P_value")]
  
  latex_table <- xtable(results_df, 
                        caption = caption_text, 
                        align = c("l","l", "r", "r"))
  
  return(print(latex_table, type = "latex", caption.placement = "top"))
}

# Run the analyses and print tables
analyses <- list(
  list(var="ln_ED_LOS", caption="LOS Regression Results by Chief Complaint"),
  list(var="nEDTests", caption="NTests Regression Results by Chief Complaint"),
  list(var="RTN_72_HR", caption="72 Return Regression Results by Chief Complaint")
)

lapply(analyses, function(x) subgroup_analysis(x$var, x$caption))

# Sub Group Analysis 2
## Use tendency constructed from each complaints

subgroup_analysis2 <- function(dependent_var, caption_text){
  
  unique_complaints <- unique(final$CHIEF_COMPLAINT)
  
  results_df <- data.frame(
    CHIEF_COMPLAINT = character(),
    Coefficient = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(complaint in unique_complaints){
    
    sub_data <- final[final$CHIEF_COMPLAINT == complaint, ]
    
    
    # Construct batch tendency for the subset
    sub_data$residual_batch <- resid(
      felm(any.batch ~ 1 | dayofweekt + month, 
           data=sub_data
      )
    )
    sub_data <- sub_data %>%
      group_by(ED_PROVIDER) %>%
      mutate(Sum_Resid=sum(residual_batch, na.rm=T),
             batch.tendency = (Sum_Resid - residual_batch) / (n() - 1)) %>% 
      ungroup()
    
    result <- summary(
      felm(as.formula(
        paste(dependent_var,"~ avg_nEDTests | 
                            dayofweekt + month  | 
                            (any.batch ~ batch.tendency + avg_nEDTests) | 0")),
        data = sub_data))
    
    coeff <- result$coef["`any.batch(fit)`", "Estimate"]
    p_val <- result$coef["`any.batch(fit)`", "Pr(>|t|)"]
    
    results_df <- rbind(results_df, 
                        data.frame(CHIEF_COMPLAINT = complaint, 
                                   Coefficient = coeff, P_value = p_val))
  }
  
  # Formatting
  results_df$Coefficient <- round(results_df$Coefficient, 3)
  results_df$P_value <- round(results_df$P_value, 3)
  results_df$Significance <- ifelse(results_df$P_value < 0.01, "***", 
                                    ifelse(results_df$P_value < 0.05, "**", 
                                           ifelse(results_df$P_value < 0.1, "*", "")))
  results_df$Coefficient_Starred <- paste0(results_df$Coefficient, 
                                           results_df$Significance)
  results_df <- results_df[, c("CHIEF_COMPLAINT", "Coefficient_Starred", "P_value")]
  
  latex_table <- xtable(results_df, 
                        caption = caption_text, 
                        align = c("l","l", "r", "r"))
  
  return(print(latex_table, type = "latex", caption.placement = "top"))
}

# Run the analyses and print tables
analyses <- list(
  list(var="ln_ED_LOS", caption="LOS Regression Results by Chief Complaint"),
  list(var="nEDTests", caption="NTests Regression Results by Chief Complaint"),
  list(var="RTN_72_HR", caption="72 Return Regression Results by Chief Complaint")
)

lapply(analyses, function(x) subgroup_analysis2(x$var, x$caption))

