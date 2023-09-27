#=========================================================================
# Purpose: Construct main IV and Validate
# Author: Jacob Jameson 
#=========================================================================

source('src/clean.R')
library(lfe)
library(stargazer)
library(texreg)
library(xtable)

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

##########################################################################
#=========================================================================
# Exclusion --------------------------------------------------------------
## Placebo Check
#=========================================================================
##########################################################################

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

