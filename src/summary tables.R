#=========================================================================
# Purpose: Produce Summary Tables for Manuscript
# Author: Jacob Jameson 
#=========================================================================

library(lfe)
library(stargazer)
library(texreg)
library(xtable)
library(tidyverse)

data <- read_csv('outputs/data/final.csv')

##########################################################################
##########################################################################
#=========================================================================
# Summary Table: Table 1
#=========================================================================
##########################################################################

# Identify the five most frequent chief complaints
top_chief_complaints <- data %>%
  count(CHIEF_COMPLAINT, sort = TRUE) %>%
  slice_max(n = 5, order_by = n) %>%
  pull(CHIEF_COMPLAINT)

# Calculate the mean for each of the top chief complaints
means_chief_complaints <- sapply(
  top_chief_complaints, 
  function(cc) mean(data$CHIEF_COMPLAINT == cc, na.rm = TRUE))

# Panel 1: Emergency department characteristics ===================================
panel1 <- data %>%
  reframe(
    Variable = c('Patients in ED', "Tachycardic", "Tachypneic", 
                 "Febrile", "Hypotensive", "ESI", 
                 paste("Complaint: ", top_chief_complaints)),
    Mean = c(mean(patients_in_hospital, na.rm = TRUE), 
             mean(tachycardic, na.rm = TRUE), 
             mean(tachypneic, na.rm = TRUE),
             mean(febrile, na.rm = TRUE), mean(hypotensive, na.rm = TRUE),
             mean(ESI, na.rm = TRUE), means_chief_complaints),
    Q1 = c(quantile(patients_in_hospital, 0.25, na.rm = TRUE),NA, NA, NA, NA, 
           quantile(ESI, 0.25, na.rm = TRUE), rep(NA, 5)),
    Median = c(median(patients_in_hospital, na.rm = TRUE), NA, NA, NA, NA, 
               median(ESI, na.rm = TRUE), rep(NA, 5)),
    Q3 = c(quantile(patients_in_hospital, 0.75, na.rm = TRUE), NA, NA, NA, NA, 
           quantile(ESI, 0.75, na.rm = TRUE), rep(NA, 5))
  ) %>%
  mutate(Category = "Emergency Department Characteristics")

# Panel 2: Patient characteristics  =================================================
panel2 <- data %>%
  summarise(
    Variable = c("Race: White", "Race: Black", "Race: Asian", 
                 "Gender: Female", "Arrival Age"),
    Mean = c(mean(race == "white", na.rm = TRUE), 
             mean(race == "black", na.rm = TRUE), 
             mean(race == "asian", na.rm = TRUE), 
             mean(GENDER == "Female", na.rm = TRUE),
             mean(ARRIVAL_AGE_DI, na.rm = TRUE)),
    Q1 = c(NA, NA, NA, NA, quantile(ARRIVAL_AGE_DI, 0.25, na.rm = TRUE)),
    Median = c(NA, NA, NA, NA, median(ARRIVAL_AGE_DI, na.rm = TRUE)),
    Q3 = c(NA, NA, NA, NA, quantile(ARRIVAL_AGE_DI, 0.75, na.rm = TRUE))
  ) %>%
  mutate(Category = "Patient Characteristics")

# Panel 3: Tests  ==================================================================
panel3 <- data %>%
  summarise(
    Variable = c("X-Ray", "Ultrasound", "Non-Contrast CT", 
                 "Contrast CT", "Lab", "Tests were Batch Ordered"),
    Mean = c(mean(XR_PERF, na.rm = TRUE), mean(US_PERF, na.rm = TRUE),
             mean(NON_CON_CT_PERF, na.rm = TRUE),
             mean(CON_CT_PERF, na.rm = TRUE),
             mean(LAB_PERF, na.rm = TRUE), 
             mean(any.batch, na.rm = TRUE)),
    Q1 = NA,
    Median = NA,
    Q3 = NA
  ) %>%
  mutate(Category = "Tests")

# Combine panels into one table  =================================================
summary_table <- bind_rows(panel1, panel2, panel3) %>%
  select(Category, Variable, Mean, Q1, Median, Q3)

# Export summary table to a text file
write.table(summary_table, "outputs/tables/Summary Statistics.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

##########################################################################
##########################################################################







