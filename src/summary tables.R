#=========================================================================
# Purpose: Produce Summary Tables for Manuscript
# Author: Jacob Jameson 
#=========================================================================

library(lfe)
library(stargazer)
library(texreg)
library(xtable)
library(tidyverse)
library(caret)
library(sandwich)
library(lmtest)

data <- read_csv('outputs/data/all_clean.csv')

##########################################################################
##########################################################################
#=========================================================================
# Summary Table: Table 1
#=========================================================================
##########################################################################

# Adjust the dataset and calculate the time differences
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


# Identify the five most frequent chief complaints
top_chief_complaints <- data %>%
  count(CHIEF_COMPLAINT, sort = TRUE) %>%
  top_n(5, wt = n) %>%
  pull(CHIEF_COMPLAINT)

# Calculate the mean proportion for each of the top chief complaints
means_chief_complaints <- sapply(
  top_chief_complaints, 
  function(cc) mean(data$CHIEF_COMPLAINT == cc, na.rm = TRUE))

# Panel 1: Emergency department characteristics ====================================
panel1 <- data %>%
  summarise(
    Variable = c(
      'Total Patients', 'Patients Admitted', 'Average ED LOS (min)', 'Patients Revisited within 72 Hours',
      'Patients with IV Fluids', 'Patients with IV Meds', 
      'Patients in ED', "Tachycardic", "Tachypneic", 
      "Febrile", "Hypotensive", "ESI", 
      paste("Complaint: ", top_chief_complaints),
      "Time from Arrival to Triage (mins)", "Time from Triage to First Contact (mins)"
    ),
    Mean = c(
      n(), mean(admit, na.rm = TRUE), mean(ED_LOS, na.rm = TRUE), mean(RTN_72_HR, na.rm = TRUE),
      mean(IV_FLUIDS == "Y", na.rm = TRUE), mean(IV_MEDS == "Y", na.rm = TRUE),
      mean(patients_in_hospital, na.rm = TRUE), 
      mean(tachycardic, na.rm = TRUE), 
      mean(tachypneic, na.rm = TRUE),
      mean(febrile, na.rm = TRUE), mean(hypotensive, na.rm = TRUE),
      mean(ESI, na.rm = TRUE), means_chief_complaints,
      mean(Time_arrival_to_triage, na.rm = TRUE), mean(Time_triage_to_firstcontact, na.rm = TRUE)
    ),
    Q1 = c(rep(NA, 16), NA, quantile(data$Time_arrival_to_triage, 0.25, na.rm = TRUE), quantile(data$Time_triage_to_firstcontact, 0.25, na.rm = TRUE)),
    Median = c(rep(NA, 16), NA, median(data$Time_arrival_to_triage, na.rm = TRUE), median(data$Time_triage_to_firstcontact, na.rm = TRUE)),
    Q3 = c(rep(NA, 16), NA, quantile(data$Time_arrival_to_triage, 0.75, na.rm = TRUE), quantile(data$Time_triage_to_firstcontact, 0.75, na.rm = TRUE))
  ) %>%
  mutate(Category = "Emergency Department Characteristics")

# Panel 2: Patient demographics ==================================================
panel2 <- data %>%
  summarise(
    Variable = c("Percent Male", "Race: White", "Race: Black", "Race: Asian", 
                 "Gender: Female", "Arrival Age"),
    Mean = c(
      mean(GENDER == "Male", na.rm = TRUE),
      mean(race == "white", na.rm = TRUE), 
      mean(race == "black", na.rm = TRUE), 
      mean(race == "asian", na.rm = TRUE),
      mean(GENDER == "Female", na.rm = TRUE),
      mean(ARRIVAL_AGE_DI, na.rm = TRUE)
    ),
    Q1 = c(rep(NA, 5), quantile(ARRIVAL_AGE_DI, 0.25, na.rm = TRUE)),
    Median = c(rep(NA, 5), median(ARRIVAL_AGE_DI, na.rm = TRUE)),
    Q3 = c(rep(NA, 5), quantile(ARRIVAL_AGE_DI, 0.75, na.rm = TRUE))
  ) %>%
  mutate(Category = "Patient Demographics")

# Panel 3: Diagnostic tests and outcomes ==========================================
panel3 <- data %>%
  summarise(
    Variable = c(
      "X-Rays Performed", "Ultrasounds Performed", "CTs Performed", "Labs Ordered", 
      "Patients Discharged", "Patients Admitted", "Contrast CT Performed",
      "Time to Result: X-Ray (mins)", "Time to Result: Ultrasound (mins)", 
      "Time to Result: Contrast CT (mins)", "Time to Result: Non-Contrast CT (mins)", "Time to Result: Lab (mins)"
    ),
    Mean = c(
      mean(XR_PERF > 0, na.rm = TRUE), mean(US_PERF > 0, na.rm = TRUE),
      mean(NON_CON_CT_PERF > 0 | CON_CT_PERF > 0, na.rm = TRUE),
      mean(LAB_PERF > 0, na.rm = TRUE), 
      mean(ED_DISPOSITION == "Discharge", na.rm = TRUE),
      mean(ED_DISPOSITION == "Admit", na.rm = TRUE),
      mean(CON_CT_PERF > 0, na.rm = TRUE),
      mean(Time_to_Result_XRay, na.rm = TRUE), mean(Time_to_Result_Ultrasound, na.rm = TRUE),
      mean(Time_to_Result_CTcon, na.rm = TRUE), mean(Time_to_Result_CTnon, na.rm = TRUE),
      mean(Time_to_Result_Lab, na.rm = TRUE)
    ),
    Q1 = c(rep(NA, 7), quantile(data$Time_to_Result_XRay, 0.25, na.rm = TRUE), quantile(data$Time_to_Result_Ultrasound, 0.25, na.rm = TRUE), quantile(data$Time_to_Result_CTcon, 0.25, na.rm = TRUE), quantile(data$Time_to_Result_CTnon, 0.25, na.rm = TRUE), quantile(data$Time_to_Result_Lab, 0.25, na.rm = TRUE)),
    Median = c(rep(NA, 7), median(data$Time_to_Result_XRay, na.rm = TRUE), median(data$Time_to_Result_Ultrasound, na.rm = TRUE), median(data$Time_to_Result_CTcon, na.rm = TRUE), median(data$Time_to_Result_CTnon, na.rm = TRUE), median(data$Time_to_Result_Lab, na.rm = TRUE)),
    Q3 = c(rep(NA, 7), quantile(data$Time_to_Result_XRay, 0.75, na.rm = TRUE), quantile(data$Time_to_Result_Ultrasound, 0.75, na.rm = TRUE), quantile(data$Time_to_Result_CTcon, 0.75, na.rm = TRUE), quantile(data$Time_to_Result_CTnon, 0.75, na.rm = TRUE), quantile(data$Time_to_Result_Lab, 0.75, na.rm = TRUE))
  ) %>%
  mutate(Category = "Diagnostic Tests and Outcomes")

# Combine panels into one table
summary_table <- bind_rows(panel1, panel2, panel3) %>%
  select(Category, Variable, Mean, Q1, Median, Q3)

# Export summary table to a text file
write.table(summary_table, "outputs/tables/Summary Statistics.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

##########################################################################

patient_data <- data %>%
  filter(imaging == 1)

# One-hot encoding for CHIEF_COMPLAINT
dummy_encoder <- dummyVars(" ~ CHIEF_COMPLAINT", data = patient_data)
one_hot_encoded_data <- data.frame(predict(dummy_encoder, newdata = patient_data))
relevant_vars <- setdiff(names(one_hot_encoded_data), 'CHIEF_COMPLAINTDROP')

# Add the one-hot encoded data to the main dataset
patient_data <- cbind(patient_data, one_hot_encoded_data)

# Create binary flags for ESI values
patient_data <- patient_data %>%
  mutate(ESI1 = as.integer(ESI == 1),
         ESI2 = as.integer(ESI == 2),
         ESI3 = as.integer(ESI == 3),
         ESI4 = as.integer(ESI == 4),
         ESI5 = as.integer(ESI == 5))

# Update the relevant variables list
relevant_vars <- c(relevant_vars, 'ESI1', 'ESI2', 'ESI3', 'ESI4', 'ESI5', 
                   "tachycardic", "tachypneic","febrile","hypotensive")


# Initialize balance dataframe
balance_stats <- data.frame(Df = numeric(), 
                            F = numeric(),
                            Pr_F = numeric(), 
                            dummy_var = character())

# Loop through each variable and collect balance statistics
for (dummy_var in relevant_vars) {
  baseline_model <- lm(as.formula(paste(dummy_var, '~ 1 + dayofweekt + month_of_year')), patient_data)
  extended_model <- lm(as.formula(paste(dummy_var, '~ ED_PROVIDER + dayofweekt + month_of_year')), patient_data)
  
  wald_test_result <- lmtest::waldtest(baseline_model, extended_model, 
                                       vcov = vcovHC(extended_model, type = "HC1"))
  
  temp_result <- data.frame(wald_test_result)[2, c(2,3,4)]
  temp_result$dummy_var <- dummy_var
  
  balance_stats <- rbind(balance_stats, temp_result)
}

row.names(balance_stats) <- NULL

balance_stats$Adjusted_p_value <- p.adjust(balance_stats$Pr..F., method = "bonferroni")

print(xtable(balance_stats[,c(2,3,4,5)], 
             caption = "Wald Test Results", 
             digits = c(0,3,3,0, 3)), type = "latex")


balance_stats



data %>%
  group_by(CHIEF_COMPLAINT) %>%
  summarize(avgimg = mean(imaging), 
            avgxr = mean(XR_PERF), 
            avgus = mean(US_PERF), 
            avgnoncon = mean(NON_CON_CT_PERF), 
            avgcon = mean(CON_CT_PERF), 
            avgtests = mean(imgTests)) %>%
  arrange(desc(avgtests)) 











all_tests <- c('US_ORDER_DTTM_REL', 'NON_CON_CT_ORDER_REL', 
               'CON_CT_ORDER_REL', 'XR_ORDER_REL')
library(dplyr)
library(tidyr)
library(ggalluvial)
library(ggplot2)

test_columns <- c("US_ORDER_DTTM_REL", "NON_CON_CT_ORDER_REL", "CON_CT_ORDER_REL", "XR_ORDER_REL")

# Prepare the data
prepared_data <- data %>%
  select(all_of(test_columns)) %>%
  mutate(patient_id = row_number()) %>%
  pivot_longer(cols = all_of(test_columns), 
               names_to = "test_type", 
               values_to = "order_time") %>%
  filter(!is.na(order_time)) %>%
  group_by(patient_id) %>%
  arrange(patient_id, order_time) %>%
  mutate(test_order = row_number()) %>%
  ungroup()

# Create flow data for all transitions
flow_data <- prepared_data %>%
  group_by(patient_id) %>%
  mutate(next_test = lead(test_type, default = "End"),
         next_order = lead(test_order, default = max(test_order) + 1)) %>%
  ungroup() %>%
  mutate(flow = paste(test_order, next_order, sep = "_"))

# Create the alluvial plot
ggplot(flow_data, 
       aes(x = test_order, stratum = test_type, alluvium = patient_id,
           y = 1, fill = test_type, label = test_type)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "white") +
  geom_stratum(alpha = 0.5) +
  geom_text(stat = "stratum", size = 3) +
 scale_x_continuous(breaks = 1:max(flow_data$test_order)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_minimal() +
  labs(x = "Test Sequence", y = "Frequency", fill = "Test Type") +
  ggtitle("Common Flows Between Tests") +
  theme(legend.position = "bottom")



ggplot(flow_data, 
       aes(x = test_order, stratum = test_type, alluvium = patient_id,
           y = 1, fill = test_type, label = test_type)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", 
             alpha = 0.8, width = 0.4) +
  geom_stratum(alpha = 1, width = 0.4) +
  geom_text(stat = "stratum", size = 3, color = "black") +
  scale_x_continuous(breaks = 1:max(flow_data$test_order)) +
          #           labels = c(paste("Test", 1:max(flow_data$test_order)), "End")) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  theme_minimal() +
  labs(x = "Test Sequence", y = "Frequency", fill = "Test Type") +
  ggtitle("Common Flows Between Tests") +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )



patient_test_sequence <- data %>%
  mutate(Patient_ID = row_number()) %>%
  select(Patient_ID, all_of(test_columns)) %>%
  pivot_longer(cols = all_of(test_columns), 
               names_to = "test_type", 
               values_to = "order_time") %>%
  filter(!is.na(order_time)) %>%
  mutate(test_type = case_when(
    test_type == "US_ORDER_DTTM_REL" ~ "US",
    test_type == "NON_CON_CT_ORDER_REL" ~ "CT Non-Con",
    test_type == "CON_CT_ORDER_REL" ~ "CT Con",
    test_type == "XR_ORDER_REL" ~ "Xray"
  )) %>%
  arrange(Patient_ID, order_time) %>%
  group_by(Patient_ID) %>%
  summarise(
    test_sequence = paste(test_type, collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    test_1 = sapply(strsplit(test_sequence, ", "), function(x) if(length(x) >= 1) x[1] else "None"),
    test_2 = sapply(strsplit(test_sequence, ", "), function(x) if(length(x) >= 2) x[2] else "None"),
    test_3 = sapply(strsplit(test_sequence, ", "), function(x) if(length(x) >= 3) x[3] else "None"),
    test_4 = sapply(strsplit(test_sequence, ", "), function(x) if(length(x) >= 4) x[4] else "None")
  ) %>%
  select(Patient_ID, test_1, test_2, test_3, test_4)

# Prepare the data for the alluvial diagram
flow_data <- patient_test_sequence %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = starts_with("test_"), 
               names_to = "stage", 
               values_to = "test") %>%
  filter(test != "None") %>%
  group_by(id) %>%
  mutate(next_test = lead(test, default = "End"),
         stage = factor(stage, levels = c("test_1", "test_2", "test_3", "test_4"))) %>%
  ungroup()

flow_data <- patient_test_sequence %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = starts_with("test_"), 
               names_to = "stage", 
               values_to = "test") %>%
  filter(test != "None") %>%
  group_by(id) %>%
  mutate(next_test = lead(test, default = "End"),
         stage = factor(stage, levels = c("test_1", "test_2", "test_3", "test_4"))) %>%
  ungroup() %>%
  group_by(stage, test) %>%
  mutate(freq = n()) %>%
  ungroup()

# Create the alluvial diagram
ggplot(flow_data,
       aes(x = stage, stratum = test, alluvium = id,
           y = freq, fill = test, label = test)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            alpha = 0.7) +
  geom_stratum(alpha = 0.5) +
  geom_text(stat = "stratum", size = 3, aes(label = after_stat(n))) +
  scale_x_discrete(labels = c("First Test", "Second Test", "Third Test", "Fourth Test"),
                   expand = c(0.1, 0.1)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  theme_minimal() +
  labs(x = NULL, y = "Number of Patients", 
       title = "Flow of Patients Through Test Sequences",
       subtitle = "Showing progression from first to fourth test") +
  theme(legend.position = "bottom")


patient_test_sequence %>%
  group_by(test_1, test_2, test_3, test_4) %>%
  summarise(freq = n()) %>%
  ungroup() %>%
  mutate(test_2 = ifelse(test_2 == "None", NA, test_2),
         test_3 = ifelse(test_3 == "None", NA, test_3),
         test_4 = ifelse(test_4 == "None", NA, test_4)) %>%
ggplot(aes(axis1 = test_1, axis2 = test_2, axis3 =test_3, axis4=test_4,  y = freq)) +
  geom_alluvium(aes(fill = test_1),
                curve_type = "quintic") +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete() +
  theme_void()

vaccination
