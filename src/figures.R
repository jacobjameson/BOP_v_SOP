#=========================================================================
# Purpose: Main R file for Generating Figures and Tables
# Author: Jacob Jameson 
#=========================================================================
rm(list = ls()) 

library(caret)
library(lmtest)
library(sandwich)
library(ggsci)
library(xtable) # Output to LaTeX table format

data <- read_csv('outputs/data/final.csv')


##########################################################################


library(ggplot2)
library(lfe)


# Run your fixed effects model
FS <- lm(any.batch ~ batch.tendency + as.factor(dayofweekt) + 
           as.factor(month_of_year), data)

data$predlm <- predict(FS)
predslm <- predict(FS, interval = "confidence")
data <- cbind(data, predslm)

predslm
# Extract the fitted values
data$fitted_values <- first_stage_final$fitted.values
data$predlm = predict(first_stage_final)


# Create the histogram with the regression curve
ggplot(data, aes(x = batch.tendency)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.05, fill = "blue", alpha = 0.7) +
  geom_ribbon( aes(ymin = lwr, ymax = upr, color = NULL), alpha = .15) +
  geom_line( aes(y = fit), size = 1) +
  labs(x = "Physician Leniency", y = "Density (Histogram) / Fitted Probability (Curve)") +
  theme_minimal() +
  theme(axis.title.y.right = element_text(angle = 0)) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Fitted Probability of Prescribed Opioids"))

# Note: Adjust binwidth in geom_histogram and span in geom_smooth as needed




##########################################################################
#=========================================================================
# Figure []
#   - Show that frequency of batching that occurs across complaints
#=========================================================================

data_for_plot <- final 

collapsed_df <- data_for_plot %>%
  group_by(CHIEF_COMPLAINT) %>%
  summarise(total_cases = n()) %>%
  mutate(freq = total_cases / sum(total_cases)) %>%
  bind_rows(
    final %>%
      pivot_longer(cols = c(lab_image_batch, image_image_batch),
                   names_to = "type",
                   values_to = "value") %>%
      group_by(CHIEF_COMPLAINT, type) %>%
      summarise(freq = mean(value), .groups = "drop")
  ) %>%
  arrange(CHIEF_COMPLAINT, type) %>%
  mutate(type = case_when(
    type == 'image_image_batch' ~ "Image + Image Batching Rate",
    type == 'lab_image_batch' ~ "Lab + Image Batching Rate"))


temp <- data_for_plot %>% 
  group_by(CHIEF_COMPLAINT) %>%
  summarise(total_cases = n()) %>%
  mutate(freq = total_cases / sum(total_cases),
         type = "Fraction of Complaint") %>%
  arrange(desc(freq)) %>%
  head(10)

collapsed_df <- rbind(merge(collapsed_df, select(temp, CHIEF_COMPLAINT)), temp)
collapsed_df$type <-  ifelse(is.na(collapsed_df$type) == T, 
                             "Fraction of Complaint", collapsed_df$type)


collapsed_df$type <- factor(collapsed_df$type, 
                            levels = c("Fraction of Complaint",
                                       "Lab + Image Batching Rate",
                                       "Image + Image Batching Rate"))

collapsed_df <- unique(collapsed_df)

collapsed_df %>%
  filter(type != "Fraction of Complaint") %>%
  ggplot(., aes(x = reorder(CHIEF_COMPLAINT, -total_cases), y = freq, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Chief Complaint\n",
       y = '\nFrequency',
       title = 'Frequency of Batched Cases by Complaint') +
  scale_fill_manual(values = c("Image + Image Batching Rate" = "#0E65A3",
                               "Lab + Image Batching Rate" = "#F79500")) +
  theme_minimal() +  scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(color = "black", size = 11, angle=0),
        axis.text.y = element_text(color = "black", size = 18),
        axis.title =element_blank(),
        plot.title = element_text(color = "black", size = 20, face = "bold", hjust=0.5),
        legend.title = element_blank(),
        legend.position = 'top',
        legend.text = element_text(color = "black", size = 14),
        plot.caption =element_text(color = "black", size = 14, hjust=0.5))

ggsave("manuscript/figures/frequency_batches.pdf", width = 14, height = 8)
ggsave("manuscript/figures/frequency_batches.png", width = 14, height = 8, bg = 'white')

##########################################################################

##########################################################################
#=========================================================================
# Figure []
#   - Show that frequency of batching that occurs across complaints
#=========================================================================

data_for_plot <- final

selected_complaints <- c('Upper Respiratory Symptoms',
                         'Abdominal Complaints',
                         'Back or Flank Pain',
                         'Gastrointestinal Issues')

data_for_plot$CHIEF_COMPLAINT <- ifelse(data_for_plot$CHIEF_COMPLAINT %in% selected_complaints, 
                                        data_for_plot$CHIEF_COMPLAINT, 
                                        "DROP")

complaint_data <- data_for_plot %>%
  filter(CHIEF_COMPLAINT != 'DROP') %>%
  group_by(ED_PROVIDER, CHIEF_COMPLAINT) %>%
  summarize(batch_rate = mean(any.batch), .groups = 'drop') 

# A function to check if the high and low groups meet the criteria.
check_criteria <- function(high_group, low_group, data) {
  results <- data %>%
    group_by(CHIEF_COMPLAINT) %>%
    summarize(lowest_high_group = min(batch_rate[ED_PROVIDER %in% high_group]),
              highest_low_group = max(batch_rate[ED_PROVIDER %in% low_group]))
  
  all(results$lowest_high_group > results$highest_low_group)
}

# Initial group selection based on median batch rates
provider_medians <- complaint_data %>%
  group_by(ED_PROVIDER) %>%
  summarize(median_rate = median(batch_rate)) %>%
  arrange(-median_rate)

# Initial top 4 and bottom 4 providers
high_group_providers <- provider_medians$ED_PROVIDER[1:4]
low_group_providers <- provider_medians$ED_PROVIDER[(nrow(provider_medians) - 3):nrow(provider_medians)]


while(!check_criteria(high_group_providers, low_group_providers, complaint_data)) {
  # Remove the last of the high group and the first of the low group
  high_group_providers <- high_group_providers[-4]
  low_group_providers <- low_group_providers[-1]
  
  # Add the next providers to each group
  next_high <- setdiff(provider_medians$ED_PROVIDER, high_group_providers)[1]
  next_low <- setdiff(provider_medians$ED_PROVIDER, low_group_providers)[1]
  
  high_group_providers <- c(high_group_providers, next_high)
  low_group_providers <- c(low_group_providers, next_low)
}

# Now high_group_providers and low_group_providers should meet the criteria
complaint_data <- complaint_data %>%
  mutate(category = case_when(
    ED_PROVIDER %in% high_group_providers ~ "High Group",
    ED_PROVIDER %in% low_group_providers ~ "Low Group",
    TRUE ~ "Vary"
  ))

# Plot

complaint_data %>%
  ggplot(aes(x = ED_PROVIDER, y = batch_rate, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~CHIEF_COMPLAINT, nrow=1) +
  scale_fill_manual(values = c("Vary" = "grey50",
                               "High Group" = "#F79500",
                               "Low Group" = "#0E65A3")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = '',
       y = 'Batch-Ordering Frequency\n\n',
       fill = '') +
  theme_classic() +
  theme(
    axis.text.y  = element_text(size=14, color='black'), 
    axis.text.x  = element_blank(),
    panel.grid.major = element_line(color = 'grey85', size = 0.3),
    legend.position = 'top',
    axis.title.y = element_text(color = 'black', size = 14),
    strip.text.x = element_text(color = 'black', size = 12, face = "bold"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16, face = 'bold')
  )

##########################################################################

##########################################################################
#=========================================================================
# Figure []
#   - Determine when it is optimal to batch or sequence
#=========================================================================

data_for_plot <- final %>%
  mutate(crowdedness_category = cut(overlap_per_min, 
                                    breaks = quantile(overlap_per_min, 
                                                      probs = seq(0, 1, by = 0.2), 
                                                      na.rm = TRUE),
                                    labels = c("Very Low", "Low", "Medium", 
                                               "High", "Very High"),
                                    include.lowest = TRUE)) %>%
  filter(ESI != 5, nEDTests >0)



# Get unique combinations of ESI and crowdedness_category
unique_combinations <- distinct(data_for_plot, ESI, crowdedness_category)

# List of outcomes
outcomes <- c("nEDTests", "ln_ED_LOS", "RTN_72_HR")

# Initialize an empty dataframe to store results
results_df <- data.frame(Outcome = character(),
                         ESI = character(), 
                         crowdedness_category = character(),
                         coef_any_batch = numeric(),
                         p_value_any_batch = numeric())

# Loop through each outcome
for(outcome in outcomes) {
  
  # Loop through each unique combination and run the regression
  for(i in 1:nrow(unique_combinations)) {
    
    current_ESI <- unique_combinations$ESI[i]
    current_crowdedness <- unique_combinations$crowdedness_category[i]
    
    subset_data <- subset(data_for_plot, ESI == current_ESI & crowdedness_category == current_crowdedness)
    
    model_formula <- as.formula(paste(outcome, "~ any.batch + CHIEF_COMPLAINT + age_groups"))
    model <- lm(model_formula, data = subset_data)
    coef_summary <- summary(model)$coefficients
    
    current_results <- data.frame(Outcome = outcome,
                                  ESI = current_ESI, 
                                  crowdedness_category = current_crowdedness,
                                  coef_any_batch = coef_summary["any.batch", "Estimate"],
                                  p_value_any_batch = coef_summary["any.batch", "Pr(>|t|)"]) 
    
    results_df <- rbind(results_df, current_results)
  }
}

results_df %>%
  mutate(batch = ifelse(coef_any_batch < 0, 'Batch', 'Sequence')) %>%
  ggplot(aes(y=ESI, x=crowdedness_category, fill=batch)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  facet_wrap(~Outcome, nrow=1, labeller = as_labeller(c(nEDTests = "Number of Tests",
                                                        ln_ED_LOS = "Length of Stay",
                                                        RTN_72_HR = "72 Hour Return"))) +
  scale_fill_manual(values = c("Batch" = "#F79500",
                               "Sequence" = "#0E65A3"), 
                    name = "Optimal Testing Strategy") +
  theme_minimal() + coord_equal() +
  theme(
    axis.text.y  = element_text(size=16, color='black'), 
    axis.text.x  = element_text(size=16, color='black'), 
    panel.grid.major = element_line(color = 'grey85', size = 0.3),
    legend.position = 'bottom',
    axis.title.y = element_text(color = 'black', size = 16),
    axis.title.x = element_text(color = 'black', size = 16),
    strip.text.x = element_text(color = 'black', size = 14, face = "bold"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16, face = 'bold')
  ) +
  labs(
    x = "\nCrowdedness",
    y = "Emergency Severity Index\n",
  )
  

ggsave("manuscript/figures/decision.pdf", width = 14, height = 8)
ggsave("manuscript/figures/decision.png", width = 14, height = 8, bg = 'white')



  
