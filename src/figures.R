#=========================================================================
# Purpose: Main R file for Generating Figures and Tables
# Author: Jacob Jameson 
#=========================================================================
rm(list = ls()) 
source('src/clean.R')

library(caret)
library(lmtest)
library(sandwich)
library(ggsci)
library(xtable) # Output to LaTeX table format

##########################################################################
#=========================================================================
# Table []
#   - Show that across physician, there is random assignment and patient
#     characteristics are balanced
#=========================================================================

fig.1 <- final

chief_complaint_freq <- table(fig.1$CHIEF_COMPLAINT)

top_10_chief_complaints <- names(chief_complaint_freq)[order(chief_complaint_freq, 
                                                             decreasing = TRUE)][1:10]

fig.1$CHIEF_COMPLAINT <- ifelse(fig.1$CHIEF_COMPLAINT %in% top_10_chief_complaints, 
                                fig.1$CHIEF_COMPLAINT, "DROP")

dummy <- dummyVars(" ~ CHIEF_COMPLAINT", data = fig.1)
one.hot <- data.frame(predict(dummy, newdata = fig.1))
vars <- setdiff(names(one.hot), 'CHIEF_COMPLAINTDROP')

fig.1 <- cbind(fig.1, one.hot)

fig.1$ESI1.2 <- ifelse(fig.1$ESI == 1 | fig.1$ESI == 2, 1, 0)
fig.1$ESI3.4.5 <- ifelse(fig.1$ESI %in% c(3, 4, 5), 1, 0)

vars <- c(vars, 'ESI1.2', 'ESI3.4.5')

# Model and collect balance statistics
balance <- data.frame(Df = numeric(), 
                      F = numeric(),
                      Pr..F. = numeric(), 
                      dummy = character())

for (dummy in vars) {
  model_pos_1 <- lm(as.formula(paste(dummy, '~ 1')), fig.1)
  model_pos_2 <- lm(as.formula(paste(dummy, '~ ED_PROVIDER')), fig.1)
  
  wald_pos <- waldtest(model_pos_1, model_pos_2, 
                       vcov = vcovHC(model_pos_2, type = "HC1"))
  
  temp <- data.frame(wald_pos)[2, c(2,3,4)]
  temp$dummy <- dummy
  
  balance <- rbind(balance, temp)
}

print(xtable(balance, 
             caption = "Wald Test Results"), type = "latex")

##########################################################################

##########################################################################
#=========================================================================
# Figure []
#   - Show that frequency of batching that occurs across complaints
#=========================================================================

fig.2 <- final 

collapsed_df <- fig.2 %>%
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


temp <- fig.2 %>% 
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








