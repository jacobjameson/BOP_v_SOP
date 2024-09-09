#=========================================================================
# Purpose: Main R file for Preparing/Cleaning Data
# Author: Jacob Jameson 
#=========================================================================
rm(list = ls()) 

library(tidyverse)
library(stringr)
library(lfe)
library(lubridate)

#=========================================================================
# Determine test times
#=========================================================================
path <- '~/Sue Goldie Dropbox/Jacob Jameson/Batch vs sequential testing/Data/'

df <- read.csv(paste0(path, 'deidentified_FINAL.csv'))

test_columns = c("US_PERF", "NON_CON_CT_PERF", "CON_CT_PERF", 
                 "LAB_PERF", "XR_PERF")


colnames(df)[colnames(df) == "PLAIN_XRAY"] = "XR_PERF"
colnames(df)[colnames(df) == "US_ORDER_DTTM_REL"] ="US_ORDER_DTTM_REL"
colnames(df)[colnames(df) == "CT_WITHOUT_CONTR_ORDER_DTTM_REL"] = "NON_CON_CT_ORDER_REL"
colnames(df)[colnames(df) == "CT_WITH_CONTR_ORDER_DTTM_REL"] ="CON_CT_ORDER_REL"
colnames(df)[colnames(df) == "LAB_ORDER_DTTM_REL"] ="LAB_ORDER_REL"
colnames(df)[colnames(df) == "PLAIN_XRAY_ORDER_DTTM_REL"] ="XR_ORDER_REL"

df$CT_PERF = ifelse(df$NON_CON_CT_PERF=='Y' | df$CON_CT_PERF=='Y', 1, 0)

for (i in test_columns){
  df[[i]] = ifelse(df[[i]] =='Y', 1, 0) 
}

df$nEDTests = rowSums(df[test_columns])

df$imgTests = df$nEDTests - df$LAB_PERF

# Identify columns with *_REL suffix
rel_cols <- grep("_REL$", names(df), value = TRUE)

# Apply the transformation to each *_REL column
for (col in rel_cols) {
  df <- df %>%
    separate(col, c(paste0(col, "_hours"), paste0(col, "_minutes")), sep = ":") %>%
    mutate(!!col := as.numeric(get(paste0(col, "_hours"))) * 60 + 
             as.numeric(get(paste0(col, "_minutes"))))
}

df <- df %>%
  select(-matches("_hours$|_minutes$"))


#=========================================================================
# Determine batching
#   - criteria: ordered within 5 minutes of each other
#               first tests ordered were in a batch
#=========================================================================


# Columns of interest
cols_of_interest <- c('US_ORDER_DTTM_REL', 'NON_CON_CT_ORDER_REL', 
                      'CON_CT_ORDER_REL', 'XR_ORDER_REL')

cols_res <- c("CT_WITHOUT_CONTR_RESULT_DTTM_REL",
              "CT_WITH_CONTR_RESULT_DTTM_REL", "PLAIN_XRAY_RESULT_DTTM_REL" ,     
              "US_RESULT_DTTM_REL" )

# determine earliest test order to latest test result

df$min_time <- apply(df[cols_of_interest], 1, min, na.rm = TRUE)

df$max_time <- apply(df[cols_res], 1, max, na.rm = TRUE)

df$total_testing_time <- df$max_time - df$min_time

# Function to determine batch
determine_batch <- function(row) {
  
  min_time <- min(row, na.rm = TRUE)
  tests <- names(row)
  
  batch <- tests[which(row <= min_time + 5 & !is.na(row))]
  
  return(paste(batch, collapse = ","))
}

# Apply function to each row for the columns of interest
df$batch <- apply(df[cols_of_interest], 1, determine_batch)

# Determine all tests that were ordered
df$all_tests <- apply(df[cols_of_interest], 1, function(row) {
  tests <- names(row)
  return(paste(tests[!is.na(row)], collapse = ","))
})

# get the difference between batched tests and all tests
df$diff <- apply(df, 1, function(row) {
  batch <- unlist(strsplit(as.character(row["batch"]), ","))
  all_tests <- unlist(strsplit(as.character(row["all_tests"]), ","))
  
  diff <- setdiff(all_tests, batch)
  
  if (length(diff) > 0) {
    return(paste(diff, collapse = ","))
  } else {
    return(NA)
  }
})

df$batched <- ifelse(str_count(df$batch, ",") > 0, 1, 0)

# Function to determine the sequenced tests
get_sequenced <- function(batch, row) {
  all_tests <- c('US_ORDER_DTTM_REL', 'NON_CON_CT_ORDER_REL', 
                 'CON_CT_ORDER_REL', 'XR_ORDER_REL')
  
  batch_tests <- unlist(strsplit(batch, ","))
  
  # Identify tests not in the batch
  non_batch_tests <- setdiff(all_tests, batch_tests)
  
  # Check if any of these tests have a non-NA value in the current row
  sequenced_tests <- non_batch_tests[!is.na(row[non_batch_tests])]
  
  if (length(sequenced_tests) > 0) {
    return(paste(sequenced_tests, collapse = ","))
  } else {
    return(NA)
  }
}

df$sequenced <- apply(df, 1, function(row) get_sequenced(as.character(row["batch"]), row))

#=========================================================================
# Clean Vars -------------------------------------------------------------
#=========================================================================

df$PATIENT_RACE <- str_to_lower(df$PATIENT_RACE)

df <- df %>%
  mutate(race = case_when(
    grepl('black', PATIENT_RACE, fixed = TRUE) ~ "black",
    grepl('african', PATIENT_RACE, fixed = TRUE) ~ "black",
    grepl('asian', PATIENT_RACE, fixed = TRUE) ~ "asian",
    grepl('pacific islander', PATIENT_RACE, fixed = TRUE) ~ "asian",
    grepl('native', PATIENT_RACE, fixed = TRUE)~ "native",
    grepl('samoan', PATIENT_RACE, fixed = TRUE) ~ "other",
    grepl('guamanian or chamorro', PATIENT_RACE, fixed = TRUE) ~ "other",
    grepl('white', PATIENT_RACE, fixed = TRUE) ~ "white",
    grepl('unknown', PATIENT_RACE, fixed = TRUE) ~ "unknown",
    grepl('choose not to disclose', PATIENT_RACE, fixed = TRUE) ~ "unknown",
    grepl('unable to provide', PATIENT_RACE, fixed = TRUE) ~ "unknown",
    grepl('other', PATIENT_RACE, fixed = TRUE) ~ "other",
    grepl('', PATIENT_RACE, fixed = TRUE) ~ "unknown",
    TRUE ~ PATIENT_RACE))

df$ARRIVAL_AGE_DI <- ifelse(
  df$ARRIVAL_AGE_DI == '85+', '85', df$ARRIVAL_AGE_DI
)
df$ARRIVAL_AGE_DI <- as.numeric(df$ARRIVAL_AGE_DI)

#=========================================================================
# Clean Chief Complaint --------------------------------------------------
#=========================================================================

complaints <- list(
  "Abdominal Complaints" = 
    c('Abdominal Cramping', 'Abdominal Distention', 'Dyspepsia',
      'Abdominal Pain', 'Ascites', 'Hernia', 
      'Abdominal Aortic Aneurysm', 'Abdominal Injury', "Pancreatitis",
      'Umbilical Hernia'),
  'Abnormal Test Results' = 
    c('Abnormal Lab', 'Abnormal Potassium', 'Abnormal Calcium', 
      'ECG Changes', 'Abnormal ECG', 'Abnormal Test Result', 
      'Blood Infection', 'Acute Renal Failure', 'Hypocalcemia',
      'Chronic Renal Failure', 'Pulmonary Embolism', 'Abnormal X-ray', 
      'Hypoglycemic Unawareness', 'Elevated Blood Pressure', 
      'Abnormal Sodium', 'Hyperglycemia', 'Hyponatremia', 
      'Platelet Disorders', 'Anemia', 'Hypoglycemia', 'Hypertension', 
      'Hypotension', 'Abnormal Chest Imaging', 'Abnormal Oximetry', 
      'Abnormal Stress Test', 'Blood Sugar Problem',
      'Hypocalcemia', 'Hyponatremia'),
  'Allergic Reaction' = 
    c('Allergic Reaction', 'Anaphylaxis'),
  'Back or Flank Pain' = 
    c('Back Pain', 'Back Problem', 'Flank Pain', 
      'Sciatica', 'Back Injury', 'Disc Disorder'),
  'Breast Complaints' = 
    c('Breast Mass', 'Breast Pain', 'Breast Problem', 'Breast Discharge',
      'Breast Cancer', 'Breast Discharge', 'Breast Inflammation'),
  'Cardiac Arrhythmias' = 
    c('Atrial Fibrillation', 'Atrial Flutter', 'Cardiac Valve Problem',
      'Bradycardia', 'Irregular Heart Beat', 'Palpitations', 'POTS', 
      'Ventricular Tachycardia','Rapid Heart Rate', 'Heart Problem', 
      'Cardiac Arrest', 'Congestive Heart Failure', 'Circulatory Problem',  
      "Transient Ischemic Attack", 'Ventricular Tachycardia'),
  'Chest Pain' = 
    c('Chest Injury', 'Chest Pain', 'Chest Wall Pain', 'Angina',
      'Collarbone Injury', 'Rib Injury', 'Heart Pain'),
  'Dizziness/Lightheadedness/Syncope' = 
    c('Dizziness', 'Near Syncope', 'Syncope', 'Vertigo', 'Spells',
      'Hypotension', 'Paroxysmal Positional Vertigo', 
      'Paroxysmal Positional Vertig'),
  'Ear Complaints' = 
    c('Cerumen Impaction', 'Ear Drainage', 'Ear Fullness',
      'Ear Laceration', 'Ear Problem', 'Earache',
      'Hearing Problem', 'Tinnitus', 'Ear Injury', 'Hearing Loss',
      'Nasal Trauma'),
  'Epistaxis' = 
    c('Epistaxis', 'Epistaxis (Nose Bleed)', 'Nose Problem'),
  'Exposures, Bites, and Envenomations' =
    c('Animal Bite', 'Body Fluid Exposure', 'Chemical Exposure', 
      'Poisoning', 'Exposure to STD', 'Insect Bite', 'Smoke Inhalation', 
      'Radiation', 'Snake Bite', 'Toxic Inhalation'),
  'Extremity Complaints' = 
    c('Ankle Injury', 'Ankle Pain', 'Arm Injury', 'Arm Pain', 
      'Cold Extremity', 'Arm Swelling', 'Arthritis', 'Elbow Injury', 
      'Elbow Pain', 'Pseudogout', 'Extremity Pain', 'Extremity Weakness', 
      'Finger Injury', 'Hip Injury', 'Extremity Weakness', 
      'Finger Injury', 'Finger Pain', 'Dislocation',
      'Foot Infection', 'Foot Injury', 'Foot Numbness', 'Foot Pain',
      'Foot Swelling', 'Foot Ulcer', 'Foot Wound Check', 'Hand Injury',
      'Hand Pain', 'Hip Pain', 'Joint Pain', 'Joint Swelling', 
      'Knee Injury', 'Knee Pain', 'Knee Problem', 'Leg Injury', 
      'Leg Pain', 'Leg Swelling', 'Lower Extremity Issue', 
      'Pain in Limb', 'Shoulder Injury', 'Shoulder Pain', 'Toe Injury', 
      'Upper Extremity Issue', 'Wrist Injury',
      'Wrist Pain', 'Hand Problem', 'Leg Cramps', 'Arm Problem',
      'Foot Problem', 'Pain In Limb', 'Toe Pain', 'Hand Burn',
      'Foot Burn', 'Neck Injury', 'Leg Problem', 'Deep Vein Thrombosis',
      "Varicose Veins",  "Ingrown Toenail", 'Hip Injury', 'Wound Care', 
      'Venous Thromboembolic Diseas'),
  'Eye Complaints' = 
    c('Blurred Vision', 'Decreased Visual Acuity', 'Diplopia', 
      'Detached Retina', 'Eye Drainage', 'Eye Exposure', 'Eye Pain', 
      'Eye Problem', 'Eye Swelling', 'Eye Trauma', 'Foreign Body Eye',
      'Flashes/Light', 'Loss of Vision', 'Red Eye', 'Visual Field Change', 
      'Eyelid Problem', 'Itchy Eye', 'Eye Exam', 'Burning Eyes', 
      'Eye Twitching', "Eyelid/brow Lift Evaluation",
      'Strabismus', 'Glaucoma',  "Spots/Floaters"),
  'Falls, Motor Vehicle Crashes, Assaults, and Trauma' =
    c('Assault Victim', 'Concussion', 'Facial Injury', 'Fall', 
      'Nasal Trauma', 'Head Injury', 'Head Laceration', 
      'Motor Vehicle Crash', 'Puncture Wound', 
      'Sexual Assault', 'Trauma', 'Domestic Violence', 'Gun Shot Wound',
      'Work Related Injury', 'Motorcycle Crash', 'Injury', 
      'Bicycle Accident', 'Near Drowning', 'Lip Laceration'),
  'Fatigue and Weakness' =
    c('Difficulty Walking', 'Fatigue', 'Gait Problem', 
      'Weakness-Generalized',
      'Chronic Fatique', 'Weakness - Generalized'),
  'Fevers, Sweats or Chills' =
    c('Chills', 'Diaphoresis', 'Fever', 'Night Sweats', 'Diaphoretic', 
      'Diapohresis', 'Hoarseness', 'Laryngitis'),
  'Foreign Body' =
    c('Food Bolus', 'Foreign Body', 'Foreign Body in Ear',
      'Foreign Body in Skin', 'Foreign Body in Vagina',
      'Swallowed Foreign Body', 'Foreign Body in Nose', 
      'Foreign Body', 'FB eye', 'Foreign Body in Rectum'),
  'Gastrointestinal Issues' =
    c('Anal Fissure', 'Black or Bloody Stool', 'Constipation', 'GERD', 
      'Anal Fistula', 'Diarrhea', 'Dysphagia', 'Fecal Impaction',
      'Fistula Follow Up', 'GIbleeding', 'GI Problem', 'Hemorrhoids', 
      'Morning Sickness', 'Nausea', 'Ostomy Care', 'Rectal Bleeding',
      'Rectal Pain', 'Vomiting', 'Vomiting Blood', 
      'Vomiting During Pregnancy',
      'GI Bleeding', 'Fecal Incontinence',
      'Bloated', 'Hematochezia', 'Urine Leakage', 'Heartburn',
      'Rectal Discharge', 'Urolithiasis', "Ulcerative Colitis",
      "Irritable Bowel Syndrome", "Rectal Prolapse","Fistula Evaluation",
      "Rectal Problems", 'Perianal Abscess', 'Fisula Evaluation', 
      'Stoma Dysfunction'),
  'Genital Complaints' =
    c('Groin Burn', 'Groin Pain', 'Groin Swelling', 'Inguinal Hernia', 
      'Menstrual Problem', 'Pelvic Pain', 'Penis Pain', 'Priapism', 
      'Testicle Pain', 'Menorrhagia',
      'Vaginal Bleed', 'Vaginal Bleeding', 'Vaginal Itching', 
      "Bartholin's Cyst", 'Genital Warts', 'Groin Injury', 
      'Vaginal Bleeding - Pregnant', 
      'Vag Bleed Pregnant', 'Female Genital Issue', 'Penis Injury',
      'Vaginal Discharge', 'Vaginal Pain', 'Erectile Dysfunction',
      'Vaginal Prolapse', 'Urethral Stricture', 'Penile Discharge',
      'Menorrhagia', "Gynecologic Exam", "Menstrual Problem", 
      "Vaginitis/Bacterial Vaginosis",
      'Ovarian Cyst', 'Vaginitis/Bacterial Vaginosi'),
  'Medical Device or Treatment Issue' =
    c('Cast Problem', 'Device Check', 'Dressing Change', 'Feeding Tube', 
      'AICD Problem', 'Insulin Pump Visit', 'Gastrostomy Tube Change', 
      'Medication Reaction', 'Shunt', 'Appliance Removal', 'Tube Problem', 
      'Urinary Catheter Change', 'Vascular Access Problem', 
      'Enteral Nutrition Evaluation', 'Device Malfunction', 
      'Pacemaker Problem', 
      'Removal / Exchange Catheter', 'Drain Removal', 
      'Outpatient Infusion', 
      'Treatment', 'Heart Assist Device', 'Stoma Dysfunction', 
      'Tracheostomy Tube Change',
      'Ureteral Stent Exchange'),
  'Medication Request' =
    c('Immunizations', 'Infusion/Injection Administration', 
      'IV Medication',  'Infusion/ Injection Administ', 'Med Refill',
      'Medication Visit', 
      'Pain Management', 'Blood Product Administration', 'Labs Only', 
      'Tetanus (Td & Tdap)', 'Wound Care'),
  'Neurological Issue' =
    c('Altered Mental Status', 'Cognitive Concerns', 'Facial Droop',
      'Pre Syncope', 'Focal Weakness', 'Headache', 'Memory Loss', 
      'Migraine', 'Dementia', 'Dysphasia',
      'Neuro Problem', 'Numbness', 'Paralysis', 'Seizures',
      'Slurred Speech', 'Spasms', 'Stroke Like Symptoms', 'Tingling',
      'Tremors', 'Trigeminal Neuralgia', 'Unable to Speak', 
      'Seizure Disorder', 'Insomnia', "Parkinson's Disease",
      'Loss of Consciousness', 'Neuropathy', 'Ataxia', 'Unable to speak',
      'Peripheral Neuropathy',
      'Stroke', 'Cerebrovascular Accident', 'Speech Problem', 
      'Acute Neurological Problem',
      'Flashes, Light', 'Unresponsive', "Multiple Sclerosis", 
      "Parkinson's Disease",
      "Febrile Seizure", 'Paresthesia', 'Peripheral Neuropathy', 
      'Hydrocephalus',
      'Spasticity', 'Neuroendocrine Tumor'),
  'Other' =
    c('Dehydration', 'Fisula Evaluation', 'Follow-Up', 'Illness', 
      'Letter for School/Work', 'Aneurysm',
      'Lung Eval', 'Error', 'Mass', 'Oral Swelling', 'Other', 
      'Advice Only', 'Deformity', 'Electric Shock',
      'Personal Problem', 'Shaking', 'Swelling', 'Swollen Glands', 
      'Adenopathy', 'Adrenal Problem',
      'Thrombophilia', 'Weight Gain', 'Weight Loss', 'Hiccups', '', 
      'Chemo Related Symptoms', 'Hot Flashes',
      'Follow-up', 'Non Healing Wound', '(Other)', 'Mouth Injury', 
      'Xerostomia', 'Prostate Check',
      'Suture / Staple Removal', 'Wellness', 'Voice Changes', 
      'Vital Sign Check', 'Coagulation Disorder',
      'Cold Exposure', 'Consult', 'Dental Problem', 
      'Tetanus (Td & Tdap)', "Infusion/ Injection Administ",
      "Tracheostomy Tube Change", 'Medical Information', 
      'Neutropenic Fever', 'Infection',
      'Leukemia',"Heat Exposure", "Poor Appetite", 'Gingivitis',
      "Pre-op Exam",
      'gingivitis', "Loss of appetite", "Failure To Thrive", 
      'Referral', 'Lymphoma',
      "Hot Flashes", 'Neutropenia', 'Radiation', 'Ingestion', 
      "TB Test", 'Fussy',
      'Lupus', 'Toxic Inhalation', 'Lung Screening', 
      'Leakage/Loss of Fluid',
      'Liver Eval', 'Hepatic Cancer', 'Lung Mass', 
      'Venous Thromboembolic Disease',
      'Insulin Pump Visit', 'Preventive Visit', 
      'Avulsion', 'Peripheral Edema',
      'Hypoglycemic Unawareness', 'Immobility', 
      'Giant Cell Arteritis', 'Polydipsia',
      'Platelet Disorders', 'Post-procedure', 'Lung Follow-up', 
      'Poisoning',
      'Injections', 'POTS', "Insulin Reaction", 
      'Liver Transplant', 'Labs Only'),
  'Other Pain' =
    c('Dental Pain', 'Facial Pain', 'Generalized Body Aches', 
      'Myalgia', 'Dental Injury',
      'Jaw Pain', 'Muscle Pain', 'Neck Pain', 'Pain', 
      'Sickle Cell Pain Crisis', 'Paresthesia',
      'Torticollis', 'Chronic Pain', 'Cancer Pain', 
      'Incisional Pain', 'Bone Pain',
      'Tailbone Pain', 'Gout', "Muscle pain/Weakness", 'Pseudogout'),
  'Post-Op Issue' =
    c('Post-Op', 'Post-Procedure', 'Post-Op Problem', 'Post-op',
      'Post-Op Issue', 'Wound Dehiscence', 'Post-op Problems',
      "Post-op Problem"),
  'Psychiatric Complaints' =
    c('Anxiety', 'Auditory Hallucinations', 'Depression', 'Panic Attack',
      'Homicidal', 'PTSD (Post-Traumatic Stress', 'Delusional', 'Fussy',
      'Paranoia', 'Suicide Attempt', 'Hallucinations', 'Manic Behavior',
      'Eating Disorder', 'Suicidal', 'Agitation', 'Psychiatric Evaluation',
      'Aggressive Behavior', 'Mental Health Problem', 
      'Inappropriate Words'),
  'Shortness of Breath' =
    c('Airway Obstruction', 'Aspiration', 'Pain With Breathing', 
      'Near Drowning',
      'Respiratory Distress', 'Shortness of Breath', 'Wheezing',
      'Increased Work Of Breathing',
      'Difficulty Breathing', 'Choking', "Oxygen Dependence",
      'Hyperventilating', 'Orthopnea'),
  'Skin Complaints' =
    c('Abrasion','Abscess', 'Bleeding/Bruising', 'Blister',
      'Angioedema', 'Lip Laceration',
      'Burn', 'Cellulitis', 'Cyst', 'Drainage from Incision', 
      'Disturb of Skin Sens',
      'Edema', 'Extremity Laceration', 'Facial Burn', 'Cyanosis',
      'Impetigo', 'Facial Laceration', 'Facial Swelling',
      'Finger Laceration', 'Leg Rash',
      'Herpes Zoster', 'Hives', 'Itching', 'Jaundice', 
      'Diabetic Ulcer', 'Diabetic Wound',
      'Laceration', 'Mouth Lesions', 'Non-Healing Wound', 'Rash', 
      'Recurrent Skin Infections', 'Skin Problem', 'Sore', 'Scabies',
      'Suture/Staple Removal', 'Wound Check', 'Wound Infection',
      'Lesion', 'Skin Check', 'Minor Skin Infection', 'Skin Ulcer',
      'Skin Discoloration', 'Sunburn', "Head Lice", 'Scabies',
      "Fungal Infection", "Leg Rash", 'Impetigo'),
  'Substance Abuse Issues' =
    c('Alcohol Intoxication', 'Alcohol Problem', 'Withdrawal',
      'Drug Overdose', 'Drug / Alcohol Dependency', 'Addiction Problem',
      'Addiction Assessment', 'Delirium Tremens (DTS)'),
  'Upper Respiratory Symptoms' =
    c('Congestion', 'Cough', 'Coughing Up Blood', 
      'Flu Symptoms', 'Enlarged Tonsils', 'Peritonsillar Abscess',
      'Nasal Congestion', 'Sinus Symptoms', 'Sinusitis',
      'Sore Throat', 'Hoarseness',
      'Throat Problem', 'Upper Respiratory Infection', 
      'Influenza', 'Laryngitis',
      'Respiratory Arrest', 'Pneumonia', 'Pleural Effusion', 'Asthma',
      'Croup', 'URI', 'Peritonsillar Abscess'),
  'Pregnancy Related' = 
    c("Pregnancy Problem", "Miscarriage", "Contractions", 
      'Ectopic Pregnancy',
      'Laboring',"Possible Pregnancy", "Pregnancy Related"), 
  'Renal' = 
    c('Av Fistula', 'Kidney Transplant', 'Elevated Serum Creatinine', 
      'End-Stage Liver Disease', "Hemodialysis Access", 'Nephritis',
      'Ureteral Stent Exchange'),
  'Urinary Complaints' =
    c('Bladder Problem', 'Blood in Urine', 'Cystitis',
      'Difficulty Urinating',
      'Dysuria', 'Gross Hematuria', 'Painful Urination',
      'Urinary Frequency', 'Urinary Symptom',
      'Urinary Incontinence', 'Urinary Problem', 
      'Urinary Retention', 'Slowing Urinary Stream',
      'Urinary Tract Infection', 'Urinary Urgency', 
      'Voiding Dysfunction',"Hesitancy Urinary")
)

df$complaint <- df$CHIEF_COMPLAINT

for (i in seq(1,length(complaints))){
  name <- names(complaints[i])
  complaint <- complaints[[i]]
  
  df$CHIEF_COMPLAINT <- ifelse(
    df$CHIEF_COMPLAINT %in% complaint, name, df$CHIEF_COMPLAINT
  )
}

df <- df %>%
  mutate(complaint_esi  = paste(ESI, CHIEF_COMPLAINT),
         complaint_esi = factor(complaint_esi))

#=========================================================================
# Categorize Vital Signs -------------------------------------------------
#=========================================================================

df$tachycardic <- ifelse(
  is.na(df$TRIAGE_PULSE) == F & df$TRIAGE_PULSE > 100, 1, 0
)
df$tachypneic <- ifelse(
  is.na(df$TRIAGE_RR)  == F  & df$TRIAGE_RR > 20, 1, 0
)
df$febrile <- ifelse(
  is.na(df$TRIAGE_TEMP)  == F  & df$TRIAGE_TEMP > 38, 1, 0
)
df$hypotensive <- ifelse(
  is.na(df$TRIAGE_SBP)  == F  & df$TRIAGE_SBP < 90, 1, 0
)

#=========================================================================
# Create Time FE ---------------------------------------------------------
#=========================================================================

df <- df %>%
  mutate(rel_minutes_arrival = ARRIVAL_DTTM_REL,
         rel_minutes_depart = rel_minutes_arrival + ED_LOS)

# Define arbitrary start date
start_date <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")

# Convert relative minutes to datetime
df$datetime <- start_date + minutes(df$rel_minutes_arrival)

# Extract hour of day, day of week, and month of year
df$hour_of_day <- hour(df$datetime)
df$hour_of_day <- cut(df$hour_of_day, breaks = c(-Inf, 6, 12, 18, Inf), 
                      labels = c("Block 1", "Block 2", "Block 3", "Block 4"))

df$day_of_week <- weekdays(df$datetime)
df$month_of_year <- month(df$datetime, label = TRUE)
df$dayofweekt <- paste(df$day_of_week, df$hour_of_day)

# Calculate the number of patients in the hospital at the time of 
# each patient's arrival
df$patients_in_hospital <- sapply(df$rel_minutes_arrival, 
                                  function(arrival_time) {
                                    sum(df$rel_minutes_arrival <= arrival_time & 
                                          df$rel_minutes_depart > arrival_time
                                    ) - 1
                                 })

#=========================================================================
# Create Final Dataset ---------------------------------------------------
#=========================================================================

df$ln_ED_LOS <- log(df$ED_LOS)

final <- df %>% 
  mutate(RTN_72_HR = ifelse(RTN_72_HR == 'Y', 1, 0),
         RTN_72_HR_ADMIT = ifelse(RTN_72_HR_ADMIT == 'Y', 1, 0)) 

final$admit = ifelse(final$ED_DISPOSITION == 'Admit', 1, 0)
final$discharge = ifelse(final$ED_DISPOSITION == 'Discharge', 1, 0)
final$observation = ifelse(final$ED_DISPOSITION == 'Observation', 1, 0)

final <- final %>%
  mutate(dispo = case_when(
    admit == 1 ~ 'admit',
    discharge == 1 ~ 'discharge',
    observation == 1 ~ 'observeration',
    TRUE ~ 'other')) 

# Limit dataset to only physicians that had more than 520 encounters
provider_counts <- table(final$ED_PROVIDER)
providers_less_than_500 <- names(provider_counts[provider_counts < 520])
final <- final[!(final$ED_PROVIDER %in% providers_less_than_500), ]
final$complaint_esi <- paste(final$CHIEF_COMPLAINT, final$ESI)
final <- filter(final, !is.na(final$ESI))

rm(list = setdiff(ls(), "final"))


final <- final %>%
  filter(ED_LOS < 6000, ED_LOS > 0) 



# add experience variable for physicians 
final <- final %>%
  mutate(EXPERIENCE = case_when(
    ED_PROVIDER == 'JUDSON, KURTIS A' ~ 2006,
    ED_PROVIDER == 'KOMARA, JAMES S' ~ 1985,
    ED_PROVIDER == 'RAPPAPORT, DOUGLAS E' ~ 2016,
    ED_PROVIDER == 'MONAS, JESSICA' ~ 2011,
    ED_PROVIDER == 'KELLEY, JAMES' ~ 1994,
    ED_PROVIDER == 'DRECHSEL, KEVIN M' ~ 1998,
    ED_PROVIDER == 'DIETRICH, BOB D' ~ 1999,
    ED_PROVIDER == 'URUMOV, ANDREJ' ~ 2005,
    ED_PROVIDER == 'MACY, CHERYL' ~ 2011,
    ED_PROVIDER == 'BRAND, SHARI I.' ~ 1999,
    ED_PROVIDER == 'TRAUB, STEPHEN J' ~ 1998,
    ED_PROVIDER == 'CHANTLER, EDMUNDO L' ~ 2005,
    ED_PROVIDER == 'HAY-ROE, NEIL' ~ 1987,
    ED_PROVIDER == 'GAUHAROU, ERIK SHAWN' ~ 1999,
    ED_PROVIDER == 'HODGSON, NICOLE R' ~ 2017,
    ED_PROVIDER == 'STEWART, CHRISTOPHER F' ~ 1999,
    ED_PROVIDER == 'TORRES, MARCELLA' ~ 1995,
    ED_PROVIDER == 'MAHER, STEVEN A' ~ 2004,
    ED_PROVIDER == 'ARNOLD, RICKY R' ~ 1997,
    ED_PROVIDER == 'VINCIJANOVIC, LISA M' ~ 2007,
    ED_PROVIDER == 'MORRO, DAVID C' ~ 2003,
    ED_PROVIDER == 'PETRI, ROLAND W' ~ 1988,
    ED_PROVIDER == 'KOZAK, PAUL A' ~ 1993,
    ED_PROVIDER == 'LINDOR, RACHEL A' ~ 2017,
  ))

final$EXPERIENCE <- 2019 - final$EXPERIENCE

final <- final %>%
  mutate(PROVIDER_SEX = case_when(
    ED_PROVIDER == 'JUDSON, KURTIS A' ~ 'M',
    ED_PROVIDER == 'KOMARA, JAMES S' ~ 'M',
    ED_PROVIDER == 'RAPPAPORT, DOUGLAS E' ~ 'M',
    ED_PROVIDER == 'MONAS, JESSICA' ~ 'F',
    ED_PROVIDER == 'KELLEY, JAMES' ~ 'M',
    ED_PROVIDER == 'DRECHSEL, KEVIN M' ~ 'M',
    ED_PROVIDER == 'DIETRICH, BOB D' ~ 'M',
    ED_PROVIDER == 'URUMOV, ANDREJ' ~ 'M',
    ED_PROVIDER == 'MACY, CHERYL' ~ 'F',
    ED_PROVIDER == 'BRAND, SHARI I.' ~ 'F',
    ED_PROVIDER == 'TRAUB, STEPHEN J' ~ 'M',
    ED_PROVIDER == 'CHANTLER, EDMUNDO L' ~ 'M',
    ED_PROVIDER == 'HAY-ROE, NEIL' ~ 'M',
    ED_PROVIDER == 'GAUHAROU, ERIK SHAWN' ~ 'M',
    ED_PROVIDER == 'HODGSON, NICOLE R' ~ 'F',
    ED_PROVIDER == 'STEWART, CHRISTOPHER F' ~ 'M',
    ED_PROVIDER == 'TORRES, MARCELLA' ~ 'F',
    ED_PROVIDER == 'MAHER, STEVEN A' ~ 'M',
    ED_PROVIDER == 'ARNOLD, RICKY R' ~ 'M',
    ED_PROVIDER == 'VINCIJANOVIC, LISA M' ~ 'F',
    ED_PROVIDER == 'MORRO, DAVID C' ~ 'M',
    ED_PROVIDER == 'PETRI, ROLAND W' ~ 'M',
    ED_PROVIDER == 'KOZAK, PAUL A' ~ 'M',
    ED_PROVIDER == 'LINDOR, RACHEL A' ~ 'F',
  ))


final <- final %>% 
  mutate(date = as.Date(datetime))

final <- final %>% 
  group_by(ED_PROVIDER, date) %>% 
  mutate(patients_seen = n(),
         patient_order_of_day = row_number(),
         patients_tbs = patients_seen - patient_order_of_day) %>%
  ungroup()


final <- final %>%
  mutate(
    Time_arrival_to_triage = TRIAGE_COMPLETED_REL - ARRIVAL_DTTM_REL,
    Time_triage_to_firstcontact = FIRST_CONTACT_DTTM_REL - TRIAGE_COMPLETED_REL,
    Time_to_Result_Lab = LAB_RESULT_DTTM_REL - LAB_ORDER_REL,
    Time_to_Result_XRay = PLAIN_XRAY_RESULT_DTTM_REL - XR_ORDER_REL,
    Time_to_Result_Ultrasound = US_RESULT_DTTM_REL - US_ORDER_DTTM_REL,
    Time_to_Result_CTcon = CT_WITH_CONTR_RESULT_DTTM_REL - CON_CT_ORDER_REL,
    Time_to_Result_CTnon = CT_WITHOUT_CONTR_RESULT_DTTM_REL - NON_CON_CT_ORDER_REL
  )

final$ln_ttr <- log(final$total_testing_time)

#=========================================================================
##########################################################################
#=========================================================================
# IV Construction --------------------------------------------------------
#=========================================================================
##########################################################################
final$imaging <- ifelse(final$imgTests > 0, 1, 0)

final$residual_batch <- resid(
  felm(batched ~ ARRIVAL_AGE_DI + patients_in_hospital | dayofweekt + month_of_year + complaint_esi + LAB_PERF |0|ED_PROVIDER, data=final)
  )

final$residual_ntests <- resid(
  felm(imgTests ~ ARRIVAL_AGE_DI + patients_in_hospital | dayofweekt + month_of_year + complaint_esi + LAB_PERF |0|ED_PROVIDER, data=final))

# Step 2: get batch tendency for each provider
final <- final %>%
  group_by(ED_PROVIDER) %>%
  mutate(Sum_Resid=sum(residual_batch, na.rm=T),
         batch.tendency = (Sum_Resid - residual_batch) / (n() - 1),
         
         Sum_Resid_ntests=sum(residual_ntests, na.rm=T),
         test.inclination = (Sum_Resid_ntests - residual_ntests) / (n() - 1)) %>%
  ungroup()


write.csv(final, 'outputs/data/all_clean.csv')


# keep complaints that appear more than 500 times
complaint_counts <- table(final$CHIEF_COMPLAINT)
complaints_less_than_500 <- names(complaint_counts[complaint_counts < 1000])
final <- final[!(final$CHIEF_COMPLAINT %in% complaints_less_than_500), ]



write.csv(final, 'outputs/data/final.csv')

