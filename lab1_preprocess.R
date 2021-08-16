library(tidyverse)
source("filling.R")

# read data ---------------------------------------------------------------
adms = read_csv("Data/ADMISSIONS.csv.gz")
dgns = read_csv("Data/DIAGNOSES_ICD.csv.gz")
prsc = read_csv("Data/PRESCRIPTIONS.csv.gz")
labt = read_csv("Data/LABEVENTS.csv.gz")

# preprocess --------------------------------------------------------------
dgns_new = dgns %>% 
  mutate(event = paste0("ICD_", ICD9_CODE)) %>% 
  select(SUBJECT_ID, HADM_ID, event) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(event) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n>=50) %>% 
  select(-n)
rm(dgns)
gc()

prsc_new = prsc %>% 
  mutate(event = paste0("drug_", FORMULARY_DRUG_CD)) %>% 
  select(SUBJECT_ID, HADM_ID, event) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(event) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n>=50) %>% 
  select(-n)
rm(prsc)
gc()

labt_new = labt %>% 
  mutate(event = paste0("lab_", ITEMID)) %>% 
  select(SUBJECT_ID, HADM_ID, event) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(event) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n>=50) %>% 
  select(-n)
rm(labt)
gc()

work_data = rbind(dgns_new, prsc_new, labt_new)

saveRDS(work_data, "Data/MIMIC_pre.rds")
# work_data = readRDS("Data/MIMIC_pre.rds")

# subj with 2 visits and more ---------------------------------------------
adms = read_csv("Data/ADMISSIONS.csv.gz")

adms_new = adms %>% 
  select(SUBJECT_ID, HADM_ID, ADMITTIME) %>% 
  filter(HADM_ID %in% work_data$HADM_ID) %>% 
  na.omit() %>% 
  distinct()

adms_new = adms_new %>% 
  group_by(SUBJECT_ID) %>% 
  arrange(ADMITTIME) %>% 
  mutate(n_adms = n(),
         time_rank = row_number()) %>% 
  ungroup() %>% 
  filter(n_adms>=2)

work_data = work_data %>% 
  left_join(adms_new) %>% 
  na.omit()

# most recent visit -------------------------------------------------------
work_data_outcome = work_data %>% 
  group_by(SUBJECT_ID) %>% 
  filter(time_rank==max(time_rank)) %>% 
  ungroup()

outcome = work_data_outcome %>% 
  filter(str_detect(event, "ICD")) %>% 
  select(SUBJECT_ID, event)

# previous visits ---------------------------------------------------------
work_data_predictor = work_data %>% 
  group_by(SUBJECT_ID) %>% 
  filter(time_rank!=max(time_rank)) %>% 
  ungroup()

predictor = work_data_predictor %>% 
  select(SUBJECT_ID, event)

length(unique(work_data_predictor$SUBJECT_ID)) #7519
length(unique(work_data_predictor$HADM_ID)) #12431
length(unique(work_data_predictor$event)) #3129

# matrix ------------------------------------------------------------------
all_patient = unique(work_data_outcome$SUBJECT_ID)

outcome = filling_01_bycol(person_id = outcome$SUBJECT_ID,
                           concept_id = outcome$event,
                           all_person = all_patient)

predictor = filling_count(person_id = predictor$SUBJECT_ID,
                          concept_id = predictor$event,
                          all_person = all_patient)

identical(rownames(outcome), rownames(predictor))

saveRDS(outcome, "Data/outcome.rds")
saveRDS(predictor, "Data/predictor.rds")




