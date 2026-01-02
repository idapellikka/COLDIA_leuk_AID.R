##PROJECT COLDIA: Analyses for article: Association of childhood acute leukemia and autoimmune diseases

#Load necessary packages
library(lubridate)
library(dplyr)
library(survival)
library(tidyverse)
library(tidyr)
library(ggplot2)

#--- Autoimmune datasets ---

#Import datasets containing information about autoimmune disease diagnoses
#from Care Register for Health Care

hilmo1 <- read.csv("FD_2022_1096_THL_THL2022_1096_avohilmo_auto.csv", sep = ";")
hilmo2 <- read.csv("FD_2022_1096_THL_THL2022_1096_Hilmo_auto.csv", sep = ";")
hilmo3 <- read.csv("FD_2022_1096_THL_THL2022_1096_hilmo72_86_auto.csv", sep = ";")
hilmo4 <- read.csv("FD_2022_1096_THL_THL2022_1096_hilmo87_93_auto.csv", sep = ";")
hilmo5 <- read.csv("FD_2022_1096_THL_THL2022_1096_hilmo94_95_auto.csv", sep = ";")

#Select releavant columns and rename columns in each dataset for consistency
hilmo1 <- hilmo1 %>%
  select(FID, KAYNTI_ALKOI, ICD10) %>% #select relevant columns
  rename(ID0 = FID, dg_date = KAYNTI_ALKOI) %>% #rename
  mutate(dg_date = as.Date(dg_date, format = "%d.%m.%Y"))

hilmo2 <- hilmo2 %>%
  select(ID0 = FID, dg_date = TUPVA, ICD10 = KOODI) %>%
  mutate(dg_date = as.Date(dg_date, format = "%d.%m.%Y"))

#For third dataset, modify from wide form to long form
hilmo3 <- hilmo3 %>%
  mutate(across(c(starts_with("DG")), #Make sure that the columns containing ICD9
                as.character)) %>%    #codes are in same format
  pivot_longer(cols = starts_with("DG"), #Mutate from wide to long form
               values_to = "ICD9",
               names_to = NULL) %>%
  select(ID0 = FID, dg_date = TULOPV,ICD9) %>%
  mutate(dg_date = as.Date(dg_date, format = "%d.%m.%Y"))

#Do same for the 4th dataset
hilmo4 <- hilmo4 %>%
  mutate(across(c("PDG", starts_with("SDG")),
                as.character)) %>%
  pivot_longer(cols = c("PDG", starts_with("SDG")), #Combine columns containing ICD9
               values_to = "ICD9", #Name combined column ICD9
               names_to = NULL) %>% #Keep original values
  select(ID0 = FID, dg_date = TUPVA,ICD9) %>% #Select and rename relevant columns
  mutate(dg_date = as.Date(dg_date, format = "%d.%m.%Y"))

#5th dataset
#Columns PDGO, PDGE, SDG1O, SDG1E, SDG2O, SDG2E are ICD10 codes
#Columns PDG, SDG1 and SDG2 are ICD9 codes
hilmo6 <- hilmo5 %>%
  mutate(across(c(starts_with("PDG"), starts_with("SDG")),
                as.character)) %>%
  #Combine ICD10 code columns
  pivot_longer(cols = c(PDGO, PDGE, SDG1O, SDG1E, SDG2O, SDG2E),
               values_to = "ICD10", #Name combined column ICD10
               names_to = NULL) %>% #Keep original values
  select(ID0 = FID, dg_date = TUPVA,ICD10) %>% #Select and rename relevant columns
  mutate(dg_date = as.Date(dg_date, format = "%d.%m.%Y"))

#ICD9 codes
hilmo5 <- hilmo5 %>%
  pivot_longer(c(PDG, SDG1, SDG2),
               values_to = "ICD9",
               names_to = NULL) %>%
  select(ID0 = FID, dg_date = TUPVA, ICD9) %>% #Select and rename relevant columns
  mutate(dg_date = as.Date(dg_date, format = "%d.%m.%Y"))

#Merge all datasets from the Care Register for Healthcare
AID <- bind_rows(hilmo1, hilmo2, hilmo3, hilmo4, hilmo5, hilmo6)
#Remove dots from ICD-codes for consistency
AID <- AID %>%
  mutate(ICD10 = gsub("\\.", "", ICD10))

#Create a table with autoimmune diseases and their ICD10 and ICD9 patterns
disease_codes <- data.frame(
  #list AIDs
  disease_name = c("addisons_disease", "aiha", "hepatitis", "thyroiditis", "AS",
                   "basedows_disease", "bechets_disease", "crohns_disease",
                   "celiac_disease", "dermatomyocitis", "GPA", "ITP", "juvenile_arthritis",
                   "ms_disease", "pemphigoid", "pemphigus", "pernicious_anemia",
                   "polyarteris_nodosa", "PMR", "psoriasis", "RF", "RA", "sarcoidosis",
                   "sjogren", "SLE", "systemic_sclerosis", "colitis_ulcerosa"),
  
  #List of corresponding ICD10 codes
  ICD10_pattern = c("^E27[124]", "^D591", "^K754", "^E063", "^M45", "^E050", "M352", "^K50[0189]",
                    "^K900", "^M33[0129]", "^M313", "^D693", "^M08[01348]", "^G35", "^L12",
                    "^L10[0123489]", "^D51[01289]", "M30[01238]", "^M315|^M353", 
                    "^L40", "^I0[012]", "^M0[56]", "^D86[012389]", "^M350",
                    "^M32", "^M34[0189]", "^K51[013589]"),
  
  #List of corresponding ICD9 codes
  ICD9_pattern = c("2554", "2830", "^57142", "^2452", "^720", "^2420", "^1361", "^555", "^5790", "^710[34]",
                   "^4460", "^2873", "^7143", "^340", "^6945", "^6944", "^2810", "^4464", "^725",
                   "^696", "^39[012]", "714", "^135", "^7102", "^7100", "^7101", "^556")
)


#Create a function to filter AID data to rows that match either ICD10 or ICD9 pattern
process_disease <- function(data, disease_name, ICD10_pattern, ICD9_pattern) {
  subset_data <- data %>%
    filter(str_detect(ICD10, ICD10_pattern) |
             str_detect(ICD9, ICD9_pattern) #Keep only the rows where ICD10 or ICD9 matches
    ) %>%
    arrange(ID0, dg_date) %>% #Sort so that the earliest diagnosis date appears first
    group_by(ID0) %>%         #Group by individual ID
    slice(which.min(dg_date)) %>% #Keep only the first diagnosis date per person
    ungroup() %>%
    rename(!!paste0(disease_name, "_date") := dg_date) %>% #Rename date column specific to sorted AID
    mutate(!!disease_name := 1) %>% #Add binary column = 1 for AID present
    select(ID0, ends_with("_date"), all_of(disease_name)) %>%
    distinct(ID0, .keep_all = TRUE) #Ensure only one row per individual per AID
  
  return(subset_data)
}

#Apply the function to all autoimmune diseases, this creates a list of tiny data frames
AID_diseases <- lapply(
  seq_len(nrow(disease_codes)),
  function(i) {
    process_disease(AID,
                    disease_codes$disease_name[i],
                    disease_codes$ICD10_pattern[i],
                    disease_codes$ICD9_pattern[i]
    )
  }
)

#Combine into one wide dataset
AID_wide <- purrr::reduce(
  AID_diseases,
  .init = AID %>% distinct(ID0),
  .f = ~ full_join(.x, .y, by = "ID0")) %>%
  #Replace NA flags (no disease) with 0
  mutate(across(disease_codes$disease_name, ~ replace_na(.,0))
  )

# Localized scleroderma is in different dataset: repeat steps to identify individuals
# diagnosed with localized scleroderma

#Import datasets
LS1 <- read.csv("FD_5869_THL2023_5869_AVOHILMO.csv", sep = ";")
LS2 <- read.csv("FD_5869_THL2023_5869_HILMO.csv", sep = ";")
LS3 <- read.csv("FD_5869_THL2023_5869_HILMO72_86.csv", sep = ";")
LS4 <- read.csv("FD_5869_THL2023_5869_HILMO87_93.csv", sep = ";")
LS5 <- read.csv("FD_5869_THL2023_5869_HILMO94_95.csv", sep = ";")

#Select releavant columns and rename columns in each dataset for consistency
LS1 <- LS1 %>%
  select(FID, KAYNTI_ALKOI, ICD10) %>% #select relevant columns
  rename(ID0 = FID, dg_date = KAYNTI_ALKOI) %>% #rename
  mutate(dg_date = as.Date(dg_date, format = "%d.%m.%Y"))

LS2 <- LS2 %>%
  select(ID0 = FID, dg_date = TUPVA, ICD10 = KOODI) %>%
  mutate(dg_date = as.Date(dg_date, format = "%d.%m.%Y"))

#For third dataset, modify from wide form to long form
LS3 <- LS3 %>%
  mutate(across(c(starts_with("DG")), #Make sure that the columns containing ICD9
                as.character)) %>%    #codes are in same format
  pivot_longer(cols = starts_with("DG"), #Mutate from wide to long form
               values_to = "ICD9",
               names_to = NULL) %>%
  select(ID0 = FID, dg_date = TULOPV,ICD9) %>%
  mutate(dg_date = as.Date(dg_date, format = "%d.%m.%Y"))

#Do same for the 4th dataset
LS4 <- LS4 %>%
  mutate(across(c("PDG", starts_with("SDG")),
                as.character)) %>%
  pivot_longer(cols = c("PDG", starts_with("SDG")), #Combine columns containing ICD9
               values_to = "ICD9", #Name combined column ICD9
               names_to = NULL) %>% #Keep original values
  select(ID0 = FID, dg_date = TUPVA,ICD9) %>% #Select and rename relevant columns
  mutate(dg_date = as.Date(dg_date, format = "%d.%m.%Y"))

#5th dataset
#Columns PDGO, PDGE, SDG1O, SDG1E, SDG2O, SDG2E are ICD10 codes
#Columns PDG, SDG1 and SDG2 are ICD9 codes
LS6 <- LS5 %>%
  mutate(across(c(starts_with("PDG"), starts_with("SDG")),
                as.character)) %>%
  #Combine ICD10 code columns
  pivot_longer(cols = c(PDGO, PDGE, SDG1O, SDG1E, SDG2O, SDG2E),
               values_to = "ICD10", #Name combined column ICD10
               names_to = NULL) %>% #Keep original values
  select(ID0 = FID, dg_date = TUPVA,ICD10) %>% #Select and rename relevant columns
  mutate(dg_date = as.Date(dg_date, format = "%d.%m.%Y"))

#ICD9 codes
LS5 <- LS5 %>%
  pivot_longer(c(PDG, SDG1, SDG2),
               values_to = "ICD9",
               names_to = NULL) %>%
  select(ID0 = FID, dg_date = TUPVA, ICD9) %>% #Select and rename relevant columns
  mutate(dg_date = as.Date(dg_date, format = "%d.%m.%Y"))

#Merge all datasets from the Care Register for Healthcare
LS <- bind_rows(LS1, LS2, LS3, LS4, LS5, LS6)
#Remove dots from ICD-codes for consistency
LS <- LS %>%
  mutate(ICD10 = gsub("\\.", "", ICD10))

#Manually collect individuals diagnosed with localized scleroderma

LS <- LS %>%
  filter(str_detect(ICD10, "L940")
         | str_detect(ICD9, "7010")) %>% #Detect cases with a diagnosis of LS
  arrange(ID0, dg_date) %>% #Sort the data by ID0 and autoimmune diagnosis date
  group_by(ID0) %>% #Group data by individual ID
  slice(which.min(dg_date)) %>% #Select the earliest autoimmune diagnosis for each individual
  ungroup() #Remove grouping

LS <- LS %>%
  mutate(localized_scleroderma = ifelse( #Create new variable for LS
    grepl("L940", ICD10) | grepl("7010", ICD9),
    1, #mark as 1 if LS is present
    0 #Otherwise mark 0
  )) %>%
  rename(localized_scleroderma_date = dg_date) %>% #Rename date column
  select(ID0, localized_scleroderma_date, localized_scleroderma) %>% #Select relevant columns
  arrange(ID0, localized_scleroderma_date) %>% #Sort data by ID0 and date
  distinct(ID0, .keep_all = TRUE) #Ensure only one row per individual keeping the earliest dg

#Merge to AID_wide dataset
AID_wide <- bind_rows(AID_wide, LS)

# --- Type 1 diabetes ---

#Collect type 1 diabetes cases from the medication reimbursement data 
#from the Social insurance Institution of Finland
#Identify cases by using codes 102, 171, 177, 371, 382
#ICD10 codes E10 
#ICD9 codes 250, 2500, 2501, 2500B, 2502B, 2507B

dm_reimbursements <- read.csv("FD_2022_1096_Kela_94_522_2022_Laakekorvausoikeudet.csv", sep = ";")

#Rename columns for consistency and keep only relevant columns
dm_reimbursements <- dm_reimbursements %>%
  rename(ID0 = FID, T1DM_date = KORVAUSOIKEUS_ALPV) %>%
  select(ID0, T1DM_date, KORVAUSOIKEUS_KOODI, DIAGNOOSI_KOODI)

#Restrict to medical reimbursement codes associated with diabetes and
t1dm_data <- dm_reimbursements %>%
  filter(
    str_detect(KORVAUSOIKEUS_KOODI, "103|171|177|371|382")
  )

#Restrict to diagnosis codes associated with type 1 diabetes and leave empty dg codes
t1dm_data <- t1dm_data %>%
  filter(str_detect(DIAGNOOSI_KOODI, "^E10") |
           str_detect( DIAGNOOSI_KOODI, "250") |
           str_detect(DIAGNOOSI_KOODI, "2500") |
           str_detect(DIAGNOOSI_KOODI, "2501") |
           str_detect(DIAGNOOSI_KOODI, "2500B") |
           str_detect(DIAGNOOSI_KOODI, "2502B") |
           str_detect(DIAGNOOSI_KOODI, "2507B") |
           str_detect(DIAGNOOSI_KOODI, " ")
  )

#Create variable where t1dm = 1 if type 1 diabetes mellitus is present and otherwise 0
t1dm_data$T1DM <- 1

#Select only relevant columns and remove duplicates
t1dm_data <- t1dm_data %>%
  select(ID0, T1DM_date, T1DM) %>%
  distinct(ID0, .keep_all = TRUE)


#--- CONFOUNDING FACTORS ---

#Information about Down syndrome
#Import Down syndrome registries
down1 <- read.csv("FD_2022_1096_THL_thl2022_1096_epamuodostuma.csv", sep = ";")
down2 <- read.csv("FD_2022_1096_THL_THL2022_1096_hilmo87_92_down.csv", sep = ";")

#Harmonize ID column name across datasets and 
#mark Down syndrome cases with indicating variable
down1 <- down1 %>%
  rename(ID0 = FID) %>%
  select(ID0)
down1$down <- 1

down2 <- down2 %>%
  rename(ID0 = FID) %>%
  select(ID0)
down2$down <- 1

#Combine the two Down syndrome datasets and remove any duplicates
down <- rbind(down1, down2) %>%
  distinct(down, ID0, .keep_all = TRUE)

#Information about acute pancreatitis

#Import datasets containing information about acute pancreatitis diagnoses
panc1 <- read.csv("FD_5869_THL2023_5869_AVOHILMO.csv", sep = ";")
panc2 <- read.csv("FD_5869_THL2023_5869_hilmo_2.csv", sep = ";")
panc3 <- read.csv("FD_5869_THL2023_5869_HILMO72_86.csv", sep = ";")
panc4 <- read.csv("FD_5869_THL2023_5869_HILMO87_93.csv", sep = ";")
panc5 <- read.csv("FD_5869_THL2023_5869_HILMO94_95.csv", sep = ";")

#Identify acute pancreatitis cases based on ICD codes in each dataset

#Dataset 1
values_panc1 <- c("K85", "K85.3#", "K85.8", "K85.9") #ICD10 codes indicating acute panc
panc1 <- panc1[panc1$ICD10 %in% values_panc1, ] %>%
  rename(ID0 = FID) %>%
  select(ID0)

#Dataset 2
values_panc2  <- c("K85", "K850", "K853", "K858", "K859") #ICD10 codes indicating acute panc
panc2 <- panc2[panc2$KOODI %in% values_panc2, ] %>%
  rename(ID0 = FID) %>%
  select(ID0)

#Dataset 3
columns_panc3 <- c(panc3$DG1, panc3$DG2, panc3$DG3, panc3$DG4)
values_panc3 <- c("5770") #ICD9 code indicating acute panc
panc3 <- panc3[columns_panc3 %in% values_panc3, ] %>%
  rename(ID0 = FID) %>%
  select(ID0)

#Dataset 4
columns_panc4 <- c(panc4$PDG, panc4$SDG1, panc4$SDG2, panc4$SDG3)
values_panc4 <- c("5570E")
panc4 <- panc4[columns_panc4 %in% values_panc4, ] %>%
  rename(ID0 = FID) %>%
  select(ID0)

#Dataset 5
#Columns PDGO, PDGE, SDG1O, SDG1E, SDG2O, SDG2E are ICD10 codes
#Columns PDG, SDG1 and SDG2 are ICD9 codes
columns_panc5 <- c(panc5$PDG0, panc5$PDGE, panc5$SDG1O, panc5$SDG1E,
                   panc5$SDG2O, panc5$SDG2E)
values_panc5 <- c("K85", "K850", "K853", "K858", "K859") #ICD10 codes indicating acute panc
panc5_ICD10 <- panc5[columns_panc5 %in% values_panc5, ] %>%
  rename(ID0 = FID) %>%
  select(ID0)

columns_panc5_ICD9 <- c(panc5$PDG, panc5$SDG1, panc5$SDG2)
values_panc5_ICD9 <- c("5770A") #ICD9 codes for acute pancreatitis
panc5 <- panc5[columns_panc5_ICD9 %in% values_panc5_ICD9, ] %>%
  rename(ID0 = FID) %>%
  select(ID0)

#Combine pancreatitis datasets and remove duplicates 
pancreatitis <- rbind(panc1, panc2, panc5) %>% #panc3 and panc4 had no relevant columns
  distinct(ID0, .keep_all = TRUE) 

#Mark pancreatitis cases with indicator variable
pancreatitis$panc <- 1

#Information from the Birth Registry
synre <- read.csv("FD_2022_1096_THL_THL2022_1096_synre.csv", sep = ";")

#Rename the identifying column for consistency and select only relevant columns
synre <- synre %>%
  rename(ID0 = LAPSI_FID) %>%
  #Large for gestational age (SGAC) and smoking (TUPAKOINTITUNNUS) were identified 
  #as confounding factors in our previous study
  select(ID0, SGAC, TUPAKOINTITUNNUS)

#Recode SGAC variable: 1 and 2 = not LGA -> 0, 3 = LGA -> 1, 9 = missing information -> NA
synre$SGAC <- with(synre, ifelse(SGAC %in% c(1,2), 0,
                                 ifelse(SGAC == 3, 1,
                                        ifelse(SGAC == 9, NA, SGAC))))
names(synre)[names(synre) == "SGAC"] <- "LGA" #Rename as LGA for clarity

#Recode maternal smoking variable
# 2,3 and 4 = smoked -> 1, 1 = non smoking -> 1, 9 = missing information -> NA
synre$TUPAKOINTITUNNUS <- with(synre, ifelse(TUPAKOINTITUNNUS %in% c(2,3,4), 1,
                                             ifelse(TUPAKOINTITUNNUS == 1, 0,
                                                    ifelse(TUPAKOINTITUNNUS == 9, NA, TUPAKOINTITUNNUS))))

# --- LEUKEMIA DATASET ---

#Import dataset about leukemia cases and controls
leuk_cases_controls <- read.csv("leuk_cases_controls.csv", sep = ",")

#Merge information about autoimmune diseases and confounding factors
#Autoimmune diseases
leuk_AID <- merge(leuk_cases_controls, AID_wide, by="ID0", all.x = TRUE)
#Type 1 diabetes mellitus
leuk_AID <- merge(leuk_AID, t1dm_data, by="ID0", all.x = TRUE)
#Down syndrome and pancreatitis
leuk_AID <- merge(leuk_AID, down, by="ID0", all.x = TRUE)
leuk_AID <- merge(leuk_AID, pancreatitis, by="ID0", all.x = TRUE)
#Variables from the Birth Registry
leuk_AID <- merge(leuk_AID, synre, by="ID0", all.x = TRUE)

#Create a list of autoimmune disease variables and dates
disease_list = c("addisons_disease", "aiha", "hepatitis", "thyroiditis", "AS",
                 "basedows_disease", "bechets_disease", "crohns_disease",
                 "celiac_disease", "dermatomyocitis", "GPA", "ITP", "juvenile_arthritis",
                 "ms_disease", "pemphigoid", "pemphigus", "pernicious_anemia",
                 "polyarteris_nodosa", "PMR", "psoriasis", "RF", "RA", "sarcoidosis",
                 "sjogren", "SLE", "systemic_sclerosis", "colitis_ulcerosa",
                 "localized_scleroderma", "T1DM")

disease_list_wo_T1DM = c("addisons_disease", "aiha", "hepatitis", "thyroiditis", "AS",
                         "basedows_disease", "bechets_disease", "crohns_disease",
                         "celiac_disease", "dermatomyocitis", "GPA", "ITP", "juvenile_arthritis",
                         "ms_disease", "pemphigoid", "pemphigus", "pernicious_anemia",
                         "polyarteris_nodosa", "PMR", "psoriasis", "RF", "RA", "sarcoidosis",
                         "sjogren", "SLE", "systemic_sclerosis", "colitis_ulcerosa",
                         "localized_scleroderma")

AID_date_list = paste0(disease_list, "_date")

AID_date_list_wo_T1DM = paste0(disease_list_wo_T1DM, "_date")

#Across autoimmune disease variables,replace NA flags (no disease) with 0
leuk_AID <- leuk_AID %>%
  mutate(across(all_of(disease_list), ~ replace(., is.na(.), 0)))

#Create variables indicating if any autoimmune disease is present and diagnosis date for that
leuk_AID$any_ai_disease <- as.integer(rowSums(leuk_AID[disease_list] == 1, na.rm = TRUE) > 0)

#Create any AID date variable
#Ensure all variables are in date form
leuk_AID[AID_date_list] <- lapply(leuk_AID[AID_date_list], as.Date)
leuk_AID <- leuk_AID %>%
  rowwise() %>%
  mutate(any_ai_disease_date = {
    vals <- c_across(all_of(AID_date_list))
    if (all(is.na(vals))) NA else min(vals, na.rm = TRUE)
  }) %>%
  ungroup()

#Any autoimmune disease except T1DM
leuk_AID$any_ai_disease_wo_T1DM <- as.integer(rowSums(leuk_AID[disease_list_wo_T1DM] == 1, na.rm = TRUE) > 0)

leuk_AID <- leuk_AID %>%
  rowwise() %>%
  mutate(any_ai_disease_date_wo_T1DM = {
    vals <- c_across(all_of(AID_date_list_wo_T1DM))
    if (all(is.na(vals))) NA else min(vals, na.rm = TRUE)
  }) %>%
  ungroup()

#Remove duplicates
leuk_AID <- leuk_AID %>%
  distinct(ID0, .keep_all = TRUE)

#Create variables for leukemia subtype
#AML
AML_icd = c("920", "C923", "C924", "C925", "C928", "C930", "C940", "C942")
AML_morpho = c(9840, 9861, 9866, 9867, 9873, 9874, 9891, 9910)

leuk_AID <- leuk_AID %>%
  mutate(
    AML = if_else(morpho %in% AML_morpho | iarccrgtools_icd10 %in% AML_icd, 1, 0)
  )

#ALL
ALL_icd = c("C910")
ALL_morpho = c(9811, 9812, 9816, 9820, 9835, 9836, 9837)

leuk_AID <- leuk_AID %>%
  mutate(
    ALL = if_else(morpho %in% ALL_morpho | iarccrgtools_icd10 %in% ALL_icd, 1, 0)
  )

# --- CREATE SUMMARY TABLE ---

#Create a list of all variables we want to include in the table
table1_vars <- c("ALL", "AML", "any_ai_disease_wo_T1DM","any_ai_disease", disease_list, "down", "panc",
                 "LGA", "TUPAKOINTITUNNUS")

#Get total counts for cases and controls
total_cases <- sum(leuk_AID$DV == 1)
total_controls <- sum(leuk_AID$DV == 0)

#Create summary table
table1 <- lapply(table1_vars, function(var) {
  case_count <- sum(leuk_AID[[var]] == 1 & leuk_AID$DV == 1,
                    na.rm = TRUE)
  control_count <- sum(leuk_AID[[var]] == 1 & leuk_AID$DV == 0,
                       na.rm = TRUE)
  data.frame(
    variable = var,
    cases_with_disease = case_count,
    cases_percent = round(case_count / total_cases * 100,
                          1),
    controls_with_disease = control_count,
    controls_percent = round(control_count / total_controls * 100, 1)
  )
})

#Combine into one data frame
table1_leuk_AID <- do.call(rbind, table1)

#Calculate age at reference date
leuk_AID$bd <- as.Date(leuk_AID$bd, format = "%Y-%m-%d")
str(leuk_AID)
leuk_AID <- leuk_AID %>%
  mutate(
    age_years = as.numeric(difftime(leuk_AID$ref.date, leuk_AID$bd, units =
                                      "days")) / 365.25,
    age_group = case_when(
      age_years >= 0 & age_years < 1 ~ "0-1",
      age_years >= 1 & age_years < 10 ~ "1-10",
      age_years >= 1.5 & age_years < 6 ~ "1.5-6",
      age_years >= 10 & age_years <= 18 ~ "10-18",
      TRUE ~ NA_character_
    ))

#Check distribution of sex among individuals diagnosed with AID
AID_cases <- leuk_AID[leuk_AID$any_ai_disease_wo_T1DM == 1, ]
AID_cases <- AID_cases[AID_cases$DV == 1, ]
table(AID_cases$sex) #48 male, 55 female

# --- EXCLUDING CONFOUNDING FACTORS ---
#Assess the number of Down syndrome cases among leukemia cases and controls
leuk_AID_down <- leuk_AID[which(leuk_AID$down==1), ]
table(leuk_AID_down$DV) #63 cases and 4 controls

#Removal of individuals diagnosed with Down syndrome
leuk_AID <- leuk_AID[-which(leuk_AID$down==1), ]

#Remove controls matched to cases with Down syndrome
leuk_AID_down_cases <- leuk_AID_down[leuk_AID_down$DV == 1, ] #Cases with Down syndrome
group_down <- leuk_AID_down_cases$group #Groups with cases diagnosed with Down syndrome
leuk_AID <- leuk_AID[!leuk_AID$group %in% group_down, ] #Remove controls matched to cases with Down

# --- CREATE PLOTS ---
#Create dataset for the plots
#Ensure reference date is in the date form
leuk_AID$ref.date <- as.Date(leuk_AID$ref.date)

#Create long format dataset
date_data <- leuk_AID %>%
  mutate(date_difference = as.numeric(any_ai_disease_date - ref.date)) %>%
  filter(leuk_AID$any_ai_disease == 1) %>%
  pivot_longer(
    cols = c(ref.date, all_of(AID_date_list)),
    names_to = "date_variable",
    values_to = "date"
  ) %>%
  select(ID0, DV, date_variable, date, date_difference) 


#Rename objects date variable for clarity
date_data <- date_data %>%
  mutate(date_variable = case_when(
    date_variable == "ref.date" ~ "Reference date",
    date_variable == "addisons_disease_date" ~ "Addison's disease",
    date_variable == "aiha_date" ~ "Autoimmune hemolytic anemia",
    date_variable == "AS_date" ~ "Ankylosing spondylitis",
    date_variable == "basedows_disease_date" ~ "Basedow's disease",
    date_variable == "celiac_disease_date" ~ "Celiac disease",
    date_variable == "colitis_ulcerosa_date" ~ "Colitis ulcerosa",
    date_variable == "crohns_disease_date" ~ "Crohn's disease",
    date_variable == "dermatomyocitis_date" ~ "Dermatomyocitis",
    date_variable == "ITP_date" ~ "Immune thrombocytopenia",
    date_variable == "juvenile_arthritis_date" ~ "Juvenile idiopathic arthritis",
    date_variable == "localized_scleroderma_date" ~ "Localized scleroderma",
    date_variable == "ms_disease_date" ~ "Multiple sclerosis (MS)",
    date_variable == "pemphigoid_date" ~ "Pemphigoid",
    date_variable == "pemphigus_date" ~ "Pemphigus",
    date_variable == "pernicious_anemia_date" ~ "Pernicious anemia",
    date_variable == "polyarteris_nodosa_date" ~ "Polyarteris nodosa",
    date_variable == "psoriasis_date" ~ "Psoriasis",
    date_variable == "RA_date" ~ "Rheumatoid arthritis",
    date_variable == "RF_date" ~ "Rheumatic fever",
    date_variable == "sarcoidosis_date" ~ "Sarcoidosis",
    date_variable == "sjogren_date" ~ "Sjögren's syndrome",
    date_variable == "SLE_date" ~ "Systemic lupus erythematosus",
    date_variable == "systemic_sclerosis_date" ~ "Systemic sclerosis",
    date_variable == "T1DM_date" ~ "Type 1 diabetes mellitus",
    date_variable == "thyroiditis_date" ~ "Autoimmune thyroiditis",
  ))

#Include only objects with diagnosis of autoimmune diseas
date_data <- date_data[!is.na(date_data$date), ] 

#Create a list of diseases that actually exist in the dataset
disease_list_plots <- unique(date_data$date_variable)

#Filter into two datasets
date_data_cases <- date_data %>%
  filter(date_data$DV == 1) #Dataset with leukemia cases

date_data_controls <- date_data %>%
  filter(date_data$DV == 0) #Dataset with controls matched to leukemia cases


# -- Creating the plot ---

#Plots for cases

#Apply signed log transformation, this keep negative values negative
#and compresses large values with base 5
date_data_cases$log_diff <- sign(date_data_cases$date_difference) *
  log10(abs(date_data_cases$date_difference) + 1)

#Sort the data by log_diff
date_data_cases <- date_data_cases[order(date_data_cases$log_diff), ]

#Reorder ID0 factor
date_data_cases$ID0 <- factor(date_data_cases$ID0, 
                              levels = 
                                unique(date_data_cases$ID0))


#Define custom breaks for the axis
breaks <- c(-4,-3,-2,-1,0,1,2,3,4)
labels <- c("-10,000", "-1,000", "-100", "-10", "0", "10", "100", "1,000",
            "10,000")

#Create the plot
p <- ggplot(date_data_cases) +
  geom_segment(aes(y = ID0, yend = ID0,
                   x = 0, xend = log_diff),
               color = "grey", linewidth = 0.3) +
  geom_point(aes(y = ID0, x = log_diff, color = "Autoimmune disease"),
             size = 1) + 
  scale_x_continuous(breaks = breaks, labels = labels) + 
  theme_minimal() +
  labs(
    title = "AID diagnosis date compared to leukemia diagnosis",
    x = "Difference in days",
    y = "Person",
    color = "Autoimmune disease"
  ) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

#Save the plot as a PDF
ggsave("plot_cases_all.pdf", plot = p, device = cairo_pdf(), width = 9, height = 6)

#Loop through each disease and create individual plots
for (disease in disease_list_plots) {
  #filter data for the current disease
  plot_data_cases <- date_data_cases[date_data_cases$date_variable == disease, ]
  
  #Skip if there is no data for this disease
  if(nrow(plot_data_cases) == 0) {
    message("Skipping", disease, "(no data)")
    next
  }
  #Sort the data by date_difference (ascending)
  plot_data_cases <- plot_data_cases[order(plot_data_cases$date_difference), ]
  
  #Reset factor levels for ID0 so ggplot respects the sorted order
  plot_data_cases$ID0 <- factor(plot_data_cases$ID0, levels = plot_data_cases$ID0)
  
  #Apply signed log transformation, this keep negative values negative
  #and compresses large values with base 5
  plot_data_cases$log_diff <- sign(plot_data_cases$date_difference) *
    log10(abs(plot_data_cases$date_difference) + 1)
  
  #Define custom breaks for the axis
  breaks <- c(-4,-3,-2,-1,0,1,2,3,4)
  labels <- c("-10,000", "-1,000", "-100", "-10", "0", "10", "100", "1,000",
              "10,000")
  
  #Create the plot
  p <- ggplot(plot_data_cases) +
    #Gray lines from reference (0) to the transformed value of AID 
    geom_segment(aes(
      y = ID0, yend = ID0,
      x = 0, xend = log_diff),
      color = "grey", linewidth = 0.3) +
    
    #Black poin at reference (0)
    geom_point(aes(y = ID0, x = 0), color = "black", size = 2) +
    
    #Colored points for AID differences
    geom_point(aes(y = ID0, x = log_diff,
                   color = date_variable), size = 2) +
    
    #Apply color palette
    scale_color_manual(values = color_palette) +
    
    #Custom x-axis with original day difference labels
    scale_x_continuous(breaks = breaks, 
                       labels = labels) +
    
    theme_minimal() +
    labs(
      title = paste(disease, "diagnosis date compared to leukemia diagnosis"),
      x = "Difference in days",
      y = "Person",
      color = "Autoimmune disease"
    ) + 
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  #Save the plot as a PDF
  file_name <- paste0("plot_cases", disease, ".pdf")
  ggsave(filename = file_name, plot = p, height = 6, width = 9)
}

#Repeat same for controls

#Apply signed log transformation, this keep negative values negative
#and compresses large values with base 5
date_data_controls$log_diff <- sign(date_data_controls$date_difference) *
  log10(abs(date_data_controls$date_difference) + 1)

#Sort the data by log_diff
date_data_controls <- date_data_controls[order(date_data_controls$log_diff), ]

#Reorder ID0 factor
date_data_controls$ID0 <- factor(date_data_controls$ID0, 
                                 levels = 
                                   unique(date_data_controls$ID0))


#Define custom breaks for the axis
breaks <- c(-4,-3,-2,-1,0,1,2,3,4)
labels <- c("-10,000", "-1,000", "-100", "-10", "0", "10", "100", "1,000",
            "10,000")

#Create the plot
p2 <- ggplot(date_data_controls) +
  geom_segment(aes(y = ID0, yend = ID0,
                   x = 0, xend = log_diff),
               color = "grey", linewidth = 0.3) +
  geom_point(aes(y = ID0, x = log_diff, color = "Autoimmune disease"),
             size = 1) + 
  scale_x_continuous(breaks = breaks, labels = labels) + 
  theme_minimal() +
  labs(
    title = "AID diagnosis date compared to control's reference date",
    x = "Difference in days",
    y = "Person",
    color = "Autoimmune disease"
  ) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

ggsave("plot_controls_all.pdf", plot = p2, device = cairo_pdf(), width = 9, height = 6)

#Loop through each disease and create plot
for (disease in disease_list_plots) {
  #filter data for the current disease
  plot_data_controls <- date_data_controls[date_data_controls$date_variable == disease, ]
  
  #Skip if there is no data for this disease
  if(nrow(plot_data_controls) == 0) {
    message("Skipping", disease, "(no data)")
    next
  }
  #Sort the data by date_difference (ascending)
  plot_data_controls <- plot_data_controls[order(plot_data_controls$date_difference), ]
  
  #Reset factor levels for ID0 so ggplot respects the sorted order
  plot_data_controls$ID0 <- factor(plot_data_controls$ID0, levels = plot_data_controls$ID0)
  
  #Apply signed log transformation, this keep negative values negative
  #and compresses large values with base 5
  plot_data_controls$log_diff <- sign(plot_data_controls$date_difference) *
    log10(abs(plot_data_controls$date_difference) + 1)
  
  #Define custom breaks for the axis
  breaks <- c(-4,-3,-2,-1,0,1,2,3,4)
  labels <- c("-10,000", "-1,000", "-100", "-10", "0", "10", "100", "1,000",
              "10,000")
  
  #Create the plot
  p <- ggplot(plot_data_controls) +
    #Gray lines from reference (0) to the transformed value of AID 
    geom_segment(aes(
      y = ID0, yend = ID0,
      x = 0, xend = log_diff),
      color = "grey", linewidth = 0.3) +
    
    #Black poin at reference (0)
    geom_point(aes(y = ID0, x = 0), color = "black", size = 2) +
    
    #Colored points for AID differences
    geom_point(aes(y = ID0, x = log_diff,
                   color = date_variable), size = 2) +
    
    #Apply color palette
    scale_color_manual(values = color_palette) +
    
    #Custom x-axis with original day difference labels
    scale_x_continuous(breaks = breaks, 
                       labels = labels) +
    
    theme_minimal() +
    labs(
      title = paste(disease, "diagnosis date compared to reference date"),
      x = "Difference in days",
      y = "Person",
      color = "Autoimmune disease"
    ) + 
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  #Save the plot as a PDF
  file_name <- paste0("plot_controls", disease, ".pdf")
  ggsave(filename = file_name, plot = p, height = 6, width = 9)
}

# --- RISK ANALYSES WITHOUT TYPE 1 DIABETES ---

# - CREATING SUBSETS -

#Stratified datasets based on sex
#1 = male
leuk_AID_male <- leuk_AID[leuk_AID$sex == 1, ]
table(leuk_AID_male$DV) #849 cases, 2543 controls

#2 = female
leuk_AID_female <- leuk_AID[leuk_AID$sex == 2, ]
table(leuk_AID_female$DV) #714 cases, 2141 controls

#Stratified by age at diagnosis
#0-0.99 years old
leuk_AID_0_1 <- leuk_AID[leuk_AID$age_group == "0-1", ]

#1-9.99 years old
leuk_AID_1_10 <- leuk_AID[leuk_AID$age_group == "1-10", ]

#10-17.99 years old
leuk_AID_10_18 <- leuk_AID[leuk_AID$age_group == "10-18", ]

#AML
leuk_AID_AML <- leuk_AID[leuk_AID$AML == 1, ]
groups_AML <- unique(leuk_AID_AML$group)
leuk_AID_AML <- leuk_AID[leuk_AID$group %in% groups_AML, ]

#ALL
leuk_AID_ALL <- leuk_AID[leuk_AID$ALL == 1, ]
groups_ALL <- unique(leuk_AID_ALL$group)
leuk_AID_ALL <- leuk_AID[leuk_AID$group %in% groups_ALL, ]

#Stratify ALL cases by age and sex
ALL_female <- leuk_AID_ALL[leuk_AID_ALL$sex == 2, ]
ALL_male <- leuk_AID_ALL[leuk_AID_ALL$sex == 1, ]
ALL_0_1 <- leuk_AID_ALL[leuk_AID_ALL$age_group == "0-1", ]
ALL_1.5_6 <- leuk_AID_ALL[leuk_AID_ALL$age_years >= 1.5 & leuk_AID_ALL$age_years <6, ]
ALL_1_10 <- leuk_AID_ALL[leuk_AID_ALL$age_group == "1-10", ]
ALL_10_18 <- leuk_AID_ALL[leuk_AID_ALL$age_group == "10-18", ]

#Create list of datasets and corresponding names
datasets = list(leuk_AID, leuk_AID_female, leuk_AID_male, leuk_AID_0_1, leuk_AID_1_10,
                leuk_AID_10_18, leuk_AID_ALL, ALL_female, ALL_male, ALL_0_1, ALL_1.5_6, ALL_1_10, ALL_10_18, 
                leuk_AID_AML)
names = list("all", "female", "male", "0-0.99", "1-9.99", "10-17.99", "ALL", "ALL_female", "ALL_male", "ALL 0-0.99", "ALL 1.5-5.99",
             "ALL 1-9.99", "ALL 10-17.99", "AML")

# -- UNIVARIATE MODEL --
# --- Conditional logistic regression model ---
res_leuk <- clogit(DV ~ any_ai_disease_wo_T1DM + strata(leuk_AID$group), data = leuk_AID)
summary(res_leuk)

#Define a function to exctract and format ORs and 95% CIs from a conditional logistic model
calculate_conf <- function(logmodel, name) {
  conf <- round(exp(confint(logmodel)), 2)
  summary_model <- summary(logmodel)
  odds_ratio <- round(exp(summary_model$coefficients[1, 1]), 2)
  temp_result <- c(name, as.character(odds_ratio),
                   as.character(conf[1, 1]), as.character(conf[1, 2]))
  return(temp_result)
}

#Any autoimmune disease (excluding T1DM) 

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk <- clogit(DV ~ any_ai_disease_wo_T1DM + strata(dataset$group), data = dataset)
  conf_results = rbind(conf_results, calculate_conf(res_leuk, name))
  
  #Model assuming AID precedes leukemia
  dataset$var <- dataset$any_ai_disease_wo_T1DM
  dataset$var[which(dataset$any_ai_disease_wo_T1DM == 1 &
                      (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order <- clogit(DV ~ var + strata(dataset$group), data = dataset)
  conf_results = rbind(conf_results, calculate_conf(res_leuk_order, paste0(name, "_AID_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2 <- dataset$any_ai_disease_wo_T1DM
  dataset$var2[which(dataset$any_ai_disease_wo_T1DM == 1 &
                       (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2 <- clogit(DV ~ var2 + strata(dataset$group), data = dataset)
  conf_results = rbind(conf_results, calculate_conf(res_leuk_order2, paste0(name, "_leuk_prec")))
  
  colnames(conf_results) <- column_names
}

#Inflammatory bowel diseases
#Create list of inflammatory bowel diseases
IBD_disease <- c("crohns_disease", "colitis_ulcerosa")

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_IBD = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Create variable for IBD
  dataset$IBD_var <- as.integer(rowSums(dataset[IBD_disease] == 1, na.rm = TRUE) > 0)
  
  #Basic model
  res_leuk_IBD <- clogit(DV ~ IBD_var + strata(dataset$group), data = dataset)
  conf_results_IBD = rbind(conf_results_IBD, calculate_conf(res_leuk_IBD, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_IBD <- dataset$IBD_var
  dataset$var_IBD[which(dataset$IBD_var == 1 &
                          (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_IBD <- clogit(DV ~ var_IBD + strata(dataset$group), data = dataset)
  conf_results_IBD = rbind(conf_results_IBD, calculate_conf(res_leuk_order_IBD, paste0(name, "_IBD_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_IBD <- dataset$IBD_var
  dataset$var2_IBD[which(dataset$IBD_var == 1 &
                           (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_IBD <- clogit(DV ~ var2_IBD + strata(dataset$group), data = dataset)
  conf_results_IBD = rbind(conf_results_IBD, calculate_conf(res_leuk_order2_IBD, paste0(name, "_leuk_prec")))
  
  colnames(conf_results_IBD) <- column_names
}

#Rheumatic disease
#Create list of rheumatic diseases
rheumatic_disease <- c("AS", "juvenile_arthritis", "RA")

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_rheumatic = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Create variable for rheumatic disease
  dataset$rheumatic_var <- as.integer(rowSums(dataset[rheumatic_disease] == 1, na.rm = TRUE) > 0)
  
  #Basic model
  res_leuk_rheumatic <- clogit(DV ~ rheumatic_var + strata(dataset$group), data = dataset)
  conf_results_rheumatic = rbind(conf_results_rheumatic, calculate_conf(res_leuk_rheumatic, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_rheumatic <- dataset$rheumatic_var
  dataset$var_rheumatic[which(dataset$rheumatic_var == 1 &
                                (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_rheumatic <- clogit(DV ~ var_rheumatic + strata(dataset$group), data = dataset)
  conf_results_rheumatic = rbind(conf_results_rheumatic, calculate_conf(res_leuk_order_rheumatic, paste0(name, "_rheumatic_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_rheumatic <- dataset$rheumatic_var
  dataset$var2_rheumatic[which(dataset$rheumatic_var == 1 &
                                 (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_rheumatic <- clogit(DV ~ var2_rheumatic + strata(dataset$group), data = dataset)
  conf_results_rheumatic = rbind(conf_results_rheumatic, calculate_conf(res_leuk_order2_rheumatic, paste0(name, "_leuk_prec")))
  
  colnames(conf_results_rheumatic) <- column_names
}

#Thyroid disease
#Create list of rheumatic diseases
thyroid_disease <- c("thyroiditis", "basedows_disease")

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_thyroid = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Create variable for thyroid disease
  dataset$thyroid_var <- as.integer(rowSums(dataset[thyroid_disease] == 1, na.rm = TRUE) > 0)
  
  #Basic model
  res_leuk_thyroid <- clogit(DV ~ thyroid_var + strata(dataset$group), data = dataset)
  conf_results_thyroid = rbind(conf_results_thyroid, calculate_conf(res_leuk_thyroid, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_thyroid <- dataset$thyroid_var
  dataset$var_thyroid[which(dataset$thyroid_var == 1 &
                              (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_thyroid <- clogit(DV ~ var_thyroid + strata(dataset$group), data = dataset)
  conf_results_thyroid = rbind(conf_results_thyroid, calculate_conf(res_leuk_order_thyroid, paste0(name, "_thyroid_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_thyroid <- dataset$thyroid_var
  dataset$var2_thyroid[which(dataset$thyroid_var == 1 &
                               (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_thyroid <- clogit(DV ~ var2_thyroid + strata(dataset$group), data = dataset)
  conf_results_thyroid = rbind(conf_results_thyroid, calculate_conf(res_leuk_order2_thyroid, paste0(name, "_leuk_prec")))
  
  colnames(conf_results_thyroid) <- column_names
}

#Sensitivity analyses

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_limited = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results) <- column_names


for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk <- clogit(DV ~ any_ai_disease_wo_T1DM + strata(dataset$group), data = dataset)
  conf_results_limited = rbind(conf_results_limited, calculate_conf(res_leuk, name))
  
  #Model excluding individuals diagnosed with AID under 90 days prior to leukemia diagnosis
  dataset$var3 <- dataset$any_ai_disease_wo_T1DM
  dataset$var3 <- dataset$any_ai_disease_wo_T1DM
  dataset$var3[which(dataset$any_ai_disease_wo_T1DM == 1 &
                       (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) < 0 &
                       (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) >= -90)] <- 0 
  res_leuk_lim <- clogit(DV ~ var3 + strata(dataset$group), data = dataset)
  conf_results_limited = rbind(conf_results_limited, calculate_conf(res_leuk_lim, paste0(name, "_AID_<30d_prec_leuk_excl")))
  
  
  #Model excluding individuals diagnosed with AID under one year (365 days) after leukemia diagnosis
  dataset$var4 <- dataset$any_ai_disease_wo_T1DM
  dataset$var4[which(dataset$any_ai_disease_wo_T1DM == 1 &
                       (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) > 0 &
                       (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) <= 365)] <- 0
  
  res_leuk_lim <- clogit(DV ~ var4 + strata(dataset$group), data = dataset)
  conf_results_limited = rbind(conf_results_limited, calculate_conf(res_leuk_lim, paste0(name, "_AID_<1y_after_leuk_excl")))
  
  #Model excluding individuals diagnosed with AID under 2,5 years (913 days) after leukemia diagnosis
  dataset$var5 <- dataset$any_ai_disease_wo_T1DM
  dataset$var5[which(dataset$any_ai_disease_wo_T1DM == 1 &
                       (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) > 0 &
                       (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) <= 913)] <- 0
  
  res_leuk_lim <- clogit(DV ~ var5 + strata(dataset$group), data = dataset)
  conf_results_limited = rbind(conf_results_limited, calculate_conf(res_leuk_lim, paste0(name, "_AID_<2.5y_after_leuk_excl")))
  
  #Model excluding individuals diagnosed with AID either AID under 90 days prior to or under 1 year (913 days) after leukemia diagnosis
  dataset$var6 <- dataset$any_ai_disease_wo_T1DM
  dataset$var6[which(dataset$any_ai_disease_wo_T1DM == 1 &
                       (((dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) < 0 &
                           (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) >= -90) |    
                          ((dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) > 0 &
                             (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) <= 365)))] <- 0
  
  res_leuk_lim <- clogit(DV ~ var6 + strata(dataset$group), data = dataset)
  conf_results_limited = rbind(conf_results_limited, calculate_conf(res_leuk_lim, paste0(name, "_AID_<90d_ prec_leuk_or_<1y_after_leuk_excl")))
  
  #Model excluding individuals diagnosed with AID either AID under 90 days prior to or under 2,5 years (913 days) after leukemia diagnosis
  dataset$var7 <- dataset$any_ai_disease_wo_T1DM
  dataset$var7[which(dataset$any_ai_disease_wo_T1DM == 1 &
                       (((dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) < 0 &
                           (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) >= -90) |    
                          ((dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) > 0 &
                             (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) <= 913)))] <- 0
  
  res_leuk_lim <- clogit(DV ~ var7 + strata(dataset$group), data = dataset)
  conf_results_limited = rbind(conf_results_limited, calculate_conf(res_leuk_lim, paste0(name, "_AID_<90d_ prec_leuk_or_<2.5y_after_leuk_excl")))
  
  colnames(conf_results_limited) <- column_names
}

#Special limitations: create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_limited_special = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results) <- column_names

for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Model excluding individuals diagnosed with ITP, Crohn's disease or colitis ulcerosa < 30 days apart from leukemia
  #or Addison's disease or JIA < 60 days apart from leukemia diagnosis
  #Or Rheumatoid arthritis < 90 days prior to leukemia diagnosis
  dataset$var8 <- dataset$any_ai_disease_wo_T1DM
  dataset$var8[which(
    (dataset$any_ai_disease_wo_T1DM == 1 &
       abs(dataset$ITP_date - dataset$ref.date) <= 30) |
      (dataset$any_ai_disease_wo_T1DM == 1 &
         abs(dataset$crohns_disease_date - dataset$ref.date) <= 30) |
      (dataset$any_ai_disease_wo_T1DM == 1 &
         abs(dataset$colitis_ulcerosa_date - dataset$ref.date) <= 30) |
      (dataset$any_ai_disease_wo_T1DM == 1 &
         abs(dataset$addisons_disease_date - dataset$ref.date) <= 60) |
      (dataset$any_ai_disease_wo_T1DM == 1 &
         abs(dataset$juvenile_arthritis_date - dataset$ref.date) <= 60) |
      (dataset$any_ai_disease_wo_T1DM == 1 &
         (dataset$RA_date - dataset$ref.date) < 0 &
         (dataset$RA_date - dataset$ref.date) >= -90)
  )] <- 0 
  
  res_leuk_lim_special <- clogit(DV ~ var8 + strata(dataset$group), data = dataset)
  conf_results_limited_special = rbind(conf_results_limited_special, calculate_conf(res_leuk_lim_special, paste0(name, "_special limitations")))
  
  colnames(conf_results_limited_special) <- column_names
}

#Individual autoimmune diseases

#Colitis ulcerosa 

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_UC = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results_UC) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_UC <- clogit(DV ~ colitis_ulcerosa + strata(dataset$group), data = dataset)
  conf_results_UC = rbind(conf_results_UC, calculate_conf(res_leuk_UC, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_UC <- dataset$colitis_ulcerosa
  dataset$var_UC[which(dataset$colitis_ulcerosa == 1 &
                         (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_UC <- clogit(DV ~ var_UC + strata(dataset$group), data = dataset)
  conf_results_UC = rbind(conf_results_UC, calculate_conf(res_leuk_order_UC, paste0(name, "_UC_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_UC <- dataset$colitis_ulcerosa
  dataset$var2_UC[which(dataset$colitis_ulcerosa == 1 &
                          (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_UC <- clogit(DV ~ var2_UC + strata(dataset$group), data = dataset)
  conf_results_UC = rbind(conf_results_UC, calculate_conf(res_leuk_order2_UC, paste0(name, "_leuk_prec")))
  
  colnames(conf_results_UC) <- column_names
}

#Celiac disease

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_celiac = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results_celiac) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_celiac <- clogit(DV ~ celiac_disease + strata(dataset$group), data = dataset)
  conf_results_celiac = rbind(conf_results_celiac, calculate_conf(res_leuk_celiac, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_celiac <- dataset$celiac_disease
  dataset$var_celiac[which(dataset$celiac_disease == 1 &
                             (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_celiac <- clogit(DV ~ var_celiac + strata(dataset$group), data = dataset)
  conf_results_celiac = rbind(conf_results_celiac, calculate_conf(res_leuk_order_celiac, paste0(name, "_celiac_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_celiac <- dataset$celiac_disease
  dataset$var2_celiac[which(dataset$celiac_disease == 1 &
                              (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_celiac <- clogit(DV ~ var2_celiac + strata(dataset$group), data = dataset)
  conf_results_celiac = rbind(conf_results_celiac, calculate_conf(res_leuk_order2_celiac, paste0(name, "_leuk_prec")))
  
  colnames(conf_results_celiac) <- column_names
}

#Crohns disease

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_crohns = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results_crohns) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_crohns <- clogit(DV ~ crohns_disease + strata(dataset$group), data = dataset)
  conf_results_crohns = rbind(conf_results_crohns, calculate_conf(res_leuk_crohns, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_crohns <- dataset$crohns_disease
  dataset$var_crohns[which(dataset$crohns_disease == 1 &
                             (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_crohns <- clogit(DV ~ var_crohns + strata(dataset$group), data = dataset)
  conf_results_crohns = rbind(conf_results_crohns, calculate_conf(res_leuk_order_crohns, paste0(name, "_crohns_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_crohns <- dataset$crohns_disease
  dataset$var2_crohns[which(dataset$crohns_disease == 1 &
                              (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_crohns <- clogit(DV ~ var2_crohns + strata(dataset$group), data = dataset)
  conf_results_crohns = rbind(conf_results_crohns, calculate_conf(res_leuk_order2_crohns, paste0(name, "_leuk_prec")))
  
  colnames(conf_results_crohns) <- column_names
}

#Addison's disease

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_addisons = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results_addisons) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_addisons <- clogit(DV ~ addisons_disease + strata(dataset$group), data = dataset)
  conf_results_addisons = rbind(conf_results_addisons, calculate_conf(res_leuk_addisons, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_addisons <- dataset$addisons_disease
  dataset$var_addisons[which(dataset$addisons_disease == 1 &
                               (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_addisons <- clogit(DV ~ var_addisons + strata(dataset$group), data = dataset)
  conf_results_addisons = rbind(conf_results_addisons, calculate_conf(res_leuk_order_addisons, paste0(name, "_addisons_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_addisons <- dataset$addisons_disease
  dataset$var2_addisons[which(dataset$addisons_disease == 1 &
                                (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_addisons <- clogit(DV ~ var2_addisons + strata(dataset$group), data = dataset)
  conf_results_addisons = rbind(conf_results_addisons, calculate_conf(res_leuk_order2_addisons, paste0(name, "_leuk_prec")))
  
  colnames(conf_results_addisons) <- column_names
}

#Psoriasis

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_psoriasis = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results_psoriasis) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_psoriasis <- clogit(DV ~ psoriasis + strata(dataset$group), data = dataset)
  conf_results_psoriasis = rbind(conf_results_psoriasis, calculate_conf(res_leuk_psoriasis, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_psoriasis <- dataset$psoriasis
  dataset$var_psoriasis[which(dataset$psoriasis == 1 &
                                (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_psoriasis <- clogit(DV ~ var_psoriasis + strata(dataset$group), data = dataset)
  conf_results_psoriasis = rbind(conf_results_psoriasis, calculate_conf(res_leuk_order_psoriasis, paste0(name, "_psoriasis_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_psoriasis <- dataset$psoriasis
  dataset$var2_psoriasis[which(dataset$psoriasis == 1 &
                                 (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_psoriasis <- clogit(DV ~ var2_psoriasis + strata(dataset$group), data = dataset)
  conf_results_psoriasis = rbind(conf_results_psoriasis, calculate_conf(res_leuk_order2_psoriasis, paste0(name, "_leuk_prec")))
  
  colnames(conf_results_psoriasis) <- column_names
}

#Juvenile idiopathic arthritis

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_juvenile_arthritis = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results_juvenile_arthritis) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_juvenile_arthritis <- clogit(DV ~ juvenile_arthritis + strata(dataset$group), data = dataset)
  conf_results_juvenile_arthritis = rbind(conf_results_juvenile_arthritis, calculate_conf(res_leuk_juvenile_arthritis, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_juvenile_arthritis <- dataset$juvenile_arthritis
  dataset$var_juvenile_arthritis[which(dataset$juvenile_arthritis == 1 &
                                         (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_juvenile_arthritis <- clogit(DV ~ var_juvenile_arthritis + strata(dataset$group), data = dataset)
  conf_results_juvenile_arthritis = rbind(conf_results_juvenile_arthritis, calculate_conf(res_leuk_order_juvenile_arthritis, paste0(name, "_juvenile_arthritis_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_juvenile_arthritis <- dataset$juvenile_arthritis
  dataset$var2_juvenile_arthritis[which(dataset$juvenile_arthritis == 1 &
                                          (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_juvenile_arthritis <- clogit(DV ~ var2_juvenile_arthritis + strata(dataset$group), data = dataset)
  conf_results_juvenile_arthritis = rbind(conf_results_juvenile_arthritis, calculate_conf(res_leuk_order2_juvenile_arthritis, paste0(name, "_leuk_prec")))
  
  colnames(conf_results_juvenile_arthritis) <- column_names
}


#Rheumatoid arthritis

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_RA = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results_RA) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_RA <- clogit(DV ~ RA + strata(dataset$group), data = dataset)
  conf_results_RA = rbind(conf_results_RA, calculate_conf(res_leuk_RA, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_RA <- dataset$RA
  dataset$var_RA[which(dataset$RA == 1 &
                         (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_RA <- clogit(DV ~ var_RA + strata(dataset$group), data = dataset)
  conf_results_RA = rbind(conf_results_RA, calculate_conf(res_leuk_order_RA, paste0(name, "_RA_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_RA <- dataset$RA
  dataset$var2_RA[which(dataset$RA == 1 &
                          (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_RA <- clogit(DV ~ var2_RA + strata(dataset$group), data = dataset)
  conf_results_RA = rbind(conf_results_RA, calculate_conf(res_leuk_order2_RA, paste0(name, "_leuk_prec")))
  
  colnames(conf_results_RA) <- column_names
}

#Ankylosing spondylitis

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_AS = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results_AS) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_AS <- clogit(DV ~ AS + strata(dataset$group), data = dataset)
  conf_results_AS = rbind(conf_results_AS, calculate_conf(res_leuk_AS, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_AS <- dataset$AS
  dataset$var_AS[which(dataset$AS == 1 &
                         (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_AS <- clogit(DV ~ var_AS + strata(dataset$group), data = dataset)
  conf_results_AS = rbind(conf_results_AS, calculate_conf(res_leuk_order_AS, paste0(name, "_AS_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_AS <- dataset$AS
  dataset$var2_AS[which(dataset$AS == 1 &
                          (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_AS <- clogit(DV ~ var2_AS + strata(dataset$group), data = dataset)
  conf_results_AS = rbind(conf_results_AS, calculate_conf(res_leuk_order2_AS, paste0(name, "_leuk_prec")))
  
  colnames(conf_results_AS) <- column_names
}

#Immune thrombocytopenia (ITP)

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_ITP = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results_ITP) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_ITP <- clogit(DV ~ ITP + strata(dataset$group), data = dataset)
  conf_results_ITP = rbind(conf_results_ITP, calculate_conf(res_leuk_ITP, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_ITP <- dataset$ITP
  dataset$var_ITP[which(dataset$ITP == 1 &
                          (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_ITP <- clogit(DV ~ var_ITP + strata(dataset$group), data = dataset)
  conf_results_ITP = rbind(conf_results_ITP, calculate_conf(res_leuk_order_ITP, paste0(name, "_ITP_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_ITP <- dataset$ITP
  dataset$var2_ITP[which(dataset$ITP == 1 &
                           (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_ITP <- clogit(DV ~ var2_ITP + strata(dataset$group), data = dataset)
  conf_results_ITP = rbind(conf_results_ITP, calculate_conf(res_leuk_order2_ITP, paste0(name, "_leuk_prec")))
  
  colnames(conf_results_ITP) <- column_names
}


# -- MULTIVARIATE MODEL --
#Based on our previous work (Ventelä et al 2025) LGA and smoking were confounding factors
#Information from the Medical Birth Register

#Define a function to exctract and format ORs and 95% CIs from a conditional logistic model
calculate_multi <- function(logmodel, name) {
  conf <- round(exp(confint(logmodel)), 2)
  summary_model <- summary(logmodel)
  odds_ratio <- round(exp(summary_model$coefficients[1, 1]), 2)
  temp_result <- c(name, as.character(odds_ratio),
                   as.character(conf[1, 1]), as.character(conf[1, 2]))
  return(temp_result)
}

#Use same list of datasets and names


#Multivariable model
res_leuk_adj <- clogit(DV ~ any_ai_disease_wo_T1DM + LGA + TUPAKOINTITUNNUS
                       + strata(leuk_AID$group), data = leuk_AID)
summary(res_leuk_adj) #OR 2.05, 95% 1.48-2.84

#Any autoimmune disease (excluding T1DM) 

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk <- clogit(DV ~ any_ai_disease_wo_T1DM + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results = rbind(multi_results, calculate_multi(res_leuk, name))
  
  #Model assuming AID precedes leukemia
  dataset$var <- dataset$any_ai_disease_wo_T1DM
  dataset$var[which(dataset$any_ai_disease_wo_T1DM == 1 &
                      (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order <- clogit(DV ~ var + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results = rbind(multi_results, calculate_multi(res_leuk_order, paste0(name, "_AID_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2 <- dataset$any_ai_disease_wo_T1DM
  dataset$var2[which(dataset$any_ai_disease_wo_T1DM == 1 &
                       (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2 <- clogit(DV ~ var2 + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results = rbind(multi_results, calculate_multi(res_leuk_order2, paste0(name, "_leuk_prec")))
  
  colnames(multi_results) <- column_names
}

#Inflammatory bowel diseases
#Create list of inflammatory bowel diseases
IBD_disease <- c("crohns_disease", "colitis_ulcerosa")

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_IBD = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Create variable for IBD
  dataset$IBD_var <- as.integer(rowSums(dataset[IBD_disease] == 1, na.rm = TRUE) > 0)
  
  #Basic model
  res_leuk_IBD <- clogit(DV ~ IBD_var + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_IBD = rbind(multi_results_IBD, calculate_multi(res_leuk_IBD, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_IBD <- dataset$IBD_var
  dataset$var_IBD[which(dataset$IBD_var == 1 &
                          (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_IBD <- clogit(DV ~ var_IBD + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_IBD = rbind(multi_results_IBD, calculate_multi(res_leuk_order_IBD, paste0(name, "_IBD_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_IBD <- dataset$IBD_var
  dataset$var2_IBD[which(dataset$IBD_var == 1 &
                           (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_IBD <- clogit(DV ~ var2_IBD + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_IBD = rbind(multi_results_IBD, calculate_multi(res_leuk_order2_IBD, paste0(name, "_leuk_prec")))
  
  colnames(multi_results_IBD) <- column_names
}

#Rheumatic disease
#Create list of rheumatic diseases
rheumatic_disease <- c("AS", "juvenile_arthritis", "RA")

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_rheumatic = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Create variable for rheumatic disease
  dataset$rheumatic_var <- as.integer(rowSums(dataset[rheumatic_disease] == 1, na.rm = TRUE) > 0)
  
  #Basic model
  res_leuk_rheumatic <- clogit(DV ~ rheumatic_var + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_rheumatic = rbind(multi_results_rheumatic, calculate_multi(res_leuk_rheumatic, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_rheumatic <- dataset$rheumatic_var
  dataset$var_rheumatic[which(dataset$rheumatic_var == 1 &
                                (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_rheumatic <- clogit(DV ~ var_rheumatic + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_rheumatic = rbind(multi_results_rheumatic, calculate_multi(res_leuk_order_rheumatic, paste0(name, "_rheumatic_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_rheumatic <- dataset$rheumatic_var
  dataset$var2_rheumatic[which(dataset$rheumatic_var == 1 &
                                 (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_rheumatic <- clogit(DV ~ var2_rheumatic + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_rheumatic = rbind(multi_results_rheumatic, calculate_multi(res_leuk_order2_rheumatic, paste0(name, "_leuk_prec")))
  
  colnames(multi_results_rheumatic) <- column_names
}

#Thyroid disease
#Create list of rheumatic diseases
thyroid_disease <- c("thyroiditis", "basedows_disease")

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_thyroid = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Create variable for thyroid disease
  dataset$thyroid_var <- as.integer(rowSums(dataset[thyroid_disease] == 1, na.rm = TRUE) > 0)
  
  #Basic model
  res_leuk_thyroid <- clogit(DV ~ thyroid_var + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_thyroid = rbind(multi_results_thyroid, calculate_multi(res_leuk_thyroid, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_thyroid <- dataset$thyroid_var
  dataset$var_thyroid[which(dataset$thyroid_var == 1 &
                              (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_thyroid <- clogit(DV ~ var_thyroid + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_thyroid = rbind(multi_results_thyroid, calculate_multi(res_leuk_order_thyroid, paste0(name, "_thyroid_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_thyroid <- dataset$thyroid_var
  dataset$var2_thyroid[which(dataset$thyroid_var == 1 &
                               (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_thyroid <- clogit(DV ~ var2_thyroid + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_thyroid = rbind(multi_results_thyroid, calculate_multi(res_leuk_order2_thyroid, paste0(name, "_leuk_prec")))
  
  colnames(multi_results_thyroid) <- column_names
}

#Sensitivity analyses

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_limited = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results) <- column_names


for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk <- clogit(DV ~ any_ai_disease_wo_T1DM + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_limited = rbind(multi_results_limited, calculate_multi(res_leuk, name))
  
  #Model excluding individuals diagnosed with AID under 90 days prior to leukemia diagnosis
  dataset$var3 <- dataset$any_ai_disease_wo_T1DM
  dataset$var3 <- dataset$any_ai_disease_wo_T1DM
  dataset$var3[which(dataset$any_ai_disease_wo_T1DM == 1 &
                       (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) < 0 &
                       (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) >= -90)] <- 0 
  res_leuk_lim <- clogit(DV ~ var3 + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_limited = rbind(multi_results_limited, calculate_multi(res_leuk_lim, paste0(name, "_AID_<30d_prec_leuk_excl")))
  
  
  #Model excluding individuals diagnosed with AID under one year (365 days) after leukemia diagnosis
  dataset$var4 <- dataset$any_ai_disease_wo_T1DM
  dataset$var4[which(dataset$any_ai_disease_wo_T1DM == 1 &
                       (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) > 0 &
                       (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) <= 365)] <- 0
  
  res_leuk_lim <- clogit(DV ~ var4 + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_limited = rbind(multi_results_limited, calculate_multi(res_leuk_lim, paste0(name, "_AID_<1y_after_leuk_excl")))
  
  #Model excluding individuals diagnosed with AID under 2,5 years (913 days) after leukemia diagnosis
  dataset$var5 <- dataset$any_ai_disease_wo_T1DM
  dataset$var5[which(dataset$any_ai_disease_wo_T1DM == 1 &
                       (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) > 0 &
                       (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) <= 913)] <- 0
  
  res_leuk_lim <- clogit(DV ~ var5 + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_limited = rbind(multi_results_limited, calculate_multi(res_leuk_lim, paste0(name, "_AID_<2.5y_after_leuk_excl")))
  
  #Model excluding individuals diagnosed with AID either AID under 90 days prior to or under 1 year (913 days) after leukemia diagnosis
  dataset$var6 <- dataset$any_ai_disease_wo_T1DM
  dataset$var6[which(dataset$any_ai_disease_wo_T1DM == 1 &
                       (((dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) < 0 &
                           (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) >= -90) |    
                          ((dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) > 0 &
                             (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) <= 365)))] <- 0
  
  res_leuk_lim <- clogit(DV ~ var6 + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_limited = rbind(multi_results_limited, calculate_multi(res_leuk_lim, paste0(name, "_AID_<90d_ prec_leuk_or_<1y_after_leuk_excl")))
  
  #Model excluding individuals diagnosed with AID either AID under 90 days prior to or under 2,5 years (913 days) after leukemia diagnosis
  dataset$var7 <- dataset$any_ai_disease_wo_T1DM
  dataset$var7[which(dataset$any_ai_disease_wo_T1DM == 1 &
                       (((dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) < 0 &
                           (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) >= -90) |    
                          ((dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) > 0 &
                             (dataset$any_ai_disease_date_wo_T1DM - dataset$ref.date) <= 913)))] <- 0
  
  res_leuk_lim <- clogit(DV ~ var7 + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_limited = rbind(multi_results_limited, calculate_multi(res_leuk_lim, paste0(name, "_AID_<90d_ prec_leuk_or_<2.5y_after_leuk_excl")))
  
  colnames(multi_results_limited) <- column_names
}

#Special limitations: create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_limited_special = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results) <- column_names

for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  
  #Model excluding individuals diagnosed with ITP, Crohn's disease or colitis ulcerosa < 30 days apart from leukemia
  #or Addison's disease or JIA < 60 days apart from leukemia diagnosis
  #Or Rheumatoid arthritis < 90 days prior to leukemia diagnosis
  dataset$var8 <- dataset$any_ai_disease_wo_T1DM
  dataset$var8[which(
    (dataset$any_ai_disease_wo_T1DM == 1 &
       abs(dataset$ITP_date - dataset$ref.date) <= 30) |
      (dataset$any_ai_disease_wo_T1DM == 1 &
         abs(dataset$crohns_disease_date - dataset$ref.date) <= 30) |
      (dataset$any_ai_disease_wo_T1DM == 1 &
         abs(dataset$colitis_ulcerosa_date - dataset$ref.date) <= 30) |
      (dataset$any_ai_disease_wo_T1DM == 1 &
         abs(dataset$addisons_disease_date - dataset$ref.date) <= 60) |
      (dataset$any_ai_disease_wo_T1DM == 1 &
         abs(dataset$juvenile_arthritis_date - dataset$ref.date) <= 60) |
      (dataset$any_ai_disease_wo_T1DM == 1 &
         (dataset$RA_date - dataset$ref.date) < 0 &
         (dataset$RA_date - dataset$ref.date) >= -90)
  )] <- 0 
  
  res_leuk_lim_special <- clogit(DV ~ var8 + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_limited_special = rbind(multi_results_limited_special, calculate_multi(res_leuk_lim_special, paste0(name, "_special limitations")))
  
  colnames(multi_results_limited_special) <- column_names
}

#Individual autoimmune diseases

#Colitis ulcerosa 

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_UC = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results_UC) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_UC <- clogit(DV ~ colitis_ulcerosa + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_UC = rbind(multi_results_UC, calculate_multi(res_leuk_UC, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_UC <- dataset$colitis_ulcerosa
  dataset$var_UC[which(dataset$colitis_ulcerosa == 1 &
                         (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_UC <- clogit(DV ~ var_UC + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_UC = rbind(multi_results_UC, calculate_multi(res_leuk_order_UC, paste0(name, "_UC_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_UC <- dataset$colitis_ulcerosa
  dataset$var2_UC[which(dataset$colitis_ulcerosa == 1 &
                          (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_UC <- clogit(DV ~ var2_UC + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_UC = rbind(multi_results_UC, calculate_multi(res_leuk_order2_UC, paste0(name, "_leuk_prec")))
  
  colnames(multi_results_UC) <- column_names
}

#Celiac disease

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_celiac = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results_celiac) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_celiac <- clogit(DV ~ celiac_disease + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_celiac = rbind(multi_results_celiac, calculate_multi(res_leuk_celiac, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_celiac <- dataset$celiac_disease
  dataset$var_celiac[which(dataset$celiac_disease == 1 &
                             (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_celiac <- clogit(DV ~ var_celiac + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_celiac = rbind(multi_results_celiac, calculate_multi(res_leuk_order_celiac, paste0(name, "_celiac_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_celiac <- dataset$celiac_disease
  dataset$var2_celiac[which(dataset$celiac_disease == 1 &
                              (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_celiac <- clogit(DV ~ var2_celiac + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_celiac = rbind(multi_results_celiac, calculate_multi(res_leuk_order2_celiac, paste0(name, "_leuk_prec")))
  
  colnames(multi_results_celiac) <- column_names
}

#Crohns disease

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_crohns = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results_crohns) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_crohns <- clogit(DV ~ crohns_disease + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_crohns = rbind(multi_results_crohns, calculate_multi(res_leuk_crohns, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_crohns <- dataset$crohns_disease
  dataset$var_crohns[which(dataset$crohns_disease == 1 &
                             (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_crohns <- clogit(DV ~ var_crohns + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_crohns = rbind(multi_results_crohns, calculate_multi(res_leuk_order_crohns, paste0(name, "_crohns_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_crohns <- dataset$crohns_disease
  dataset$var2_crohns[which(dataset$crohns_disease == 1 &
                              (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_crohns <- clogit(DV ~ var2_crohns + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_crohns = rbind(multi_results_crohns, calculate_multi(res_leuk_order2_crohns, paste0(name, "_leuk_prec")))
  
  colnames(multi_results_crohns) <- column_names
}

#Addison's disease

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_addisons = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results_addisons) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_addisons <- clogit(DV ~ addisons_disease + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_addisons = rbind(multi_results_addisons, calculate_multi(res_leuk_addisons, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_addisons <- dataset$addisons_disease
  dataset$var_addisons[which(dataset$addisons_disease == 1 &
                               (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_addisons <- clogit(DV ~ var_addisons + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_addisons = rbind(multi_results_addisons, calculate_multi(res_leuk_order_addisons, paste0(name, "_addisons_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_addisons <- dataset$addisons_disease
  dataset$var2_addisons[which(dataset$addisons_disease == 1 &
                                (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_addisons <- clogit(DV ~ var2_addisons + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_addisons = rbind(multi_results_addisons, calculate_multi(res_leuk_order2_addisons, paste0(name, "_leuk_prec")))
  
  colnames(multi_results_addisons) <- column_names
}

#Psoriasis

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_psoriasis = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results_psoriasis) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_psoriasis <- clogit(DV ~ psoriasis + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_psoriasis = rbind(multi_results_psoriasis, calculate_multi(res_leuk_psoriasis, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_psoriasis <- dataset$psoriasis
  dataset$var_psoriasis[which(dataset$psoriasis == 1 &
                                (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_psoriasis <- clogit(DV ~ var_psoriasis + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_psoriasis = rbind(multi_results_psoriasis, calculate_multi(res_leuk_order_psoriasis, paste0(name, "_psoriasis_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_psoriasis <- dataset$psoriasis
  dataset$var2_psoriasis[which(dataset$psoriasis == 1 &
                                 (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_psoriasis <- clogit(DV ~ var2_psoriasis + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_psoriasis = rbind(multi_results_psoriasis, calculate_multi(res_leuk_order2_psoriasis, paste0(name, "_leuk_prec")))
  
  colnames(multi_results_psoriasis) <- column_names
}

#Juvenile idiopathic arthritis

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_juvenile_arthritis = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results_juvenile_arthritis) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_juvenile_arthritis <- clogit(DV ~ juvenile_arthritis + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_juvenile_arthritis = rbind(multi_results_juvenile_arthritis, calculate_multi(res_leuk_juvenile_arthritis, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_juvenile_arthritis <- dataset$juvenile_arthritis
  dataset$var_juvenile_arthritis[which(dataset$juvenile_arthritis == 1 &
                                         (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_juvenile_arthritis <- clogit(DV ~ var_juvenile_arthritis + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_juvenile_arthritis = rbind(multi_results_juvenile_arthritis, calculate_multi(res_leuk_order_juvenile_arthritis, paste0(name, "_juvenile_arthritis_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_juvenile_arthritis <- dataset$juvenile_arthritis
  dataset$var2_juvenile_arthritis[which(dataset$juvenile_arthritis == 1 &
                                          (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_juvenile_arthritis <- clogit(DV ~ var2_juvenile_arthritis + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_juvenile_arthritis = rbind(multi_results_juvenile_arthritis, calculate_multi(res_leuk_order2_juvenile_arthritis, paste0(name, "_leuk_prec")))
  
  colnames(multi_results_juvenile_arthritis) <- column_names
}


#Rheumatoid arthritis

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_RA = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results_RA) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_RA <- clogit(DV ~ RA + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_RA = rbind(multi_results_RA, calculate_multi(res_leuk_RA, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_RA <- dataset$RA
  dataset$var_RA[which(dataset$RA == 1 &
                         (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_RA <- clogit(DV ~ var_RA + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_RA = rbind(multi_results_RA, calculate_multi(res_leuk_order_RA, paste0(name, "_RA_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_RA <- dataset$RA
  dataset$var2_RA[which(dataset$RA == 1 &
                          (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_RA <- clogit(DV ~ var2_RA + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_RA = rbind(multi_results_RA, calculate_multi(res_leuk_order2_RA, paste0(name, "_leuk_prec")))
  
  colnames(multi_results_RA) <- column_names
}

#Ankylosing spondylitis

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_AS = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results_AS) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_AS <- clogit(DV ~ AS + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_AS = rbind(multi_results_AS, calculate_multi(res_leuk_AS, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_AS <- dataset$AS
  dataset$var_AS[which(dataset$AS == 1 &
                         (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_AS <- clogit(DV ~ var_AS + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_AS = rbind(multi_results_AS, calculate_multi(res_leuk_order_AS, paste0(name, "_AS_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_AS <- dataset$AS
  dataset$var2_AS[which(dataset$AS == 1 &
                          (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_AS <- clogit(DV ~ var2_AS + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_AS = rbind(multi_results_AS, calculate_multi(res_leuk_order2_AS, paste0(name, "_leuk_prec")))
  
  colnames(multi_results_AS) <- column_names
}

#Immune thrombocytopenia (ITP)

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_ITP = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results_ITP) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets)) {
  dataset = datasets[[i]]
  name = names[[i]]
  #Basic model
  res_leuk_ITP <- clogit(DV ~ ITP + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_ITP = rbind(multi_results_ITP, calculate_multi(res_leuk_ITP, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_ITP <- dataset$ITP
  dataset$var_ITP[which(dataset$ITP == 1 &
                          (dataset$ref.date < dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order_ITP <- clogit(DV ~ var_ITP + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_ITP = rbind(multi_results_ITP, calculate_multi(res_leuk_order_ITP, paste0(name, "_ITP_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_ITP <- dataset$ITP
  dataset$var2_ITP[which(dataset$ITP == 1 &
                           (dataset$ref.date > dataset$any_ai_disease_date_wo_T1DM))] <- 0
  res_leuk_order2_ITP <- clogit(DV ~ var2_ITP + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_ITP = rbind(multi_results_ITP, calculate_multi(res_leuk_order2_ITP, paste0(name, "_leuk_prec")))
  
  colnames(multi_results_ITP) <- column_names
}

# --- RISK ANALYSES INCLUDING TYPE 1 DIABETES ---
#Create dataset, same dataset as leuk_AID
leuk_pooled <- leuk_AID

#Removal of individuals with a history of pancreatitis
#Assess the number of pancreatitis cases among leukemia cases and controls
leuk_panc <- leuk_pooled[which(leuk_pooled$panc == 1), ]
table(leuk_panc$DV) #30 cases, 2 controls

#Remove individuals with pancreatitis
leuk_pooled <- leuk_pooled[-which(leuk_pooled$panc==1), ]

#Remove controls matched to cases with pancreatitis
leuk_panc_cases <- leuk_panc[leuk_panc$DV == 1, ]
group_panc <- leuk_panc_cases$group
leuk_pooled <- leuk_pooled[!leuk_pooled$group %in% group_panc, ]

# - CREATE SUBSETS -

#1 = male
pooled_male <- leuk_pooled[leuk_pooled$sex == 1, ]
table(pooled_male$DV) #833 cases, 2493 controls

#2 = female
pooled_female <- leuk_pooled[leuk_pooled$sex == 2, ]
table(pooled_female$DV) #700 cases, 2100 controls

#Stratified by age at diagnosis
#0-0.99 years old
pooled_0_1 <- leuk_pooled[leuk_pooled$age_group == "0-1", ]

#1-9.99 years old
pooled_1_10 <- leuk_pooled[leuk_pooled$age_group == "1-10", ]

#10-17.99 years old
pooled_10_18 <- leuk_pooled[leuk_pooled$age_group == "10-18", ]

#AML
pooled_AML <- leuk_pooled[leuk_pooled$AML == 1, ]
groups_AML <- unique(pooled_AML$group)
pooled_AML <- leuk_pooled[leuk_pooled$group %in% groups_AML, ]

#ALL
pooled_ALL <- leuk_pooled[leuk_pooled$ALL == 1, ]
groups_ALL <- unique(pooled_ALL$group)
pooled_ALL <- leuk_pooled[leuk_pooled$group %in% groups_ALL, ]

#Stratify ALL cases by age and sex
pooled_ALL_female <- pooled_ALL[pooled_ALL$sex == 2, ]
pooled_ALL_male <- pooled_ALL[pooled_ALL$sex == 1, ]
pooled_ALL_0_1 <- pooled_ALL[pooled_ALL$age_group == "0-1", ]
pooled_ALL_1.5_6 <- pooled_ALL[pooled_ALL$age_years >= 1.5 & pooled_ALL$age_years <6, ]
pooled_ALL_1_10 <- pooled_ALL[pooled_ALL$age_group == "1-10", ]
pooled_ALL_10_18 <- pooled_ALL[pooled_ALL$age_group == "10-18", ]

#Create list of datasets and corresponding names
datasets_pooled = list(leuk_pooled, pooled_female, pooled_male, pooled_0_1, pooled_1_10,
                       pooled_10_18, pooled_ALL, pooled_ALL_female, pooled_ALL_male, pooled_ALL_0_1,
                       pooled_ALL_1.5_6, pooled_ALL_1_10, pooled_ALL_10_18, 
                       pooled_AML)
names_pooled = list("all", "female", "male", "0-0.99", "1-9.99", "10-17.99", "ALL", "ALL_female", "ALL_male", "ALL 0-0.99", "ALL 1.5-5.99",
                    "ALL 1-9.99", "ALL 10-17.99", "AML")

# -- UNIVARIATE MODEL --
# --- Conditional logistic regression model ---
res_leuk_pooled <- clogit(DV ~ any_ai_disease + strata(leuk_pooled$group), data = leuk_pooled)
summary(res_leuk_pooled) #OR=1.80 95% CI: 1.40-2.30

#Define a function to exctract and format ORs and 95% CIs from a conditional logistic model
calculate_conf <- function(logmodel, name) {
  conf <- round(exp(confint(logmodel)), 2)
  summary_model <- summary(logmodel)
  odds_ratio <- round(exp(summary_model$coefficients[1, 1]), 2)
  temp_result <- c(name, as.character(odds_ratio),
                   as.character(conf[1, 1]), as.character(conf[1, 2]))
  return(temp_result)
}


# -- UNIVARIATE MODEL --
#Any autoimmune disease (excluding T1DM) 

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
conf_results_pooled = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(conf_results_pooled) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets_pooled)) {
  dataset = datasets_pooled[[i]]
  name = names_pooled[[i]]
  #Basic model
  res_leuk_pooled <- clogit(DV ~ any_ai_disease + strata(dataset$group), data = dataset)
  conf_results_pooled = rbind(conf_results_pooled, calculate_conf(res_leuk_pooled, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_pooled <- dataset$any_ai_disease
  dataset$var_pooled[which(dataset$any_ai_disease == 1 &
                             (dataset$ref.date < dataset$any_ai_disease_date))] <- 0
  res_leuk_order_pooled <- clogit(DV ~ var_pooled + strata(dataset$group), data = dataset)
  conf_results_pooled = rbind(conf_results_pooled, calculate_conf(res_leuk_order_pooled, paste0(name, "_AID_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_pooled <- dataset$any_ai_disease
  dataset$var2_pooled[which(dataset$any_ai_disease == 1 &
                              (dataset$ref.date > dataset$any_ai_disease_date))] <- 0
  res_leuk_order2_pooled <- clogit(DV ~ var2_pooled + strata(dataset$group), data = dataset)
  conf_results_pooled = rbind(conf_results_pooled, calculate_conf(res_leuk_order2_pooled, paste0(name, "_leuk_prec")))
  
  colnames(conf_results_pooled) <- column_names
}


# -- MULTIVARIATE MODEL --

#Define a function to exctract and format ORs and 95% CIs from a conditional logistic model
calculate_multi <- function(logmodel, name) {
  conf <- round(exp(confint(logmodel)), 2)
  summary_model <- summary(logmodel)
  odds_ratio <- round(exp(summary_model$coefficients[1, 1]), 2)
  temp_result <- c(name, as.character(odds_ratio),
                   as.character(conf[1, 1]), as.character(conf[1, 2]))
  return(temp_result)
}


#Multivariable model
res_leuk_adj_pooled <- clogit(DV ~ any_ai_disease + LGA + TUPAKOINTITUNNUS
                              + strata(leuk_pooled$group), data = leuk_pooled)
summary(res_leuk_adj_pooled) #OR 1.88, 95% CI: 1.40-2.52

#Any autoimmune disease (excluding T1DM) 

#Create empty dataframe to store results
column_names = c("name", "OR", "CI_lower", "CI_upper")
multi_results_pooled = data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(multi_results_pooled) <- column_names

#Fit model across all datasets
for (i in 1:length(datasets_pooled)) {
  dataset = datasets_pooled[[i]]
  name = names_pooled[[i]]
  #Basic model
  res_leuk_pooled <- clogit(DV ~ any_ai_disease + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_pooled = rbind(multi_results_pooled, calculate_multi(res_leuk_pooled, name))
  
  #Model assuming AID precedes leukemia
  dataset$var_pooled <- dataset$any_ai_disease
  dataset$var_pooled[which(dataset$any_ai_disease == 1 &
                             (dataset$ref.date < dataset$any_ai_disease_date))] <- 0
  res_leuk_order_pooled <- clogit(DV ~ var_pooled + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_pooled = rbind(multi_results_pooled, calculate_multi(res_leuk_order_pooled, paste0(name, "_AID_prec")))
  
  #Model assuming leukemia precedes AID
  dataset$var2_pooled <- dataset$any_ai_disease
  dataset$var2_pooled[which(dataset$any_ai_disease == 1 &
                              (dataset$ref.date > dataset$any_ai_disease_date))] <- 0
  res_leuk_order2_pooled <- clogit(DV ~ var2_pooled + LGA + TUPAKOINTITUNNUS + strata(dataset$group), data = dataset)
  multi_results_pooled = rbind(multi_results_pooled, calculate_multi(res_leuk_order2_pooled, paste0(name, "_leuk_prec")))
  
  colnames(multi_results_pooled) <- column_names
}

#Create a forest plot to describe results

forest_plot_df <- tibble::tribble(
  ~Group, ~OR, ~lower, ~upper, ~bold,
  "total", 1.84, 1.41, 2.40, TRUE,
  "female", 1.78, 1.23, 2.58, TRUE,
  "male", 1.90, 1.29, 2.79, TRUE,
  "0-0.99y", 16.12, 1.92, 135, TRUE,
  "1-9.99y", 2.01, 1.43, 2.84, TRUE,
  "10-17.99y", 1.35, 0.86, 2.13, FALSE,
  "ALL", 1.91, 1.43, 2.55, TRUE,
  "AML", 1.20, 0.34, 2.56, FALSE
)

#Order groups top-to-bottom as they appear
forest_plot_df <- forest_plot_df %>%
  mutate(Group = factor(Group, levels = rev(Group)))


#Create the plot
plot <- ggplot(forest_plot_df, aes(x = OR, y = Group)) +
  #confidence interval segments
  geom_segment(aes(x = lower, xend = upper, y = Group, yend = Group),
               linewidth = 0.9, color = "#0b3d91") +
  #point estimate
  geom_point(size = 2.8, color = "#0b3d91") +
  #vertical reference line at OR=1
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  #log scale is standard for ratios
  scale_x_log10(
    breaks = c(0.5, 1, 2, 5, 10, 15),
    labels = scales::label_number(accuracy = 0.1) 
  ) +
  labs(
    x = "Odds ratio (95 % CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(margin = margin(t = 8)),
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave("forest_plot.pdf", plot = plot, device = cairo_pdf(), width = 9, height = 6)
