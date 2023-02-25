library(survival) # load libraries
library(survminer)
library(ggplot2)
setwd("/Users/soniaqi/Documents/QBio-490/qbio_490_sonia/analysis_data") # set dir
read.csv("/Users/soniaqi/Documents/QBio-490/qbio_490_sonia/analysis_data/brca_clinical_data.csv") # import data
clinical <- read.csv("/Users/soniaqi/Documents/QBio-490/qbio_490_sonia/analysis_data/brca_clinical_data.csv")
# I added this part because I had error msg saying "object 'clin_query' not found"
clin_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Clinical",
                       file.type = "xml")
clinical_drug <- GDCprepare_clinic(query = clin_query, clinical.info = "drug") # loading drug and rad data
clinical_rad <- GDCprepare_clinic(query = clin_query, clinical.info = "radiation")
sum(is.na(clinical$age_at_initial_pathologic_diagnosis)) # check # of na
# 1. I chose age at initial pathologic diagnosis. The sum of TRUEs is 0, meaning that there is no missing data for this variable
# 2. It is a discrete numeric variable
sum(is.na(clinical_drug$drug_name)) # check # of na
# 3. I chose drug name
# 4. It's also a categorical variable
# 5. 1) The age of patients would have a signifiance correlation to the drug used in treatment
# 2) The age of patients would have a signifiance correlation to the survival in breast cancer
# 3) The drug used in treatment would have a signifiance correlation to the survival in breast cancer
clinic_drug_merge <- merge(clinical, clinical_drug) # merge clinical and clinical_drug
age_drug <- boxplot(clinic_drug_merge$age_at_initial_pathologic_diagnosis~
                      clinic_drug_merge$drug_name,xlab="drug name",ylab="age") # make boxplot
ggsave("/Users/soniaqi/Documents/QBio-490/qbio_490_sonia/analysis_data/age_drug.png")
# 1. I chose using boxplot because it works for one categorical and one numerical value, 
# while scatterplots require two numeric variables. 
# In the graph, there is no obvious association between the two variables

# KM plot for age
age_na_mask <- ifelse(is.na(clinical$age_at_initial_pathologic_diagnosis),F,T) # mask for na
age_cleaned <- clinical[age_na_mask,] # create new df with no na
young_mask <- ifelse(age_cleaned$age_at_initial_pathologic_diagnosis<=35,T,F) # convert age to catogorical data
middle_mask <- ifelse(age_cleaned$age_at_initial_pathologic_diagnosis>35&age_cleaned$age_at_initial_pathologic_diagnosis<50,T,F)
old_mask <- ifelse(age_cleaned$age_at_initial_pathologic_diagnosis>=50,T,F)
age_cleaned$age_status <- ifelse(young_mask,"Young",ifelse(middle_mask,"Middle","Old")) # create new column to store corresponding categories
age_cleaned$survival_time <- ifelse(is.na(age_cleaned$days_to_death), # create column for survival time
                                    age_cleaned$survival_time <- age_cleaned$days_to_last_followup,
                                    age_cleaned$survival_time <- age_cleaned$days_to_death)
inf_mask <- ifelse(age_cleaned$survival_time=="-Inf",F,T) # mask for -Inf values in survival time
age_cleaned <- age_cleaned[inf_mask,] # remove data with -Inf
age_cleaned$death_event <- ifelse(age_cleaned$vital_status=="Alive", # create bool column to record if the patient is alive
                                  age_cleaned$death_event <- F,
                                  age_cleaned$death_event <- T)
surv_object_age <- Surv(time = age_cleaned$survival_time, event = age_cleaned$death_event) # initialize survival objects
age_fit <- surv_fit(surv_object_age ~ age_cleaned$age_status,data=age_cleaned) # create fit objects
survplot_age <- ggsurvplot(age_fit,pval=T) # create KM plot
KM_plot_age <- survplot_age$plot # save plot to variable
KM_plot_age
ggsave("/Users/soniaqi/Documents/QBio-490/qbio_490_sonia/analysis_data/KM_age.png")

# KM plot for drug name
drug_na_mask <- ifelse(is.na(clinic_drug_merge$drug_name),F,T) # mask for na
drug_cleaned <- clinic_drug_merge[drug_na_mask,] # create new df with no na
drug_cleaned$survival_time <- ifelse(is.na(drug_cleaned$days_to_death), # create column for survival time
                                    drug_cleaned$survival_time <- drug_cleaned$days_to_last_followup,
                                    drug_cleaned$survival_time <- drug_cleaned$days_to_death)
inf_mask <- ifelse(drug_cleaned$survival_time=="-Inf",F,T) # mask for -Inf values in survival time
drug_cleaned <- drug_cleaned[inf_mask,] # remove data with -Inf
drug_cleaned$death_event <- ifelse(drug_cleaned$vital_status=="Alive", # create bool column to record if the patient is alive
                                  drug_cleaned$death_event <- F,
                                  drug_cleaned$death_event <- T)
surv_object_drug <- Surv(time = drug_cleaned$survival_time, event = drug_cleaned$death_event) # initialize survival objects
drug_fit <- surv_fit(surv_object_drug ~ drug_cleaned$drug_name,data=drug_cleaned) # create fit objects
survplot_drug <- ggsurvplot(drug_fit,pval=T) # create KM plot
KM_plot_drug <- survplot_drug$plot # save plot to variable
KM_plot_drug
ggsave("/Users/soniaqi/Documents/QBio-490/qbio_490_sonia/analysis_data/KM_drug.png")
# 4. In the age KM plot, it seems like older people have lower survival rate. 
# However, the p value of 0.38 is greater than 0.05, indicating that the difference is not statistically significant
# In the drug KM plot, the relationship between drug used and survival rate is hard to tell 
# because there are too many different drugs.
