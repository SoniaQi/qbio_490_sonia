library(survival) # load libraries
library(ggplot2)
library(survminer)
library(maftools)
library(TCGAbiolinks)

setwd("/Users/soniaqi/Documents/QBio-490/qbio_490_sonia/analysis_data") # set dir
clinical <- read.csv("/Users/soniaqi/Documents/QBio-490/qbio_490_sonia/analysis_data/brca_clinical_data.csv") # load data
maf_query <- GDCquery(project = 'TCGA-BRCA', data.category = "Simple Nucleotide Variation",
  access = "open", data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking") 
#GDCdownload(maf_query)
maf <- GDCprepare(maf_query) # create maf object
maf_object <- read.maf(maf = maf,
                       clinicalData = clinical,
                       isTCGA = TRUE)
# 1. 
table(maf_object@clinical.data$race_list) # select most common categories
black_mask <- ifelse(maf_object@clinical.data$race_list=='BLACK OR AFRICAN AMERICAN',T,F) # mask for black patients
black_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[black_mask]
black_maf <- subsetMaf(maf = maf_object, tsb = black_patient_barcodes)
white_mask <- ifelse(maf_object@clinical.data$race_list=='WHITE',T,F) # mask for white patients
white_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[white_mask]
white_maf <- subsetMaf(maf = maf_object, tsb = white_patient_barcodes)
maf_object@clinical.data$race_mod <- ifelse(maf_object@clinical.data$race_list=='BLACK OR AFRICAN AMERICAN', # create new column
                                            'Black', ifelse(maf_object@clinical.data$race_list=='WHITE',
                                                            'White', 'NA'))
# 2.
black.genes = getGeneSummary(black_maf)[10:20] # select most common genes
white.genes = getGeneSummary(white_maf)[10:20]
mdt = merge(black.genes[,.(Hugo_Symbol, MutatedSamples)], white.genes[,.(Hugo_Symbol, MutatedSamples)], by = 'Hugo_Symbol', all = TRUE)
mdt$MutatedSamples.x[is.na(mdt$MutatedSamples.x)] = 0
mdt$MutatedSamples.y[is.na(mdt$MutatedSamples.y)] = 0
mdt$max = apply(mdt[,.(MutatedSamples.x, MutatedSamples.y)], 1, max)
mdt = mdt[order(max, decreasing = TRUE)]
par(mar=c(1,1,1,1)) # set margins
coOncoplot(m1 = black_maf, m2 = white_maf, genes = mdt$Hugo_Symbol, # create cooncoplot
           m1Name = 'Black Patients', m2Name = 'White Patients', borderCol = NA)
# The FLG gene is used to make protein profilaggrin, which make up the outermost layer of skin
# and helps the hydration of the skin. The mutation rate in this gene is different for white
# and black/African American because they are originated from different areas that have
# different sunlight intensity.
ggsave("/Users/soniaqi/Documents/QBio-490/qbio_490_sonia/week7_maf/genderCooncoplot.png")
# 3. 
FLG_maf <- subsetMaf(maf = maf_object,genes = 'FLG') # mask for FLG
mut_pats_FLG <- FLG_maf@clinical.data$Tumor_Sample_Barcode
num_pats_FLG <- length(mut_pats_FLG) # num of patients with FLG mutations
mut_black_pats <- intersect(mut_pats_FLG, black_patient_barcodes) # find num of black patients with FLG mutations
num_black_FLG <- length(mut_black_pats)
mut_white_pats <- intersect(mut_pats_FLG, white_patient_barcodes) # find num of white patients with FLG mutations
num_white_FLG <- length(mut_white_pats)
num_black_noFLG <- length(black_patient_barcodes) - num_black_FLG # find num of black patients w/out FLG mutations
num_white_noFLG <- length(white_patient_barcodes) - num_white_FLG # find num of white patients w/out FLG mutations
contig <- matrix(c(num_black_FLG, num_black_noFLG,num_white_FLG,num_white_noFLG), nrow=2) # create contig matrix
fisher_test <- fisher.test(contig) # perform fisher's test
fisher_test
# The chance of getting a mutated FLG if the patient is black or African American is around
# 0.358. Because the p-value of 0.0234 is less than 0.05, we have statistically significant
# evidence to reject the null hyothesis. 
mosaicplot(contig) # create mosiacplot
ggsave("/Users/soniaqi/Documents/QBio-490/qbio_490_sonia/week7_maf/RaceFLGMosaicplot.png")
# 4. 
lollipopPlot2(m1 = black_maf, m2 = white_maf, m1_name = 'Black Patients', m2_name = 'White Patients', gene = "FLG") # create lollipop plot
ggsave("/Users/soniaqi/Documents/QBio-490/qbio_490_sonia/week7_maf/RaceFLGLollipopplot.png")
# White patients have a lot more missense mutation in FLG than Black patients. Also, frame 
# shift insertion and nonsense mutation only occurred in white patients.
# 5. 
maf_object@clinical.data$FLG_status <- ifelse(maf_object@clinical.data$Tumor_Sample_Barcode %in% mut_pats_FLG,
                                              "mutated","not mutated") # mask for FLG
maf_object@clinical.data$survival_time <- ifelse(is.na(maf_object@clinical.data$days_to_death), # create column for survival time
                                                 maf_object@clinical.data$survival_time <- maf_object@clinical.data$days_to_last_followup,
                                                 maf_object@clinical.data$survival_time <- maf_object@clinical.data$days_to_death)
inf_mask <- ifelse(maf_object@clinical.data$survival_time=="-Inf",F,T) # mask for -Inf values in survival time
FLG_cleaned <- maf_object@clinical.data[inf_mask,] # remove data with -Inf
na_maf <- ifelse(is.na(FLG_cleaned$vital_status),F,T) # mask for NA vital status
FlG_cleaned <- FLG_cleaned[na_maf,]
FLG_cleaned$death_event <- ifelse(FLG_cleaned$vital_status=="Alive", # create bool column to record if the patient is alive
                                   FLG_cleaned$death_event <- F,
                                   FLG_cleaned$death_event <- T)
FLG_cleaned$survival_time <- as.numeric(FLG_cleaned$survival_time)
surv_object_FLG <- Surv(time = FLG_cleaned$survival_time, event = FLG_cleaned$death_event) # initialize survival objects
FLG_fit <- surv_fit(surv_object_FLG ~ FLG_cleaned$FLG_status,data=FLG_cleaned) # create fit objects
survplot_FLG <- ggsurvplot(FLG_fit,pval=T) # create KM plot
KM_plot_FLG <- survplot_FLG$plot # save plot to variable
KM_plot_FLG
ggsave("/Users/soniaqi/Documents/QBio-490/qbio_490_sonia/week7_maf/FLGKMplot.png")
# In the FLG KM plot, it seems like people with mutated FLG have lower survival 
# rate during 2000 to 3000 days and higher survival rate after 4000 days. 
# However, the p value of 0.89 is greater than 0.05, indicating that the difference 
# is not statistically significant.
# Study shows that loss-of-function mutation in FLG is associated with a higher
# serum vitamin D level, which may protect against cancer through enhancing skin 
# barrier. This somewhat supports the higher survival rate of patients with mutated
# FLG. 