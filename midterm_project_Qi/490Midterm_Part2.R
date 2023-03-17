library(survival) 
library(survminer)
library(ggplot2)
library(DESeq2)
library(maftools)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(EnhancedVolcano)

setwd("/Users/soniaqi/Documents/QBio-490/qbio_490_sonia/midterm_project_Qi")
#dir.create("outputs")

clin_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Clinical",
                       file.type = "xml")
#GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query,
                            clinical.info = "patient")
sum(is.na(clinic$age_at_initial_pathologic_diagnosis)) # no na value, no mask needed
clinic$age_category <- ifelse(clinic$age_at_initial_pathologic_diagnosis < 45, "Young", "Old")

clinical_drug <- GDCprepare_clinic(query = clin_query, clinical.info = "drug") # loading drug and rad data
clinical_rad <- GDCprepare_clinic(query = clin_query, clinical.info = "radiation")

colnames(clinic)[colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", # we only have access to somatic mutations which are open access
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

#GDCdownload(maf_query)

maf <- GDCprepare(maf_query) # as long as it runs, ignore any errors
colnames(clinic)[colnames(clinic) == "bcr_patient_barcode" ] <-
  "Tumor_Sample_Barcode"
maf_object <- read.maf(maf = maf, 
                       clinicalData = clinic,
                       isTCGA = TRUE)

rna_query <- GDCquery(project ="TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

#GDCdownload(rna_query)
rna_se <- GDCprepare(rna_query)

rna_clinical <- rna_se@colData
rna_clinical <- as.data.frame(rna_clinical)

rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)

rna_counts <- rna_se@assays@data$unstranded

rownames(rna_genes) <- rna_genes$gene_id
colnames(rna_counts) <- rna_clinical$barcode
rownames(rna_counts) <- rna_genes$gene_id

# create a maf object with only young patients
young_mask <- ifelse(maf_object@clinical.data$age_category=='Young',T,F) 
young_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[young_mask]
young_maf <- subsetMaf(maf = maf_object,
                       tsb = young_patient_barcodes)
# find the top 5 most common mutations in young patients
oncoplot(maf = young_maf,
         top = 5,
         borderCol = NA) 
ggsave("outputs/young_oncoplot.png")

# create survival status column for young_maf
young_maf@clinical.data$Overall_Survival_Status <- ifelse(young_maf@clinical.data$vital_status=='Alive',T,F)
par(mar=c(1,1,1,1))
# create a survival plot for young patients with TP53 mutations
mafSurvival(maf = young_maf,
            genes = "TP53",
            time = "days_to_last_followup", 
            Status = "Overall_Survival_Status", 
            isTCGA = TRUE)
ggsave("outputs/young_mafsurvplot.png")

# create KM plot for TP53 mutation in young patients
TP53_maf <- subsetMaf(maf = young_maf,genes = 'TP53') # mask for TP53
mut_pats_TP53 <- TP53_maf@clinical.data$Tumor_Sample_Barcode
young_maf@clinical.data$TP53_status <- ifelse(young_maf@clinical.data$Tumor_Sample_Barcode %in% mut_pats_TP53,
                                              "mutated","not mutated")
young_maf@clinical.data$survival_time <- ifelse(is.na(young_maf@clinical.data$days_to_death), # create column for survival time
                                                young_maf@clinical.data$survival_time <- young_maf@clinical.data$days_to_last_followup,
                                                young_maf@clinical.data$survival_time <- young_maf@clinical.data$days_to_death)
inf_mask <- ifelse(young_maf@clinical.data$survival_time=="-Inf",F,T) # mask for -Inf values in survival time
TP53_cleaned <- young_maf@clinical.data[inf_mask,] # remove data with -Inf
na_maf <- ifelse(is.na(TP53_cleaned$vital_status),F,T) # mask for NA vital status
TP53_cleaned <- TP53_cleaned[na_maf,]
TP53_cleaned$death_event <- ifelse(TP53_cleaned$vital_status=="Alive", # create bool column to record if the patient is alive
                                   TP53_cleaned$death_event <- F,
                                   TP53_cleaned$death_event <- T)
TP53_cleaned$survival_time <- as.numeric(TP53_cleaned$survival_time)
surv_object_TP53 <- Surv(time = TP53_cleaned$survival_time, event = TP53_cleaned$death_event) # initialize survival objects
TP53_fit <- surv_fit(surv_object_TP53 ~ TP53_cleaned$TP53_status,data=TP53_cleaned) # create fit objects
survplot_TP53 <- ggsurvplot(TP53_fit,pval=T) # create KM plot
KM_plot_TP53 <- survplot_TP53$plot # save plot to variable
KM_plot_TP53
ggsave("outputs/young_TP53survplot.png")

# create boxplot on the survival time of young patients with/without TP53 mutations
survtime_mask <- !is.na(young_maf@clinical.data$survival_time)
young_maf@clinical.data <- young_maf@clinical.data[survtime_mask,]
young_maf@clinical.data$survival_time <- as.numeric(young_maf@clinical.data$survival_time)
young_maf@clinical.data$TP53_mut <- ifelse(young_maf@clinical.data$Tumor_Sample_Barcode %in% TP53_cleaned$Tumor_Sample_Barcode, T, F)
boxplot(young_maf@clinical.data$survival_time ~ young_maf@clinical.data$TP53_mut,
        xlab = "mutation status", ylab = "survival time", 
        main = "influence of TP53 mutation on survival time")
ggsave("outputs/TP53boxplot.png")

# adjust data for dds
rna_clinical$age_category <- ifelse(rna_clinical$age_at_index<45, "Young", "Old")
rna_clinical$age_category <- factor(rna_clinical$age_category)
rna_clinical$ajcc_pathologic_stage <- factor(rna_clinical$ajcc_pathologic_stage)

sum(is.na(rna_clinical$age_category))
sum(is.na(rna_clinical$ajcc_pathologic_stage))

na_mask <-  ifelse(is.na(rna_clinical$age_category)|is.na(rna_clinical$ajcc_pathologic_stage), F, T)
rna_clinical <- rna_clinical[na_mask,]
rna_counts <- rna_counts[,na_mask]

row_sums <- rowSums(rna_counts)
# create a boolean mask where genes with < 10 total counts are FALSE, and genes with >= 10 total counts are TRUE
low_counts_mask <- ifelse(row_sums<10,F,T)
# rewrite the rna_counts df, subsetting for only genes with >= 10 total counts
rna_counts <- rna_counts[low_counts_mask,]
#update rna_genes with the low_counts_mas
rna_genes <- rna_genes[low_counts_mask,]

# create dds object
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = rna_clinical,
                              design = ~ajcc_pathologic_stage + age_category)

dds_obj <- DESeq(dds) 
resultsNames(dds_obj)  # see what comparisons got run
# get the young vs. old comparison
results <- results(dds_obj, format = "DataFrame", contrast = c("age_category", "Young", "Old")) # this is case sensitive so be careful to match it with your age_category factors closely!

# rename columns
results <- data.frame(rna_genes$gene_name, results@rownames, results$log2FoldChange, results$pvalue, results$padj, -log10(results$padj))
colnames(results) <- c("gene_name", "gene_id", "log2foldchange", "pvalue", "padj", "neglog10padj") ## FIX

genename_mask <- ifelse(rna_genes$gene_id %in% results$gene_id, T, F) 
rna_genes <- rna_genes[genename_mask,]

padj_mask <- ifelse(results$padj<0.05, T, F)
sig_results <- results[padj_mask, ]

up_reg_results <- results[order(results$log2foldchange, decreasing = T),]
up_reg_mask <- ifelse(up_reg_results$log2foldchange>1, T, F)
up_reg_results <- up_reg_results[up_reg_mask, ]

down_reg_results <- results[order(results$log2foldchange),]
down_reg_mask <- ifelse(down_reg_results$log2foldchange < (-1), T, F)
down_reg_results <- down_reg_results[down_reg_mask, ]

rownames(up_reg_results) <- up_reg_results$gene_id
rownames(down_reg_results) <- down_reg_results$gene_id
# create volcano plot
EnhancedVolcano(results,
                lab = results$gene_name,
                x = 'log2foldchange',
                y = 'neglog10padj',
                title = 'Young vs. Old', 
                pCutoff = 1.3)
ggsave("outputs/agevolcanoplot.png")

