# BiocManager::install("TCGAbiolinks")

# Step 2: Load the library
library(TCGAbiolinks)
library(SummarizedExperiment) # To access assay(), colData(), etc.
library(tidyverse)


gbm_barcode_list <- readRDS("data/gbm_barcode_list.rds")
brca_barcode_list <- readRDS("data/brca_barcode_list.rds")
kirc_barcode_list <- readRDS("data/kirc_barcode_list.rds")

##########
############### GBM dataset
##########
## Gene expression
gbm_query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = gbm_barcode_list
)
GDCdownload(gbm_query)
gbm_data <- GDCprepare(gbm_query)
## clinical
gbm_clinical_query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  barcode = substr(gbm_barcode_list, 1, 12)
)
GDCdownload(gbm_clinical_query)
gbm_clinical_data <- GDCprepare(gbm_clinical_query)
gbm_clinical <- GDCquery_clinic(
  project = "TCGA-GBM",
  type = "Clinical"
)

saveRDS(gbm_clinical, "data/gbm_clinical.rds")
saveRDS(gbm_clinical_data, "data/gbm_clinical_raw.rds")
saveRDS(gbm_data, "data/gbm_gene_raw.rds")

rm(gbm_clinical_query, gbm_query, gbm_barcode_list,
   gbm_data, gbm_clinical_data)
gc()

##########
############### BRCA dataset
##########
## Gene expression
brca_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = brca_barcode_list
)
GDCdownload(brca_query)
brca_data <- GDCprepare(brca_query)
## clinical
brca_clinical_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  barcode = substr(brca_barcode_list, 1, 12)
)
GDCdownload(brca_clinical_query, files.per.chunk = 100)
brca_clinical_data <- GDCprepare(brca_clinical_query)
brca_clinical <- GDCquery_clinic(
  project = "TCGA-BRCA",
  type = "Clinical"
)

saveRDS(brca_clinical, "data/brca_clinical.rds")
saveRDS(brca_clinical_data, "data/brca_clinical_raw.rds")
saveRDS(brca_data, "data/brca_gene_raw.rds")

rm(brca_clinical, brca_clinical_query, brca_query, brca_barcode_list,
   brca_data, brca_clinical_data)
gc()

##########
############### KIRC dataset
##########
## Gene expression
kirc_query <- GDCquery(
  project = "TCGA-KIRC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = kirc_barcode_list
)
GDCdownload(kirc_query, files.per.chunk = 70)
kirc_data <- GDCprepare(kirc_query)
## clinical
kirc_clinical_query <- GDCquery(
  project = "TCGA-KIRC",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  barcode = substr(kirc_barcode_list, 1, 12)
)
GDCdownload(kirc_clinical_query)
kirc_clinical_data <- GDCprepare(kirc_clinical_query)
kirc_clinical <- GDCquery_clinic(
  project = "TCGA-KIRC",
  type = "Clinical"
)

saveRDS(kirc_clinical, "data/kirc_clinical.rds")
saveRDS(kirc_clinical_data, "data/kirc_clinical_raw.rds")
saveRDS(kirc_data, "data/kirc_gene_raw.rds")

rm(kirc_clinical_query, kirc_query, kirc_barcode_list,
   kirc_data, kirc_clinical_data)
gc()



