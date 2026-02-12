# BiocManager::install("TCGAbiolinks")

# Step 2: Load the library
library(TCGAbiolinks)
library(SummarizedExperiment) # To access assay(), colData(), etc.
library(tidyverse)


##########
############### GBM dataset
##########
## Gene expression
gbm_query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(gbm_query,
            directory = "GDCdata",
            files.per.chunk = 50)
gbm_data <- GDCprepare(gbm_query,
                       directory = "GDCdata")

## clinical
gbm_clinical_query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Clinical",
  data.type = "Clinical Supplement"
)
GDCdownload(gbm_clinical_query,
            directory = "GDCdata",
            files.per.chunk = 200)
gbm_clinical_data <- GDCprepare(gbm_clinical_query,
                                directory = "GDCdata")
gbm_clinical <- GDCquery_clinic(
  project = "TCGA-GBM",
  type = "Clinical"
)

saveRDS(gbm_clinical, "data/raw/gbm_clinical.rds")
saveRDS(gbm_clinical_data, "data/raw/gbm_clinical_raw.rds")
saveRDS(gbm_data, "data/raw/gbm_gene_raw.rds")

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
  workflow.type = "STAR - Counts"
)
GDCdownload(brca_query,
            directory = "GDCdata",
            files.per.chunk = 50)
# mem.maxVSize(vsize = 18000) # in case, max available memory reached
brca_data <- GDCprepare(brca_query,
                        directory = "GDCdata")
# mem.maxVSize(vsize = 16000) # set back to default
## clinical
brca_clinical_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement"
)
GDCdownload(brca_clinical_query,,
            directory = "GDCdata",
            files.per.chunk = 100)
brca_clinical_data <- GDCprepare(brca_clinical_query,
                                 directory = "GDCdata")
brca_clinical <- GDCquery_clinic(
  project = "TCGA-BRCA",
  type = "Clinical"
)

saveRDS(brca_clinical, "data/raw/brca_clinical.rds")
saveRDS(brca_clinical_data, "data/raw/brca_clinical_raw.rds")
saveRDS(brca_data, "data/raw/brca_gene_raw.rds")

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
  workflow.type = "STAR - Counts"
)
GDCdownload(kirc_query,
            directory = "GDCdata",
            files.per.chunk = 50)
kirc_data <- GDCprepare(kirc_query,
                        directory = "GDCdata")
saveRDS(kirc_data, "data/raw/kirc_gene_raw.rds")
rm(kirc_query, kirc_data)
gc()

## clinical
kirc_clinical_query <- GDCquery(
  project = "TCGA-KIRC",
  data.category = "Clinical",
  data.type = "Clinical Supplement"
)
GDCdownload(kirc_clinical_query,
            directory = "GDCdata",
            files.per.chunk = 200)
kirc_clinical_data <- GDCprepare(kirc_clinical_query,
                                 directory = "GDCdata")
kirc_clinical <- GDCquery_clinic(
  project = "TCGA-KIRC",
  type = "Clinical"
)

saveRDS(kirc_clinical, "data/raw/kirc_clinical.rds")
saveRDS(kirc_clinical_data, "data/raw/kirc_clinical_raw.rds")

rm(kirc_clinical_query, kirc_clinical_data, kirc_clinical)
gc()


