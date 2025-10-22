library(tidyverse)
library(magrittr)
library(TCGAbiolinks)
library(pec)
library(survival)

##########
############### GBM dataset
##########
gbm_data_all <- readRDS("data/gbm_gene_raw.rds")
gbm_geneExp <- SummarizedExperiment::assay(gbm_data_all, 5)
gbm_geneInfo <- SummarizedExperiment::elementMetadata(gbm_data_all)
gbm_geneInfo <- gbm_geneInfo[order(match(gbm_geneInfo$gene_id,
                                         rownames(gbm_geneExp))) , ]
rownames(gbm_geneExp) <- apply(cbind(gbm_geneInfo$gene_name, rownames(gbm_geneExp)),
                               1, function(x) paste(x[1], x[2], sep = "_"))
gbm_geneExp_df <- as.data.frame(t(gbm_geneExp))
gbm_data <-
  SummarizedExperiment::colData(gbm_data_all) %>%
  data.frame %>%
  select(barcode, days_to_death, vital_status, "paper_Survival..months.") %>%
  filter(vital_status %in% c("Alive", "Dead")) %>%
  mutate(vital_status = ifelse(vital_status == "Alive",
                               0, 1),
         tte = ifelse(vital_status == 0,
                      paper_Survival..months.,
                      days_to_death / 30)) %>%
  filter(! is.na(tte),
         tte > 0) %>%
  select(barcode, vital_status, tte)
gbm_geneExp_df <- gbm_geneExp_df[row.names(gbm_geneExp_df) %in%
                                   gbm_data$barcode, ]
gbm_geneExp_df <- gbm_geneExp_df[order(match(row.names(gbm_geneExp_df),
                                             gbm_data$barcode)) , ]
genes_withExp <- apply(gbm_geneExp_df, 2, function(x) (mean(x == 0) < 0.99))
gbm_geneExp_df <- gbm_geneExp_df[ , genes_withExp]
sapply(colnames(gbm_geneExp_df),
       function(x) {
         gbm_data %>%
           select(! barcode) %>%
           bind_cols(gene_tmp = gbm_geneExp_df[, x]) -> gene_df_tmp
        colnames(gene_df_tmp)[colnames(gene_df_tmp) == "gene_tmp"] <- x
        coxph(Surv(tte, vital_status) ~ ., data = gene_df_tmp,
              x = TRUE, y = TRUE) -> cox_tmp
        return(c(AIC = AIC(cox_tmp),
                 IBS = suppressMessages(ibs(pec(list(Cox = cox_tmp), Hist(tte, vital_status) ~ 1, data = gene_df_tmp))[2, 1])))
      },
      simplify = "array") %>%
  t -> gene_gbm_scores
gene_gbm_scores %>%
  data.frame %>%
  arrange(IBS) %>%
  head(n = 50) -> genes_gbm_IBS
gene_gbm_scores %>%
  data.frame %>%
  arrange(AIC) %>%
  head(n = 50) -> genes_gbm_AIC
genes_gbm_merged <-
  genes_gbm_IBS %>%
  bind_rows(genes_gbm_AIC[! rownames(genes_gbm_AIC) %in%
                            rownames(genes_gbm_IBS), ])
gbm_geneExp_df$barcode <- rownames(gbm_geneExp_df)
gbm_data %<>%
  left_join(gbm_geneExp_df[ , c("barcode", rownames(genes_gbm_merged))], by = "barcode") %>%
  mutate_at(rownames(genes_gbm_merged), ~ (scale(.) %>% as.vector))

saveRDS(gbm_data, "data/gbm_fnl.rds")


##########
############### BRCA dataset
##########
brca_data_all <- readRDS("data/brca_gene_raw.rds")
brca_geneExp <- SummarizedExperiment::assay(brca_data_all, 5)
brca_geneInfo <- SummarizedExperiment::elementMetadata(brca_data_all)
brca_data_clinic <- readRDS("data/brca_clinical.rds")
brca_geneInfo <- brca_geneInfo[order(match(brca_geneInfo$gene_id,
                                         rownames(brca_geneExp))) , ]
rownames(brca_geneExp) <- apply(cbind(brca_geneInfo$gene_name, rownames(brca_geneExp)),
                               1, function(x) paste(x[1], x[2], sep = "_"))
brca_geneExp_df <- as.data.frame(t(brca_geneExp))
brca_data <-
  SummarizedExperiment::colData(brca_data_all) %>%
  data.frame %>%
  select(barcode, vital_status, days_to_death) %>%
  mutate(patient_barcode = substr(barcode, start = 1, stop = 12)) %>%
  left_join(brca_data_clinic[ , - which(colnames(brca_data_clinic) %in%
                                          c("vital_status", "days_to_death"))],
            by = c("patient_barcode" = "bcr_patient_barcode")) %>%
  filter(vital_status %in% c("Alive", "Dead"),
         days_to_last_follow_up >= 0) %>%
  mutate(vital_status = ifelse(vital_status == "Alive",
                               0, 1),
         tte = ifelse(vital_status == 0,
                      days_to_last_follow_up / 30,
                      days_to_death / 30)) %>%
  filter(! is.na(tte),
         tte > 0) %>%
  select(barcode, vital_status, tte)
brca_geneExp_df <- brca_geneExp_df[row.names(brca_geneExp_df) %in%
                                   brca_data$barcode, ]
brca_geneExp_df <- brca_geneExp_df[order(match(row.names(brca_geneExp_df),
                                             brca_data$barcode)) , ]
genes_withExp <- apply(brca_geneExp_df, 2, function(x) (mean(x == 0) < 0.99))
brca_geneExp_df <- brca_geneExp_df[ , genes_withExp]
sapply(colnames(brca_geneExp_df),
       function(x) {
         brca_data %>%
           select(! barcode) %>%
           bind_cols(gene_tmp = brca_geneExp_df[, x]) -> gene_df_tmp
         colnames(gene_df_tmp)[colnames(gene_df_tmp) == "gene_tmp"] <- x
         coxph(Surv(tte, vital_status) ~ ., data = gene_df_tmp,
               x = TRUE, y = TRUE) -> cox_tmp
         return(c(AIC = AIC(cox_tmp),
                  IBS = suppressMessages(ibs(pec(list(Cox = cox_tmp), Hist(tte, vital_status) ~ 1, data = gene_df_tmp))[2, 1])))
       },
       simplify = "array") %>%
  t -> gene_brca_scores
gene_brca_scores %>%
  data.frame %>%
  arrange(IBS) %>%
  head(n = 50) -> genes_brca_IBS
gene_brca_scores %>%
  data.frame %>%
  arrange(AIC) %>%
  head(n = 50) -> genes_brca_AIC
genes_brca_merged <-
  genes_brca_IBS %>%
  bind_rows(genes_brca_AIC[! rownames(genes_brca_AIC) %in%
                            rownames(genes_brca_IBS), ])
brca_geneExp_df$barcode <- rownames(brca_geneExp_df)
colnames(SummarizedExperiment::colData(brca_data_all))
brca_data %<>%
  left_join(brca_geneExp_df[ , c("barcode", rownames(genes_brca_merged))], by = "barcode") %>%
  mutate_at(rownames(genes_brca_merged), ~ (scale(.) %>% as.vector))

saveRDS(brca_data, "data/brca_fnl.rds")


##########
############### KIRC dataset
##########
kirc_data_all <- readRDS("data/kirc_gene_raw.rds")
kirc_geneExp <- SummarizedExperiment::assay(kirc_data_all, 5)
kirc_geneInfo <- SummarizedExperiment::elementMetadata(kirc_data_all)
kirc_data_clinic <- readRDS("data/kirc_clinical.rds")
kirc_geneInfo <- kirc_geneInfo[order(match(kirc_geneInfo$gene_id,
                                           rownames(kirc_geneExp))) , ]
rownames(kirc_geneExp) <- apply(cbind(kirc_geneInfo$gene_name, rownames(kirc_geneExp)),
                                1, function(x) paste(x[1], x[2], sep = "_"))
kirc_geneExp_df <- as.data.frame(t(kirc_geneExp))
kirc_data <-
  SummarizedExperiment::colData(kirc_data_all) %>%
  data.frame %>%
  select(barcode, vital_status, days_to_death) %>%
  mutate(patient_barcode = substr(barcode, start = 1, stop = 12)) %>%
  left_join(kirc_data_clinic[ , - which(colnames(kirc_data_clinic) %in%
                                          c("vital_status", "days_to_death"))],
            by = c("patient_barcode" = "bcr_patient_barcode")) %>%
  filter(vital_status %in% c("Alive", "Dead"),
         days_to_last_follow_up >= 0) %>%
  mutate(vital_status = ifelse(vital_status == "Alive",
                               0, 1),
         tte = ifelse(vital_status == 0,
                      days_to_last_follow_up / 30,
                      days_to_death / 30)) %>%
  filter(! is.na(tte),
         tte > 0) %>%
  select(barcode, vital_status, tte)
kirc_geneExp_df <- kirc_geneExp_df[row.names(kirc_geneExp_df) %in%
                                     kirc_data$barcode, ]
kirc_geneExp_df <- kirc_geneExp_df[order(match(row.names(kirc_geneExp_df),
                                               kirc_data$barcode)) , ]
genes_withExp <- apply(kirc_geneExp_df, 2, function(x) (mean(x == 0) < 0.99))
kirc_geneExp_df <- kirc_geneExp_df[ , genes_withExp]
sapply(colnames(kirc_geneExp_df),
       function(x) {
         kirc_data %>%
           select(! barcode) %>%
           bind_cols(gene_tmp = kirc_geneExp_df[, x]) -> gene_df_tmp
         colnames(gene_df_tmp)[colnames(gene_df_tmp) == "gene_tmp"] <- x
         coxph(Surv(tte, vital_status) ~ ., data = gene_df_tmp,
               x = TRUE, y = TRUE) -> cox_tmp
         return(c(AIC = AIC(cox_tmp),
                  IBS = suppressMessages(ibs(pec(list(Cox = cox_tmp), Hist(tte, vital_status) ~ 1, data = gene_df_tmp))[2, 1])))
       },
       simplify = "array") %>%
  t -> gene_kirc_scores
gene_kirc_scores %>%
  data.frame %>%
  arrange(IBS) %>%
  head(n = 50) -> genes_kirc_IBS
gene_kirc_scores %>%
  data.frame %>%
  arrange(AIC) %>%
  head(n = 50) -> genes_kirc_AIC
genes_kirc_merged <-
  genes_kirc_IBS %>%
  bind_rows(genes_kirc_AIC[! rownames(genes_kirc_AIC) %in%
                             rownames(genes_kirc_IBS), ])
kirc_geneExp_df$barcode <- rownames(kirc_geneExp_df)
kirc_data %<>%
  left_join(kirc_geneExp_df[ , c("barcode", rownames(genes_kirc_merged))], by = "barcode") %>%
  mutate_at(rownames(genes_kirc_merged), ~ (scale(.) %>% as.vector))

saveRDS(kirc_data, "data/kirc_fnl.rds")

