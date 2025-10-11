## Required libraries
library(tidyverse)
library(magrittr)
library(survival)
library(glmnet)
library(pec)
# library(devtools)
# install_github(repo="trevorhastie/uniLasso")
library(uniLasso)

## set parameters
n_bins <- 50
sampleSizes <- c(300, 500, 1000, 2000, 4000)

## Load simulated data as well as results for Binacox
cuts_betas_org <- read_delim("simData/Scens1_3/Scen1_3_cutPoints_beta.csv",
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE,
                             show_col_types = FALSE)
true_cuts_betas <-
  cuts_betas_org[-1, ] %>%
  select(sample_size, iteration, cut_points_x1, cut_points_x2, beta_star) %>%
  rowwise %>%
  mutate(X1_cuts = list(as.numeric(strsplit(str_sub(cut_points_x1, start = 2, end = nchar(cut_points_x1) - 1), "\\s+")[[1]])),
         X2_cuts = list(as.numeric(strsplit(str_sub(cut_points_x2, start = 2, end = nchar(cut_points_x2) - 1), "\\s+")[[1]])),
         betas = list(as.numeric(strsplit(str_sub(beta_star, start = 2, end = nchar(beta_star) - 1), "\\s+")[[1]]))) %>%
  select(sample_size, iteration, X1_cuts, X2_cuts, betas) %>%
  ungroup %>%
  filter(sample_size %in% sampleSizes)

Bussy_cuts_betas <-
  cuts_betas_org[-1, ] %>%
  select(sample_size, iteration,
         `Processing_Time (seconds) sc1`, `Processing_Time (seconds) sc3`,
         cut_points_estimate_x1_sc1,	cut_points_estimate_x2_sc1,
         cut_points_estimate_x1_sc3,	cut_points_estimate_x2_sc3) %>%
  rowwise %>%
  mutate(X1_cuts_sc1 = list(as.numeric(strsplit(str_sub(cut_points_estimate_x1_sc1, start = 2, end = nchar(cut_points_estimate_x1_sc1) - 1), "\\s+")[[1]])),
         X2_cuts_sc1 = list(as.numeric(strsplit(str_sub(cut_points_estimate_x2_sc1, start = 2, end = nchar(cut_points_estimate_x2_sc1) - 1), "\\s+")[[1]])),
         
         X1_cuts_sc3 = list(as.numeric(strsplit(str_sub(cut_points_estimate_x1_sc3, start = 2, end = nchar(cut_points_estimate_x1_sc3) - 1), "\\s+")[[1]])),
         X2_cuts_sc3 = list(as.numeric(strsplit(str_sub(cut_points_estimate_x2_sc3, start = 2, end = nchar(cut_points_estimate_x2_sc3) - 1), "\\s+")[[1]]))) %>%
  rename(c_time_sc1 = `Processing_Time (seconds) sc1`,
         c_time_sc3 = `Processing_Time (seconds) sc3`) %>%
  select(sample_size, iteration,
         c_time_sc1, c_time_sc3,
         X1_cuts_sc1, X2_cuts_sc1,
         X1_cuts_sc3, X2_cuts_sc3) %>%
  ungroup %>%
  filter(sample_size %in% sampleSizes)

file_names <- list.files("simData/Scens1_3")
file_names <- file_names[- which(grepl("Scen1_3_cutPoints", file_names))]
file_names_filter_tmp <-
  sapply(file_names,
         function(x){
           ss_tmp <- as.numeric(strsplit(x, split = "_")[[1]][2])
           if (! ss_tmp %in% sampleSizes) ss_tmp = NA
           return(ss_tmp)
         })
file_names_filter <- names(file_names_filter_tmp[! is.na(file_names_filter_tmp)])

## Fit models using biniLasso & miniLasso
res_df <- tibble(SS = rep(NA, length(file_names_filter)),
                 iter = rep(NA, length(file_names_filter)),
                 c_time_sc1 = rep(NA, length(file_names_filter)),
                 AIC_sc1 = rep(NA, length(file_names_filter)),
                 BIC_sc1 = rep(NA, length(file_names_filter)),
                 IBS_sc1 = rep(NA, length(file_names_filter)),
                 M1_sc1 = rep(NA, length(file_names_filter)),
                 c_time_sc1_uni = rep(NA, length(file_names_filter)),
                 AIC_sc1_uni = rep(NA, length(file_names_filter)),
                 BIC_sc1_uni = rep(NA, length(file_names_filter)),
                 IBS_sc1_uni = rep(NA, length(file_names_filter)),
                 M1_sc1_uni = rep(NA, length(file_names_filter)),
                 c_time_sc3 = rep(NA, length(file_names_filter)),
                 AIC_sc3 = rep(NA, length(file_names_filter)),
                 BIC_sc3 = rep(NA, length(file_names_filter)),
                 IBS_sc3 = rep(NA, length(file_names_filter)),
                 M1_sc3 = rep(NA, length(file_names_filter)),
                 c_time_sc3_uni = rep(NA, length(file_names_filter)),
                 AIC_sc3_uni = rep(NA, length(file_names_filter)),
                 BIC_sc3_uni = rep(NA, length(file_names_filter)),
                 IBS_sc3_uni = rep(NA, length(file_names_filter)),
                 M1_sc3_uni = rep(NA, length(file_names_filter))) %>%
  mutate(x1_cuts_sc1 = list(list(NULL)),
         x2_cuts_sc1 = list(list(NULL)),
         betas_sc1 = list(list(NULL)),
         x1_cuts_sc1_uni = list(list(NULL)),
         x2_cuts_sc1_uni = list(list(NULL)),
         betas_sc1_uni = list(list(NULL)),
         x1_cuts_sc3 = list(list(NULL)),
         x2_cuts_sc3 = list(list(NULL)),
         betas_sc3 = list(list(NULL)),
         x1_cuts_sc3_uni = list(list(NULL)),
         x2_cuts_sc3_uni = list(list(NULL)),
         betas_sc3_uni = list(list(NULL)))
Bussy_cuts_betas %<>%
  mutate(AIC_sc1 = NA,
         BIC_sc1 = NA,
         IBS_sc1 = NA,
         M1_sc1 = NA,
         betas_sc1 = list(list(NULL)),
         AIC_sc3 = NA,
         BIC_sc3 = NA,
         IBS_sc3 = NA,
         M1_sc3 = NA,
         betas_sc3 = list(list(NULL)))
true_cuts_betas$AIC_sc1 <- true_cuts_betas$BIC_sc1 <- true_cuts_betas$IBS_sc1 <-
  true_cuts_betas$AIC_sc3 <- true_cuts_betas$BIC_sc3 <- true_cuts_betas$IBS_sc3 <- NA
true_cuts_betas$betas_est_sc1 <- true_cuts_betas$betas_est_sc3 <- list(list(NULL))

for (i in 1 : nrow(true_cuts_betas)) {
  file_tmp <- paste0("simData/Scens1_4/", file_names_filter[i])
  simData <- read_delim(file_tmp,
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE, show_col_types = FALSE)
  simData_sc1 <-
    simData %>%
    select(X1, X2, Y, delta)
  simData_sc3 <-
    simData %>%
    select(X1, X2, X1_new, X2_new, Y_cont, delta_cont) %>%
    rename(Y = Y_cont, delta = delta_cont) %>%
    mutate(X1 = X1 * 3,
           X2 = X2 * 3)
  rm(simData)
  
  vals_tmp <- strsplit(file_tmp, "_")[[1]]
  res_df$SS[i] <- as.numeric(vals_tmp[2])
  res_df$iter[i] <- as.numeric(str_sub(vals_tmp[3], start = 1, end = nchar(vals_tmp[3]) - 4))
  
  x1_bounds_sc1 <- quantile(simData_sc1$X1, probs = c(1 : n_bins) / (1 + n_bins))
  x2_bounds_sc1 <- quantile(simData_sc1$X2, probs = c(1 : n_bins) / (1 + n_bins))
  x1_bounds_sc3 <- quantile(simData_sc3$X1, probs = c(1 : n_bins) / (1 + n_bins))
  x2_bounds_sc3 <- quantile(simData_sc3$X2, probs = c(1 : n_bins) / (1 + n_bins))
  
  for (j in 1 : n_bins) {
    simData_sc1 %<>%
      mutate(temp1 = ifelse(X1 >= x1_bounds_sc1[j], 1, 0),
             temp2 = ifelse(X2 >= x2_bounds_sc1[j], 1, 0))
    simData_sc3 %<>%
      mutate(temp1 = ifelse(X1 >= x1_bounds_sc3[j], 1, 0),
             temp2 = ifelse(X2 >= x2_bounds_sc3[j], 1, 0))
    if (j < 10){
      colnames(simData_sc1)[colnames(simData_sc1) == "temp1"] <- paste0("X1_bin_0", j)
      colnames(simData_sc1)[colnames(simData_sc1) == "temp2"] <- paste0("X2_bin_0", j)
      colnames(simData_sc3)[colnames(simData_sc3) == "temp1"] <- paste0("X1_bin_0", j)
      colnames(simData_sc3)[colnames(simData_sc3) == "temp2"] <- paste0("X2_bin_0", j)
    }
    else {
      colnames(simData_sc1)[colnames(simData_sc1) == "temp1"] <- paste0("X1_bin_", j)
      colnames(simData_sc1)[colnames(simData_sc1) == "temp2"] <- paste0("X2_bin_", j)
      colnames(simData_sc3)[colnames(simData_sc3) == "temp1"] <- paste0("X1_bin_", j)
      colnames(simData_sc3)[colnames(simData_sc3) == "temp2"] <- paste0("X2_bin_", j)
    }
  }
  
  x_min_max_sc1 <- list(c(min(simData_sc1$X1),
                          max(simData_sc1$X1)),
                        c(min(simData_sc1$X2),
                          max(simData_sc1$X2)))
  x_min_max_sc3 <- list(c(min(simData_sc3$X1),
                          max(simData_sc3$X1)),
                        c(min(simData_sc3$X2),
                          max(simData_sc3$X2)))
  
  simData_mod_sc1 <-
    simData_sc1 %>%
    dplyr::select(Y, delta) %>%
    bind_cols(simData_sc1[ , colnames(simData_sc1)[grepl("X1_bin_", colnames(simData_sc1))]]) %>%
    bind_cols(simData_sc1[ , colnames(simData_sc1)[grepl("X2_bin_", colnames(simData_sc1))]])
  X_mat_sc1 <-
    model.matrix( ~ . -1,
                  data = simData_mod_sc1[ , colnames(simData_mod_sc1)[grepl("_bin_", colnames(simData_mod_sc1))]])
  simData_mod_sc3 <-
    simData_sc3 %>%
    dplyr::select(Y, delta) %>%
    bind_cols(simData_sc3[ , colnames(simData_sc3)[grepl("X1_bin_", colnames(simData_sc3))]]) %>%
    bind_cols(simData_sc3[ , colnames(simData_sc3)[grepl("X2_bin_", colnames(simData_sc3))]])
  X_mat_sc3 <-
    model.matrix( ~ . -1,
                  data = simData_mod_sc3[ , colnames(simData_mod_sc3)[grepl("_bin_", colnames(simData_mod_sc3))]])
  
  t1 <- Sys.time()
  glm_cv_sc1 <- cv.glmnet(x = X_mat_sc1, y = Surv(simData_mod_sc1$Y,
                                                  simData_mod_sc1$delta),
                          family = "cox", nfolds = 5, relax = FALSE)
  glm_cox_sc1 <- glmnet(x = X_mat_sc1, y = Surv(simData_mod_sc1$Y, simData_mod_sc1$delta),
                        family = "cox", lambda = glm_cv_sc1$lambda.1se)
  t2 <- Sys.time()
  res_df$c_time_sc1[i] <- as.numeric(t2 - t1)
  beta_nonZero <- names(glm_cox_sc1$beta[ , 1][glm_cox_sc1$beta[, 1] != 0])
  X1_cuts_ind <- as.numeric(substr(beta_nonZero[grepl("X1", beta_nonZero)], 8, 9))
  X2_cuts_ind <- as.numeric(substr(beta_nonZero[grepl("X2", beta_nonZero)], 8, 9))
  X1_cuts_ind <- X1_cuts_ind[c(1, which(diff(X1_cuts_ind) > 1) + 1)]
  X2_cuts_ind <- X2_cuts_ind[c(1, which(diff(X2_cuts_ind) > 1) + 1)]
  res_df$x1_cuts_sc1[i] <- list(unique(x1_bounds_sc1[X1_cuts_ind]))
  res_df$x2_cuts_sc1[i] <- list(unique(x2_bounds_sc1[X2_cuts_ind]))
  
  t1 <- Sys.time()
  unilasso_cv_sc1 <- cv.uniLasso(x = X_mat_sc1, y = Surv(simData_mod_sc1$Y,
                                                         simData_mod_sc1$delta),
                                 family = "cox", nfolds = 5, relax = FALSE)
  unilasso_cox_sc1 <- uniLasso(x = X_mat_sc1, y = Surv(simData_mod_sc1$Y, simData_mod_sc1$delta),
                               family = "cox", lambda = unilasso_cv_sc1$lambda.1se)
  t2 <- Sys.time()
  res_df$c_time_sc1_uni[i] <- as.numeric(t2 - t1)
  beta_nonZero <- names(unilasso_cox_sc1$beta[ , 1][unilasso_cox_sc1$beta[, 1] != 0])
  X1_cuts_ind <- as.numeric(substr(beta_nonZero[grepl("X1", beta_nonZero)], 8, 9))
  X2_cuts_ind <- as.numeric(substr(beta_nonZero[grepl("X2", beta_nonZero)], 8, 9))
  X1_cuts_ind <- X1_cuts_ind[c(1, which(diff(X1_cuts_ind) > 1) + 1)]
  X2_cuts_ind <- X2_cuts_ind[c(1, which(diff(X2_cuts_ind) > 1) + 1)]
  res_df$x1_cuts_sc1_uni[i] <- list(unique(x1_bounds_sc1[X1_cuts_ind]))
  res_df$x2_cuts_sc1_uni[i] <- list(unique(x2_bounds_sc1[X2_cuts_ind]))
  
  t1 <- Sys.time()
  glm_cv_sc3 <- cv.glmnet(x = X_mat_sc3, y = Surv(simData_mod_sc3$Y,
                                                  simData_mod_sc3$delta),
                          family = "cox", nfolds = 5, relax = FALSE)
  glm_cox_sc3 <- glmnet(x = X_mat_sc3, y = Surv(simData_mod_sc3$Y, simData_mod_sc3$delta),
                        family = "cox", lambda = glm_cv_sc3$lambda.1se)
  t2 <- Sys.time()
  res_df$c_time_sc3[i] <- as.numeric(t2 - t1)
  beta_nonZero <- names(glm_cox_sc3$beta[ , 1][glm_cox_sc3$beta[, 1] != 0])
  X1_cuts_ind <- as.numeric(substr(beta_nonZero[grepl("X1", beta_nonZero)], 8, 9))
  X2_cuts_ind <- as.numeric(substr(beta_nonZero[grepl("X2", beta_nonZero)], 8, 9))
  X1_cuts_ind <- X1_cuts_ind[c(1, which(diff(X1_cuts_ind) > 1) + 1)]
  X2_cuts_ind <- X2_cuts_ind[c(1, which(diff(X2_cuts_ind) > 1) + 1)]
  res_df$x1_cuts_sc3[i] <- list(unique(x1_bounds_sc3[X1_cuts_ind]))
  res_df$x2_cuts_sc3[i] <- list(unique(x2_bounds_sc3[X2_cuts_ind]))
  
  t1 <- Sys.time()
  unilasso_cv_sc3 <- cv.uniLasso(x = X_mat_sc3, y = Surv(simData_mod_sc3$Y,
                                                         simData_mod_sc3$delta),
                                 family = "cox", nfolds = 5, relax = FALSE)
  unilasso_cox_sc3 <- uniLasso(x = X_mat_sc3, y = Surv(simData_mod_sc3$Y, simData_mod_sc3$delta),
                               family = "cox", lambda = unilasso_cv_sc3$lambda.1se)
  t2 <- Sys.time()
  res_df$c_time_sc3_uni[i] <- as.numeric(t2 - t1)
  beta_nonZero <- names(unilasso_cox_sc3$beta[ , 1][unilasso_cox_sc3$beta[, 1] != 0])
  X1_cuts_ind <- as.numeric(substr(beta_nonZero[grepl("X1", beta_nonZero)], 8, 9))
  X2_cuts_ind <- as.numeric(substr(beta_nonZero[grepl("X2", beta_nonZero)], 8, 9))
  X1_cuts_ind <- X1_cuts_ind[c(1, which(diff(X1_cuts_ind) > 1) + 1)]
  X2_cuts_ind <- X2_cuts_ind[c(1, which(diff(X2_cuts_ind) > 1) + 1)]
  res_df$x1_cuts_sc3_uni[i] <- list(unique(x1_bounds_sc3[X1_cuts_ind]))
  res_df$x2_cuts_sc3_uni[i] <- list(unique(x2_bounds_sc3[X2_cuts_ind]))
  
  true_ind_tmp <- which(true_cuts_betas$sample_size == res_df$SS[i] &
                          true_cuts_betas$iteration == res_df$iter[i])
  true_vals_tmp <- true_cuts_betas[true_ind_tmp, ]
  Bussy_ind_tmp <- which(Bussy_cuts_betas$sample_size == res_df$SS[i] &
                           Bussy_cuts_betas$iteration == res_df$iter[i])
  Bussy_vals_tmp <- Bussy_cuts_betas[Bussy_ind_tmp, ]
  
  res_df$M1_sc1[i] <-
    get_m_1(cut_points_estimates = list(res_df$x1_cuts_sc1[i], res_df$x2_cuts_sc1[i]),
            cut_points = list(c(true_vals_tmp$X1_cuts[[1]]),
                              c(true_vals_tmp$X2_cuts[[1]])))
  Bussy_cuts_betas$M1_sc1[Bussy_ind_tmp] <-
    get_m_1(cut_points_estimates = list(Bussy_vals_tmp$X1_cuts_sc1[[1]],
                                        Bussy_vals_tmp$X2_cuts_sc1[[1]]),
            cut_points = list(c(true_vals_tmp$X1_cuts[[1]]),
                              c(true_vals_tmp$X2_cuts[[1]])))

  res_df$M1_sc1_uni[i] <-
    get_m_1(cut_points_estimates = list(res_df$x1_cuts_sc1_uni[i], res_df$x2_cuts_sc1_uni[i]),
            cut_points = list(c(true_vals_tmp$X1_cuts[[1]]),
                              c(true_vals_tmp$X2_cuts[[1]])))

  res_df$M1_sc3[i] <-
    get_m_1(cut_points_estimates = list(res_df$x1_cuts_sc3[i], res_df$x2_cuts_sc3[i]),
            cut_points = list(c(true_vals_tmp$X1_cuts[[1]] * 3),
                              c(true_vals_tmp$X2_cuts[[1]] * 3)))
  Bussy_cuts_betas$M1_sc3[Bussy_ind_tmp] <-
    get_m_1(cut_points_estimates = list(Bussy_vals_tmp$X1_cuts_sc3[[1]],
                                        Bussy_vals_tmp$X2_cuts_sc3[[1]]),
            cut_points = list(c(true_vals_tmp$X1_cuts[[1]] * 3),
                              c(true_vals_tmp$X2_cuts[[1]] * 3)))

  res_df$M1_sc3_uni[i] <-
    get_m_1(cut_points_estimates = list(res_df$x1_cuts_sc3_uni[i], res_df$x2_cuts_sc3_uni[i]),
            cut_points = list(c(true_vals_tmp$X1_cuts[[1]] * 3),
                              c(true_vals_tmp$X2_cuts[[1]] * 3)))

  cox_tmp_sc1 <-
    simData_sc1 %>%
    cox_cat_fit(simData = ., x1_cuts = res_df$x1_cuts_sc1[i],
                x2_cuts = res_df$x2_cuts_sc1[i])
  res_df$AIC_sc1[i] <- AIC(cox_tmp_sc1$coxModel)
  res_df$BIC_sc1[i] <- BIC(cox_tmp_sc1$coxModel)
  res_df$IBS_sc1[i] <- cox_tmp_sc1$cox_ibs
  res_df$betas_sc1[i] <- list(coef(cox_tmp_sc1$coxModel))
  cox_tmp_sc1_uni <-
    simData_sc1 %>%
    cox_cat_fit(simData = .,
                x1_cuts = res_df$x1_cuts_sc1_uni[i],
                x2_cuts = res_df$x2_cuts_sc1_uni[i])
  res_df$AIC_sc1_uni[i] <- AIC(cox_tmp_sc1_uni$coxModel)
  res_df$BIC_sc1_uni[i] <- BIC(cox_tmp_sc1_uni$coxModel)
  res_df$IBS_sc1_uni[i] <- cox_tmp_sc1_uni$cox_ibs
  res_df$betas_sc1_uni[i] <- list(coef(cox_tmp_sc1_uni$coxModel))
  
  cox_bussy_tmp_sc1 <-
    simData_sc1 %>%
    cox_cat_fit(simData = ., x1_cuts = Bussy_vals_tmp$X1_cuts_sc1,
                x2_cuts = Bussy_vals_tmp$X2_cuts_sc1)
  Bussy_cuts_betas$AIC_sc1[Bussy_ind_tmp] <- AIC(cox_bussy_tmp_sc1$coxModel)
  Bussy_cuts_betas$BIC_sc1[Bussy_ind_tmp] <- BIC(cox_bussy_tmp_sc1$coxModel)
  Bussy_cuts_betas$IBS_sc1[Bussy_ind_tmp] <- cox_bussy_tmp_sc1$cox_ibs
  Bussy_cuts_betas$betas_sc1[Bussy_ind_tmp] <- list(coef(cox_bussy_tmp_sc1$coxModel))
  cox_true_tmp_sc1 <-
    simData_sc1 %>%
    cox_cat_fit(simData = ., cont = TRUE)
  true_cuts_betas$AIC_sc1[true_ind_tmp] <- AIC(cox_true_tmp_sc1$coxModel)
  true_cuts_betas$BIC_sc1[true_ind_tmp] <- BIC(cox_true_tmp_sc1$coxModel)
  true_cuts_betas$IBS_sc1[true_ind_tmp] <- cox_true_tmp_sc1$cox_ibs
  true_cuts_betas$betas_est_sc1[true_ind_tmp] <- list(coef(cox_true_tmp_sc1$coxModel))
  
  cox_tmp_sc3 <-
    simData_sc3 %>%
    cox_cat_fit(simData = ., x1_cuts = res_df$x1_cuts_sc3[i],
                x2_cuts = res_df$x2_cuts_sc3[i])
  res_df$AIC_sc3[i] <- AIC(cox_tmp_sc3$coxModel)
  res_df$BIC_sc3[i] <- BIC(cox_tmp_sc3$coxModel)
  res_df$IBS_sc3[i] <- cox_tmp_sc3$cox_ibs
  res_df$betas_sc3[i] <- list(coef(cox_tmp_sc3$coxModel))
  cox_tmp_sc3_uni <-
    simData_sc3 %>%
    cox_cat_fit(simData = .,
                x1_cuts = res_df$x1_cuts_sc3_uni[i],
                x2_cuts = res_df$x2_cuts_sc3_uni[i])
  res_df$AIC_sc3_uni[i] <- AIC(cox_tmp_sc3_uni$coxModel)
  res_df$BIC_sc3_uni[i] <- BIC(cox_tmp_sc3_uni$coxModel)
  res_df$IBS_sc3_uni[i] <- cox_tmp_sc3_uni$cox_ibs
  res_df$betas_sc3_uni[i] <- list(coef(cox_tmp_sc3_uni$coxModel))
  
  cox_bussy_tmp_sc3 <-
    simData_sc3 %>%
    cox_cat_fit(simData = ., x1_cuts = Bussy_vals_tmp$X1_cuts_sc3,
                x2_cuts = Bussy_vals_tmp$X2_cuts_sc3)
  Bussy_cuts_betas$AIC_sc3[Bussy_ind_tmp] <- AIC(cox_bussy_tmp_sc3$coxModel)
  Bussy_cuts_betas$BIC_sc3[Bussy_ind_tmp] <- BIC(cox_bussy_tmp_sc3$coxModel)
  Bussy_cuts_betas$IBS_sc3[Bussy_ind_tmp] <- cox_bussy_tmp_sc3$cox_ibs
  Bussy_cuts_betas$betas_sc3[Bussy_ind_tmp] <- list(coef(cox_bussy_tmp_sc3$coxModel))
  cox_true_tmp_sc3 <-
    simData_sc3 %>%
    mutate(X1 = X1_new,
           X2 = X2_new) %>%
    cox_cat_fit(simData = ., cont = TRUE)
  true_cuts_betas$AIC_sc3[true_ind_tmp] <- AIC(cox_true_tmp_sc3$coxModel)
  true_cuts_betas$BIC_sc3[true_ind_tmp] <- BIC(cox_true_tmp_sc3$coxModel)
  true_cuts_betas$IBS_sc3[true_ind_tmp] <- cox_true_tmp_sc3$cox_ibs
  true_cuts_betas$betas_est_sc3[true_ind_tmp] <- list(coef(cox_true_tmp_sc3$coxModel))
  
  return(list(res_df = res_df[i, ],
              true_cuts_betas = true_cuts_betas[true_ind_tmp, ],
              Bussy_cuts_betas = Bussy_cuts_betas[Bussy_ind_tmp , ]))
}

