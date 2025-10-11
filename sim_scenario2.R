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
sample_size <- 1000
n_featr_s <- c(4, 6, 10, 20, 40, 60, 80, 100)

## Load simulated data as well as results for Binacox
for (i in 1 : length(n_featr_s)) {
  p <- n_featr_s[i]
  cuts_betas_org_tmp <-
    read.csv(paste0("simData/Scen2/Scen2_nfeatrs_", p, "_cutPoints_beta.csv"),
             header = TRUE, row.names = NULL)
  true_cuts_betas_tmp <-
    cuts_betas_org_tmp[-1, ] %>%
    select(sample_size, n_features, iteration, cut_points_x, beta_star) %>%
    filter(n_features %in% n_featr_s) %>%
    rowwise %>%
    mutate(X_cuts = list(cuts_to_list(cut_points_x)),
           betas = list(beta_star)) %>%
    select(sample_size, n_features, iteration, X_cuts, betas) %>%
    ungroup
  if (i == 1) true_cuts_betas <- cuts_betas_org_tmp
  else {
    true_cuts_betas %>%
      bind_rows(true_cuts_betas)
  }
  Bussy_cuts_betas_tmp <-
    cuts_betas_org_tmp[-1, ] %>%
    select(sample_size, n_features, iteration,
           `Processing_Time..seconds.`, cut_points_estimate_x) %>%
    filter(n_features %in% n_featr_s) %>%
    rowwise %>%
    mutate(X_cuts_est = list(cuts_to_list(cut_points_estimate_x))) %>%
    rename(c_time = `Processing_Time..seconds.`) %>%
    select(sample_size, n_features, iteration, X_cuts_est, c_time) %>%
    ungroup
  if (i == 1) Bussy_cuts_betas <- Bussy_cuts_betas_tmp
  else {
    Bussy_cuts_betas %>%
      bind_rows(Bussy_cuts_betas)
  }
  rm(cuts_betas_org_tmp)
  rm(cuts_betas_org_tmp)
  rm(Bussy_cuts_betas_tmp)
}
file_names <- list.files("simData/Scen2")
file_names <- file_names[- which(grepl("Scen2_nfeatrs_", file_names))]
file_names_filter_tmp <-
  sapply(file_names,
         function(x){
           nFeat_tmp <- as.numeric(strsplit(x, split = "_")[[1]][2])
           if (! nFeat_tmp %in% n_featr_s) nFeat_tmp = NA
           return(nFeat_tmp)
         })
file_names_filter <- names(file_names_filter_tmp[! is.na(file_names_filter_tmp)])
file_names_filter <- file_names

## Fit models using biniLasso & miniLasso
res_df <- tibble(SS = NA,
                 iter = NA,
                 n_features = NA,
                 c_time = NA,
                 AIC = NA,
                 BIC = NA,
                 IBS = NA,
                 M1 = NA,
                 c_time_uni = NA,
                 AIC_uni = NA,
                 BIC_uni = NA,
                 IBS_uni = NA,
                 M1_uni = NA) %>%
  mutate(x_cuts = list(list(NULL)),
         betas = list(list(NULL)),
         x_cuts_uni = list(list(NULL)),
         betas_uni = list(list(NULL)))
Bussy_cuts_betas %<>%
  mutate(AIC = NA,
         BIC = NA,
         IBS = NA,
         M1 = NA,
         betas = list(list(NULL)),
         betas_uni = list(list(NULL)))

true_cuts_betas$AIC <- true_cuts_betas$BIC <- true_cuts_betas$IBS <- NA
true_cuts_betas$betas_est <- list(list(NULL))

for (i in 1 : nrow(true_cuts_betas)) {
  file_tmp <- paste0("simData/Scen2/", file_names_filter[i])
  simData <- read.csv(file_tmp,
                      header = TRUE, row.names = NULL)
  simData_sc1 <-
    simData %>%
    select(Y, delta, colnames(simData)[grepl("X", colnames(simData))])
  rm(simData)
  
  vals_tmp <- strsplit(file_tmp, "_")[[1]]
  res_df$SS <- sample_size
  res_df$n_features <- as.numeric(vals_tmp[3])
  res_df$iter <- as.numeric(str_sub(vals_tmp[4], start = 1, end = nchar(vals_tmp[4]) - 4))
  
  x1_bounds_sc1 <- quantile(simData_sc1$X1, probs = c(1 : n_bins) / (1 + n_bins))
  x2_bounds_sc1 <- quantile(simData_sc1$X2, probs = c(1 : n_bins) / (1 + n_bins))
  
  for (nf in 1 : res_df$n_features) {
    x_bounds_tmp <- quantile(simData_sc1[ , paste0("X", nf - 1)],
                             probs = c(1 : n_bins) / (1 + n_bins))
    for (j in 1 : n_bins) {
      simData_sc1$tmpCol <- simData_sc1[ , paste0("X", nf - 1)]
      simData_sc1 %<>%
        mutate(temp = ifelse(tmpCol >= x_bounds_tmp[j], 1, 0))
      if (j < 10){
        colnames(simData_sc1)[colnames(simData_sc1) == "temp"] <- paste0("X", nf - 1,"_bin_0", j)
      }
      else {
        colnames(simData_sc1)[colnames(simData_sc1) == "temp"] <- paste0("X", nf - 1, "_bin_", j)
      }
    }
    
    if (nf == 1) x_min_max <- list(c(min(simData_sc1[ , paste0("X", nf - 1)]),
                                     max(simData_sc1[ , paste0("X", nf - 1)])))
    if (nf == 2) x_min_max <- list(c(x_min_max, list(c(min(simData_sc1[ , paste0("X", nf - 1)]),
                                                       max(simData_sc1[ , paste0("X", nf - 1)])))))
    if (nf > 2) x_min_max <- list(c(x_min_max[[1]], list(c(min(simData_sc1[ , paste0("X", nf - 1)]),
                                                           max(simData_sc1[ , paste0("X", nf - 1)])))))
  }
  
  simData_sc1_mod <-
    simData_sc1 %>%
    dplyr::select(Y, delta) %>%
    bind_cols(simData_sc1[ , colnames(simData_sc1)[grepl("_bin_", colnames(simData_sc1))]])
  X_mat_sc1 <-
    model.matrix( ~ . -1,
                  data = simData_sc1_mod[ , colnames(simData_sc1_mod)[grepl("_bin_", colnames(simData_sc1_mod))]])
  
  t1 <- Sys.time()
  glm_cv_sc1 <- cv.glmnet(x = X_mat_sc1, y = Surv(simData_sc1_mod$Y,
                                                  simData_sc1_mod$delta),
                          family = "cox", nfolds = 4)
  glm_cox_sc1 <- glmnet(x = X_mat_sc1, y = Surv(simData_sc1_mod$Y, simData_sc1_mod$delta),
                        family = "cox", lambda = glm_cv_sc1$lambda.1se)
  t2 <- Sys.time()
  res_df$c_time <- as.numeric(t2 - t1)
  beta_nonZero_sc1 <- names(glm_cox_sc1$beta[ , 1][glm_cox_sc1$beta[, 1] != 0])
  for (nf in 0 : (res_df$n_features - 1)) {
    X_cuts_ind_tmp <- as.numeric(unlist(lapply(beta_nonZero_sc1[grepl(paste0("X", nf), beta_nonZero_sc1)], function(x) strsplit(x, "_")[[1]][3])))
    X_cuts_ind_tmp <- X_cuts_ind_tmp[c(1, which(diff(X_cuts_ind_tmp) > 1) + 1)]
    x_bounds_tmp <- quantile(simData_sc1[ , paste0("X", nf)],
                             probs = c(1 : n_bins) / (1 + n_bins))
    if (nf == 0) res_df$x_cuts <- list(unique(x_bounds_tmp[X_cuts_ind_tmp]))
    if (nf == 1) res_df$x_cuts <- list(c(res_df$x_cuts, list(unique(x_bounds_tmp[X_cuts_ind_tmp]))))
    if (nf > 1) res_df$x_cuts <- list(c(res_df$x_cuts[[1]], list(unique(x_bounds_tmp[X_cuts_ind_tmp]))))
  }
  
  t1 <- Sys.time()
  unilasso_cv_sc1 <- cv.glmnet(x = X_mat_sc1, y = Surv(simData_sc1_mod$Y,
                                                       simData_sc1_mod$delta),
                               family = "cox", nfolds = 4)
  unilasso_cox_sc1 <- glmnet(x = X_mat_sc1, y = Surv(simData_sc1_mod$Y, simData_sc1_mod$delta),
                             family = "cox", lambda = unilasso_cv_sc1$lambda.1se)
  t2 <- Sys.time()
  res_df$c_time_uni <- as.numeric(t2 - t1)
  beta_nonZero_sc1 <- names(unilasso_cox_sc1$beta[ , 1][unilasso_cox_sc1$beta[, 1] != 0])
  for (nf in 0 : (res_df$n_features - 1)) {
    X_cuts_ind_tmp <- as.numeric(unlist(lapply(beta_nonZero_sc1[grepl(paste0("X", nf), beta_nonZero_sc1)], function(x) strsplit(x, "_")[[1]][3])))
    X_cuts_ind_tmp <- X_cuts_ind_tmp[c(1, which(diff(X_cuts_ind_tmp) > 1) + 1)]
    x_bounds_tmp <- quantile(simData_sc1[ , paste0("X", nf)],
                             probs = c(1 : n_bins) / (1 + n_bins))
    if (nf == 0) res_df$x_cuts_uni <- list(unique(x_bounds_tmp[X_cuts_ind_tmp]))
    if (nf == 1) res_df$x_cuts_uni <- list(c(res_df$x_cuts_uni, list(unique(x_bounds_tmp[X_cuts_ind_tmp]))))
    if (nf > 1) res_df$x_cuts_uni <- list(c(res_df$x_cuts_uni[[1]], list(unique(x_bounds_tmp[X_cuts_ind_tmp]))))
  }
  
  true_ind_tmp <- which(true_cuts_betas$sample_size == res_df$SS &
                          true_cuts_betas$n_features == res_df$n_features &
                          true_cuts_betas$iteration == res_df$iter)
  true_vals_tmp <- true_cuts_betas[true_ind_tmp, ]
  Bussy_ind_tmp <- which(Bussy_cuts_betas$sample_size == res_df$SS &
                           Bussy_cuts_betas$n_features == res_df$n_features &
                           Bussy_cuts_betas$iteration == res_df$iter)
  Bussy_vals_tmp <- Bussy_cuts_betas[Bussy_ind_tmp, ]
  
  res_df$M1 <-
    get_m_1(cut_points_estimates = res_df$x_cuts[[1]],
            cut_points = true_vals_tmp$X_cuts[[1]])
  Bussy_vals_tmp$M1 <-
    get_m_1(cut_points_estimates = Bussy_vals_tmp$X_cuts_est[[1]],
            cut_points = true_vals_tmp$X_cuts[[1]])

  res_df$M1_uni <-
    get_m_1(cut_points_estimates = res_df$x_cuts[[1]],
            cut_points = true_vals_tmp$X_cuts[[1]])

  cox_tmp <-
    cox_cat_fit(simData_sc1, x_cuts = res_df$x_cuts[[1]])
  res_df$AIC <- AIC(cox_tmp$coxModel)
  res_df$BIC <- BIC(cox_tmp$coxModel)
  res_df$IBS <- cox_tmp$cox_ibs
  res_df$betas <- list(coef(cox_tmp$coxModel))
  cox_tmp_uni <-
    cox_cat_fit(simData_sc1, x_cuts = res_df$x_cuts_uni[[1]])
  res_df$AIC_uni <- AIC(cox_tmp_uni$coxModel)
  res_df$BIC_uni <- BIC(cox_tmp_uni$coxModel)
  res_df$IBS_uni <- cox_tmp_uni$cox_ibs
  res_df$betas_uni <- list(coef(cox_tmp_uni$coxModel))
  
  cox_bussy_tmp <-
    cox_cat_fit(simData_sc1, x_cuts = Bussy_vals_tmp$X_cuts_est[[1]])
  Bussy_vals_tmp$AIC <- AIC(cox_bussy_tmp$coxModel)
  Bussy_vals_tmp$BIC <- BIC(cox_bussy_tmp$coxModel)
  Bussy_vals_tmp$IBS <- cox_bussy_tmp$cox_ibs
  Bussy_vals_tmp$betas <- list(coef(cox_bussy_tmp$coxModel))
  cox_true_tmp <-
    cox_cat_fit(simData_sc1, x_cuts = true_vals_tmp$X_cuts[[1]])
  true_vals_tmp$AIC <- AIC(cox_true_tmp$coxModel)
  true_vals_tmp$BIC <- BIC(cox_true_tmp$coxModel)
  true_vals_tmp$IBS <- cox_true_tmp$cox_ibs
  true_vals_tmp$betas_est <- list(coef(cox_true_tmp$coxModel))
  
  return(list(res_df = res_df,
              true_cuts_betas = true_vals_tmp,
              Bussy_cuts_betas = Bussy_vals_tmp))
}

