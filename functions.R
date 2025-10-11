cox_cat_fit <-
  function(simData,
           x1_cuts,
           x2_cuts,
           cont = FALSE) {
    
    if (cont) {
      simData %>%
        coxph(Surv(Y, delta) ~ X1 + X2, data = .,
              x = TRUE, y = TRUE) -> cox_tmp
      suppressMessages({ibs_tmp <- ibs(pec(list(Cox = cox_tmp), Hist(Y, delta) ~ 1, data = simData))[2, 1]})
      return(list(coxModel = cox_tmp,
                  cox_ibs = ibs_tmp))
    }
    
    if (any(is.na(x1_cuts[[1]]))) x1_cuts[[1]] <- numeric(0)
    if (any(is.na(x2_cuts[[1]]))) x2_cuts[[1]] <- numeric(0)
    
    if (length(x1_cuts[[1]]) != 0 & length(x2_cuts[[1]]) != 0) {
      simData %>%
        mutate(X1_bin = factor(cut(X1, breaks = c(-Inf, x1_cuts[[1]], Inf),
                                   labels = c(0 : length(x1_cuts[[1]])))),
               X2_bin = factor(cut(X2, breaks = c(-Inf, x2_cuts[[1]], Inf),
                                   labels = c(0 : length(x2_cuts[[1]]))))) -> simData_tmp
      cox_tmp <- coxph(Surv(Y, delta) ~ X1_bin + X2_bin, data = simData_tmp,
                       x = TRUE, y = TRUE)
    }
    if (length(x1_cuts[[1]]) == 0 & length(x2_cuts[[1]]) != 0) {
      simData %>%
        mutate(X2_bin = factor(cut(X2, breaks = c(-Inf, x2_cuts[[1]], Inf),
                                   labels = c(0 : length(x2_cuts[[1]]))))) -> simData_tmp
      cox_tmp <- coxph(Surv(Y, delta) ~ X2_bin, data = simData_tmp,
                       x = TRUE, y = TRUE) -> cox_tmp
    }
    if (length(x1_cuts[[1]]) != 0 & length(x2_cuts[[1]]) == 0) {
      simData %>%
        mutate(X1_bin = factor(cut(X1, breaks = c(-Inf, x1_cuts[[1]], Inf),
                                   labels = c(0 : length(x1_cuts[[1]]))))) -> simData_tmp
      cox_tmp <- coxph(Surv(Y, delta) ~ X1_bin, data = simData_tmp,
                       x = TRUE, y = TRUE)
    }
    if (length(x1_cuts[[1]]) == 0 & length(x2_cuts[[1]]) == 0) {
      simData_tmp <- simData
      cox_tmp <- coxph(Surv(Y, delta) ~ 1, data = simData_tmp,
                       x = TRUE, y = TRUE)
    }
    
    try(suppressMessages({ibs_tmp <- ibs(pec(list(Cox = cox_tmp), Hist(Y, delta) ~ 1, data = simData_tmp))[2, 1]}),
        silent = T)
    if (! exists("ibs_tmp")) ibs_tmp <- NA
    return(list(coxModel = cox_tmp,
                cox_ibs = ibs_tmp))
  }

get_H <- function(A, B) {
  max(get_E(A, B), get_E(B, A))
}

get_E <- function(A, B) {
  max(sapply(B, function(b) min(sapply(A, function(a) abs(a - b)))))
}

get_m_1 <- function(cut_points_estimates = list(X1_cuts, X2_cuts),
                    cut_points = list(c(-0.38828854, 0.57771615), c(-0.61267879, 0.60985693)),
                    S = list()) {
  m_1 <- 0
  d <- 0
  n_features <- length(cut_points)
  
  # j = 1
  for (j in setdiff(seq_len(n_features), S)) {  # Adjusts for zero-based indexing
    mu_star_j <- cut_points[[j]]
    hat_mu_star_j <- cut_points_estimates[[j]]
    
    if (all(! is.na(hat_mu_star_j)) & length(hat_mu_star_j) > 0) {
      d <- d + 1
      m_1 <- m_1 + get_H(mu_star_j, hat_mu_star_j)
    }
  }
  
  if (d == 0) {
    m_1 <- NA
  } else {
    m_1 <- m_1 / d
  }
  
  return(m_1)
}

cuts_to_list <-
  function(cuts_str = cuts_betas_org$cut_points_x[2]) {
    cuts_str <- strsplit(cuts_str, split = "array")[[1]][-1]
    lapply(cuts_str, function(x) {
      x_tmp <- strsplit(x, split = ",")[[1]][-1]
      if (any(grepl(":", x_tmp))) x_tmp <- x_tmp[-c(length(x_tmp) - 1, length(x_tmp))]
      else x_tmp <- x_tmp[-c(length(x_tmp))]
      x_tmp <- as.numeric(trimws(x_tmp))
    })
  }

cuts_finder <- function(x = X_mat_sc1,
                        y = Surv(simData_mod_sc1$Y, simData_mod_sc1$delta),
                        group_var = group_var,
                        uniLasso_flag = FALSE) {
  lambda_grid <- exp(seq(-10, 2, length.out = 1000))
  if (uniLasso_flag) {
    fit_init <- uniLasso(x, y, family = "cox", lambda = lambda_grid)
  }
  else {
    fit_init <- glmnet(x, y, family = "cox", lambda = lambda_grid)
  }
  coef_init <- as.matrix(coef(fit_init))
  nz_counts <- sapply(unique(group_var), function(g) {
    idx <- which(group_var == g)
    apply(coef_init[idx, , drop = FALSE], 2, function(col) sum(col != 0))
  })
  devs <- apply(nz_counts, 1,
                function(row) ifelse(row[1] >= 2 & row[2] >= 2,
                                     (row[1] - 2) ^ 2 + (row[2] - 2) ^ 2,
                                     Inf))
  lambda_idx <- which.min(devs)
  lambda_init <- fit_init$lambda[lambda_idx]
  if (uniLasso_flag) {
    fit_stage1 <- uniLasso(x, y , family = "cox", lambda = lambda_init,
                           penalty.factor = rep(1, ncol(x)))
  }
  else {
    fit_stage1 <- glmnet(x, y , family = "cox", lambda = lambda_init,
                         penalty.factor = rep(1, ncol(x)))
  }
  coef_stage1 <- coef(fit_stage1)
  nz_g1 <- sum(coef_stage1[group_var == unique(group_var)[1]] != 0)
  nz_g2 <- sum(coef_stage1[group_var == unique(group_var)[2]] != 0)
  
  orders_x1 <- order(abs(coef_stage1[c(1 : 50), 1]), decreasing = TRUE)
  if (nz_g1 > 2 & abs(orders_x1[1] - orders_x1[2]) == 1) orders_x1 = orders_x1[-2]
  cuts_x1 <- sapply(row.names(coef_stage1)[c(1 : 50)],
                    function(x) as.numeric(substr(x, 8, 9)))
  cuts_x1 <- cuts_x1[orders_x1[c(1 : 2)]]
  orders_x2 <- order(abs(coef_stage1[c(51 : 100), 1]), decreasing = TRUE)
  if (nz_g2 > 2 & abs(orders_x2[1] - orders_x2[2]) == 1) orders_x2 = orders_x2[-2]
  cuts_x2 <- sapply(row.names(coef_stage1)[c(51 : 100)],
                    function(x) as.numeric(substr(x, 8, 9)))
  cuts_x2 <- cuts_x2[orders_x2[c(1 : 2)]]
  
  return(list(cuts_x1 = cuts_x1, cuts_x2 = cuts_x2,
              fit_obj = fit_stage1))
}
