source('00_system.R')

load('Datasets/app_dataset.RData')

# Model B - Stationary

if(T){
  
  list_formulas <- list("mean" = as.formula("  ~  1 + wind + new_wind + BIO04 + BIO15 + cloud_c + elevation + lati + long"),
                        "std.dev" = as.formula(" ~ 1"),
                        "scale" = as.formula("  ~ 1"),
                        "aniso" = as.formula("  ~ 1"),
                        "tilt" = as.formula("  ~ 1"),
                        "smooth" = as.formula(" ~ 1"),
                        "nugget" = -Inf)
  
  test_coco_classic <- coco(type = 'dense',
                            data = all_dfs[[2]],
                            locs = as.matrix(all_dfs[[2]][, c("long","lati")]),
                            z = c(scale(all_dfs[[2]]$prec)),
                            model.list = list_formulas,
                            info = list('smooth.limits' = c(0.5, 2.0)
                                        )
                            )
  
  # Tailored boundaries
  if(T){
    
    boundaries_B <- getBoundariesV2(coco.object = test_coco_classic,
                                    mean.limits = c(-1.5, 0, 1.5),
                                    std.dev.limits = c(-2.5, 0, 2.5),
                                    scale.limits = c(-2.5, 0, 2.5),
                                    aniso.limits =  c(-2, 0, 2),
                                    tilt.limits =  c(-2, 0, 2),
                                    smooth.limits = c(-2, 0, 1.5),
                                    nugget.limits = c(-5, -3, 1))
    
    first_var <- which(names(boundaries_B$theta_init) == "std.dev.limits")[1]
    first_range <- which(names(boundaries_B$theta_init) == "scale.limits")[1]
    first_smooth <- which(names(boundaries_B$theta_init) == "smooth.limits")[1]
    
    boundaries_B$theta_upper[c(first_var, first_range)] <- c(5, 5)
    boundaries_B$theta_lower[c(first_var, first_range)] <- c(-5, -1.75)
    
    boundaries_B$theta_init[first_range] <- log(sd(c(dist(test_coco_classic@locs))))
    
    boundaries_B$theta_upper[first_smooth] <- 2
    boundaries_B$theta_lower[first_smooth] <- -3
    boundaries_B$theta_init[first_smooth] <- -3
    
    boundaries_B$theta_init[1] <- mean(test_coco_classic@z)
    boundaries_B$theta_upper[1] <- boundaries_B$theta_init[1] + 1
    boundaries_B$theta_lower[1] <- boundaries_B$theta_init[1] - 1
    
  }
  
  Time_B <- system.time({Model_B <- cocoOptim(coco.object = test_coco_classic,
                                              ncores = "auto",
                                              boundaries = boundaries_B,
                                              optim.type = 'ml')})
  
  HESS_B <- getHessian(Model_B, ncores = ncores)

  save(Time_B, Model_B, HESS_B, file = 'RData/Model_B.RData')
  
}

# Predictions

if(T){
  
  # Create holdouts
  if(T){
    
    newdataset <- all_dfs[[1]]
    newlocs <- all_dfs[[1]][ ,c("long","lati")]
    z_values <- all_dfs[[1]]$prec
    
    set.seed(100621)
    hetero_holdouts <- kmeans(as.data.frame(scale(newdataset[,c(1,2,4:9)])),centers = 100,iter.max = 100)
    groups <- as.factor(hetero_holdouts$cluster)
    quilt.plot(newlocs, hetero_holdouts$cluster, nx = 150, ny = 150)
    
    sample_to_tune_hyperparameters <- sample(1:100, 30)
    
    newdataset_final <- newdataset[!(groups %in% sample_to_tune_hyperparameters),]
    newlocs_final <- as.matrix(newdataset[!(groups %in% sample_to_tune_hyperparameters), c("long","lati")])
    z_values_final <- all_dfs[[1]]$prec[!(groups %in% sample_to_tune_hyperparameters)]
    
  }  
  
  load('RData/Model_A.RData')

  Pred_B <- cocoPredict(coco.object = Model_A, 
                         newdataset = newdataset_final, 
                         newlocs = newlocs_final,
                         type = 'pred')
  
  Pred_B <- cocons::cocoPredict(coco.object = Model_B, 
                                newdataset = newdataset_final, 
                                newlocs = newlocs_final,
                                type = 'pred')
  
  CRPS_B <- getCRPS(z_values_final,mean.pred = Pred_A$systematic + Pred_A$stochastic, sd.pred = Pred_A$sd.pred)
  CRPS_C <- getCRPS(z_values_final,mean.pred = Pred_B$systematic + Pred_B$stochastic, sd.pred = Pred_B$sd.pred)
  
  Logscore_A <- getLogRank(z_values, mean.pred = Pred_A$systematic + Pred_A$stochastic,sd.pred = Pred_A$sd.pred)
  Logscore_B <- getLogRank(z_values, mean.pred = Pred_B$systematic + Pred_B$stochastic,sd.pred = Pred_B$sd.pred)

  z_std_A <- (Pred_A$systematic + Pred_A$stochastic - z_values) / Pred_A$sd.pred
  z_std_B <- (Pred_B$systematic + Pred_B$stochastic - z_values) / Pred_B$sd.pred

  save(Pred_A, Pred_B,
       CRPS_A, CRPS_B, 
       Logscore_A, Logscore_B,
       z_std_A, z_std_B,
       file = 'RData/Pred_dense.RData')
  
}

# Metrics

load('RData/Pred_dense.RData')

if(T){
  
  # RMSPE
  
  if(T){
    
    RMSPE_A <- numeric(70)
    
    for(ii in 1:70){
      RMSPE_A[ii] <- sqrt(mean((Pred_C$mean[groups == levels(groups)[ii]] + Pred_C$trend[groups == levels(groups)[ii]] - z_values[groups == levels(groups)[ii]])^2))
    }
    
    mean(RMSPE_A)
    sd(RMSPE_A)
    
    RMSPE_B <- numeric(70)
    
    for(ii in 1:70){
      RMSPE_B[ii] <- sqrt(mean((Pred_D$mean[groups == levels(groups)[ii]] + Pred_D$trend[groups == levels(groups)[ii]] - z_values[groups == levels(groups)[ii]])^2))
    }
    
    mean(RMSPE_B)
    sd(RMSPE_B)
    
    matrix_to_show <- t(round(cbind(c(mean(RMSPE_A),mean(RMSPE_B)),
                                    c(sd(RMSPE_A),sd(RMSPE_B))),2))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS')
    
    cat('RMSPE\n')
    print(matrix_to_show)
    cat('------\n')
    
  }
  
  # CRPS
  
  if(T){
    
    CRPS_A_holdout <- numeric(70)
    
    for(ii in 1:70){
      CRPS_A_holdout[ii] <- mean(CRPS_A[groups == levels(groups)[ii]])
    }
    
    mean(CRPS_A_holdout)
    sd(CRPS_A_holdout)
    
    CRPS_B_holdout <- numeric(70)
    
    for(ii in 1:70){
      CRPS_B_holdout[ii] <- mean(CRPS_B[groups == levels(groups)[ii]])
    }
    
    mean(CRPS_B_holdout)
    sd(CRPS_B_holdout)
    
    matrix_to_show <- t(round(cbind(c(mean(CRPS_D_holdout),mean(CRPS_C_holdout)),
                                    c(sd(CRPS_D_holdout),sd(CRPS_C_holdout))), 3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS','M-NS-G')
    
    cat('CRPS\n')
    print(matrix_to_show)
    cat('------\n')
    
  }
  
  # Q_0.9 CRPS
  
  if(T){
    
    Q_CRPS_A <- numeric(70)
    
    for(ii in 1:70){
      Q_CRPS_A[ii] <- quantile(CRPS_A[groups == levels(groups)[ii]],probs = 0.95)
    }
    
    mean(Q_CRPS_A)
    sd(Q_CRPS_A)
    
    Q_CRPS_B <- numeric(70)
    
    for(ii in 1:70){
      Q_CRPS_B[ii] <- quantile(CRPS_B[groups == levels(groups)[ii]],probs = 0.95)
    }
    
    mean(Q_CRPS_B)
    sd(Q_CRPS_B)
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS')
    
    cat('0.95 Quantile CRPS\n')
    print(matrix_to_show)
    cat('---------\n')
    
  }
  
  # D_n
  
  if(T){
    
    KS_A <- numeric(70)
    
    for(ii in 1:70){
      KS_A[ii] <- ks.test(z_std_A[groups == levels(groups)[ii]],y = pnorm)$statistic
    }
    
    mean(KS_A)
    sd(KS_A)
    
    KS_B <- numeric(70)
    
    for(ii in 1:70){
      KS_B[ii] <- ks.test(z_std_B[groups == levels(groups)[ii]],
                          y = pnorm)$statistic
    }
    
    mean(KS_B)
    sd(KS_B)  
    
    matrix_to_show <- t(round(cbind(c(mean(KS_A),mean(KS_B)),
                                    c(sd(KS_A),sd(KS_B))),2))
    
    rownames(matrix_to_show) <- c('mean', 'se')
    colnames(matrix_to_show) <- c('M-STAT', 'M-NS')
    
    cat('KS statistic \n')
    print(matrix_to_show)
    cat('--------\n')
    
  }
  
  # Coverage_prob
  
  if(T){
    
    Cov_prob_A <- numeric(70)
    
    for(ii in 1:70){
      
      upper_bound <- Pred_A$trend[groups == levels(groups)[ii]] + 
        Pred_A$mean[groups == levels(groups)[ii]] +
        qnorm(1 - 0.025) * Pred_A$sd.pred[groups == levels(groups)[ii]]
      
      lower_bound <- Pred_A$trend[groups == levels(groups)[ii]] + 
        Pred_A$mean[groups == levels(groups)[ii]] -
        qnorm(1 - 0.025) * Pred_A$sd.pred[groups == levels(groups)[ii]]
      
      Cov_prob_A[ii] <-  1 - length(which(z_values[groups == levels(groups)[ii]] < lower_bound |
                                            z_values[groups == levels(groups)[ii]] > upper_bound  )) /
        length(z_values[groups == levels(groups)[ii]])
      
    }
    
    mean(Cov_prob_A)
    sd(Cov_prob_A)
    
    Cov_prob_B <- numeric(70)
    
    for(ii in 1:70){
      
      upper_bound <- Pred_B$trend[groups == levels(groups)[ii]] + 
        Pred_B$mean[groups == levels(groups)[ii]] +
        qnorm(1 - 0.025) * Pred_B$sd.pred[groups == levels(groups)[ii]]
      
      lower_bound <- Pred_B$trend[groups == levels(groups)[ii]] + 
        Pred_B$mean[groups == levels(groups)[ii]] -
        qnorm(1 - 0.025) * Pred_B$sd.pred[groups == levels(groups)[ii]]
      
      Cov_prob_B[ii] <-  1 - length(which(z_values[groups == levels(groups)[ii]] < lower_bound |
                                            z_values[groups == levels(groups)[ii]] > upper_bound  )) /
        length(z_values[groups == levels(groups)[ii]])
      
      
    }
    
    mean(Cov_prob_B)
    sd(Cov_prob_B)
    
    matrix_to_show <- t(round(cbind(c(mean(Cov_prob_A),mean(Cov_prob_B)),
                                    c(sd(Cov_prob_A),sd(Cov_prob_B))),2))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS','M-NS-G')
    
    
    cat('coverage probability\n')
    print(matrix_to_show)
    cat('--------------\n')
    
  }
  
  # log-score
  
  if(T){
    
    LS_A <- numeric(70)
    
    for(ii in 1:70){
      LS_A[ii] <- mean(Logscore_A[groups == levels(groups)[ii]])
    }
    
    mean(LS_A)
    sd(LS_A)
    
    LS_B <- numeric(70)
    
    for(ii in 1:70){
      LS_B[ii] <-  mean(Logscore_D[groups == levels(groups)[ii]])
    }
    
    matrix_to_show <- t(round(cbind(c(mean(LS_A),mean(LS_B)),
                                    c(sd(LS_A),sd(LS_B))),2))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS')
    
    cat('logscore\n')
    print(matrix_to_show)
    cat('--------\n')
    
  }
  
  # logliks
  
  vector_to_show <- c(getLoglik(Model_A),getLoglik(Model_B))
  names(vector_to_show) <- c('M-STAT','M-NS')
  
  cat('Logliks\n')
  print(vector_to_show)
  cat('--------\n')
  
  # times
  
  vector_to_show <- round(c(Time_D[3]/60,time_C[3]/60),2)
  names(vector_to_show) <- c('M-STAT','M-NS')
  
  cat('system times\n')
  print(vector_to_show)
  cat('--------\n')
  
}
