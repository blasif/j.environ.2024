
source('00_system.R')
load('Datasets/app_dataset.RData')

# First model

if(T){
  
  list_formulas <- list("mean" = as.formula("  ~ 1 + wind + new_wind + BIO04 + BIO15 + cloud_c + elevation + lati + long"),
                        "std.dev" = as.formula("  ~ 1 + log_cloud + log_elevation + log_wind + log_lati + log_long"),
                        "scale" = as.formula("  ~ 1 + log_cloud + log_elevation + log_wind + log_lati + log_long"),
                        "aniso" = 0,
                        "tilt" = 0,
                        "smooth" = as.formula("  ~ 1 + log_cloud + log_elevation + log_wind + log_lati + log_long"),
                        "nugget" = -Inf)
  
  test_coco <- coco(type = 'sparse',
                      data = all_dfs[[3]],
                      locs = as.matrix(all_dfs[[3]][, 1:2]),
                      z = all_dfs[[3]]$prec,
                      model.list = list_formulas,
                      info = list('lambda' = 0.0125 * sqrt(8),
                                  'taper' = spam::cov.wend2,
                                  'delta' = 0.23, 
                                  'smooth_limits' = c(0.5, 2.5))
  )
  
  # Tailored boundaries
  if(T){
    boundaries_A <- getBoundariesV2(coco.object = test_coco,
                                       mean.limits = c(-2, 0, 2),
                                       std.dev.limits = c(-2.5, 0, 2.5),
                                       scale.limits = c(-2, 0, 2),
                                       aniso.limits =  c(-2, 0, 2),
                                       tilt.limits =  c(-2, 0, 2),
                                       smooth.limits = c(-2, 0, 2),
                                       nugget.limits = c(NA, NA, NA))
    
    boundaries_A$theta_init[1] <- mean(test_coco@z)
    boundaries_A$theta_upper[1] <- boundaries_A$theta_init[1] + 2
    boundaries_A$theta_lower[1] <- boundaries_A$theta_init[1] - 2
    
    first_var <- which(names(boundaries_A$theta_init) == "std.dev.limits")[1]
    first_range <- which(names(boundaries_A$theta_init) == "scale.limits")[1]
    first_smooth <- which(names(boundaries_A$theta_init) == "smooth.limits")[1]
    
    boundaries_A$theta_upper[c(first_var, first_range)] <- c(2, 4)
    boundaries_A$theta_lower[c(first_var, first_range)] <- c(-5.5, -7.5)
    
    boundaries_A$theta_upper[first_smooth] <- 2
    boundaries_A$theta_lower[first_smooth] <- -2.5
    boundaries_A$theta_init[first_smooth] <- -2.5
  }
  
  time_T_A <- system.time({Model_T_A <- cocoOptim(coco.object = test_coco,
                                                   ncores = ncores,
                                                   boundaries = boundaries_A,
                                                   optim.type = 'mle')})
  
  save(Model_T_A, time_T_A, file = 'RData/Model_T_A.RData')
  
}

# Classic model

if(T){
  
  list_formulas <- list("mean" = as.formula("  ~ 1 + wind + new_wind + BIO04 + BIO15 + cloud_c + 
                                            elevation + lati + long"),
                        "std.dev" = as.formula("~ 1"),
                        "scale" = as.formula("~ 1"),
                        "aniso" = 0,
                        "tilt" = 0,
                        "smooth" = as.formula("~ 1"),
                        "nugget" = -Inf)
  
  OT_classic <- coco(type = 'sparse',
                      data = all_dfs[[3]],
                      locs = as.matrix(all_dfs[[3]][, 1:2]),
                      z = all_dfs[[3]]$prec,
                      model.list = list_formulas,
                      info = list('lambda' = 0.0125 * sqrt(8),
                                  'taper' = spam::cov.wend1,
                                  'delta' = 0.23, 
                                  'smooth_limits' = c(0.25, 1.5))
  )
  
  # Tailored boundaries
  
  if(T){
    
    Boundaries_B <- getBoundariesV2(coco.object = OT_classic,
                                    mean.limits = c(-2, 0, 2),
                                    std.dev.limits = c(-2.5, 0, 2.5),
                                    scale.limits = c(-2, 0, 2),
                                    aniso.limits =  c(-2, 0, 2),
                                    tilt.limits =  c(-2, 0, 2),
                                    smooth.limits = c(-2, 0, 2),
                                    nugget.limits = c(NA, NA, NA))
    
    Boundaries_B$theta_init[1] <- mean(OT_classic@z)
    Boundaries_B$theta_upper[1] <- Boundaries_B$theta_init[1] + 2
    Boundaries_B$theta_lower[1] <- Boundaries_B$theta_init[1] - 2
    
    first_var <- which(names(Boundaries_B$theta_init) == "std.dev.limits")[1]
    first_range <- which(names(Boundaries_B$theta_init) == "scale.limits")[1]
    first_smooth <- which(names(Boundaries_B$theta_init) == "smooth.limits")[1]
    
    Boundaries_B$theta_upper[c(first_var, first_range)] <- c(2, 4)
    Boundaries_B$theta_lower[c(first_var, first_range)] <- c(-4, -6)
    
    Boundaries_B$theta_upper[first_smooth] <- 2
    Boundaries_B$theta_lower[first_smooth] <- -2.5
    Boundaries_B$theta_init[first_smooth] <- -2.5
  }
  
  Time_T_B <- system.time({ Model_T_B <- cocoOptim(coco.object = OT_classic,
                                                    ncores = 4,
                                                    optim.type = 'mle',
                                                    boundaries = Boundaries_B)})
  
  save(Model_T_B, Time_T_B, file = 'RData/Model_T_B.RData')
  
}

# Global smoothness based on reduced model

if(T){
  
  list_formulas <- list("mean" = as.formula("  ~ 1 + wind + new_wind + BIO04 + BIO15 + cloud_c + elevation + lati + long"),
                        "std.dev" = as.formula("  ~ 1 + log_cloud + log_elevation + log_wind + log_lati + log_long"),
                        "scale" = as.formula("  ~ 1 + log_cloud + log_elevation + log_wind + log_lati + log_long"),
                        "aniso" = 0,
                        "tilt" = 0,
                        "smooth" = as.formula("  ~ 1"),
                        "nugget" = -Inf)
  
  test_coco <- coco(type = 'sparse',
                    data = all_dfs[[3]],
                    locs = as.matrix(all_dfs[[3]][, 1:2]),
                    z = all_dfs[[3]]$prec,
                    model.list = list_formulas,
                    info = list('lambda' = 0.0125 * sqrt(8),
                                'taper' = spam::cov.wend2,
                                'delta' = 0.23,
                                'smooth_limits' = c(0.5, 2.5))
  )
  
  # Boundaries
  
  if(T){
    
    boundaries_taper_C <- getBoundariesV2(coco.object = test_coco,
                                          mean.limits = c(-2, 0, 2),
                                          std.dev.limits = c(-2.5, 0, 2.5),
                                          scale.limits = c(-2, 0, 2),
                                          aniso.limits =  c(-2, 0, 2),
                                          tilt.limits =  c(-2, 0, 2),
                                          smooth.limits = c(-2, 0, 2),
                                          nugget.limits = c(NA, NA, NA))
    
    boundaries_taper_C$theta_init[1] <- mean(test_coco@z)
    boundaries_taper_C$theta_upper[1] <- boundaries_taper_C$theta_init[1] + 2
    boundaries_taper_C$theta_lower[1] <- boundaries_taper_C$theta_init[1] - 2
    
    first_var <- which(names(boundaries_taper_C$theta_init) == "std.dev.limits")[1]
    first_range <- which(names(boundaries_taper_C$theta_init) == "scale.limits")[1]
    first_smooth <- which(names(boundaries_taper_C$theta_init) == "smooth.limits")[1]
    
    boundaries_taper_C$theta_upper[c(first_var, first_range)] <- c(2, 4)
    boundaries_taper_C$theta_lower[c(first_var, first_range)] <- c(-5.5, -7.5)
    
    boundaries_taper_C$theta_upper[first_smooth] <- 2
    boundaries_taper_C$theta_lower[first_smooth] <- -2.5
    boundaries_taper_C$theta_init[first_smooth] <- -2.5
    
  }
  
  Time_T_C <- system.time({Model_T_C <- cocoOptim(coco.object = test_coco,
                                                  ncores = ncores,
                                                  boundaries = boundaries_taper_C,
                                                  optim.type = 'mle')})
  
  save(Model_T_C, Time_T_C, file = 'RData/Model_T_C.RData')
  
}

# Predictions
if(T){
  
  newdataset <- all_dfs[[1]]
  newlocs <- all_dfs[[1]][,1:2]
  
  z_values <- all_dfs[[1]]$prec
  
  # Create holdouts
  if(T){
    set.seed(100621)
    hetero_holdouts <- kmeans(as.data.frame(scale(newdataset[,c(1,2,4:9)])),centers = 100,iter.max = 100)
    groups <- as.factor(hetero_holdouts$cluster)
    quilt.plot(newlocs, hetero_holdouts$cluster, nx = 150, ny = 150)
  }

  Pred_Model_T_A <- cocoPredict(coco.object = Model_T_A, 
                                newdataset = newdataset, 
                                newlocs = newlocs,
                                type = 'pred')
  
  Pred_Model_T_B <- cocoPredict(coco.object = Model_T_B, 
                                newdataset = newdataset, 
                                newlocs = newlocs,
                                type = 'pred')
  
  Pred_Model_T_C <- cocoPredict(coco.object = Model_T_C, 
                                newdataset = newdataset, 
                                newlocs = newlocs,
                                type = 'pred')
  
  CRPS_T_A <- getCRPS(z_values,mean.pred = Pred_Model_T_A$trend + Pred_Model_T_A$mean,sd.pred = Pred_Model_T_A$sd.pred)
  CRPS_T_B <- getCRPS(z_values,mean.pred = Pred_Model_T_B$trend + Pred_Model_T_B$mean,sd.pred = Pred_Model_T_B$sd.pred)
  CRPS_T_C <- getCRPS(z_values,mean.pred = Pred_Model_T_C$trend + Pred_Model_T_C$mean,sd.pred = Pred_Model_T_C$sd.pred)
  
  Logscore_T_A <- getLogRank(z_values,mean.pred = Pred_Model_T_A$trend + Pred_Model_T_A$mean,sd.pred = Pred_Model_T_A$sd.pred)
  Logscore_T_B <- getLogRank(z_values,mean.pred = Pred_Model_T_B$trend + Pred_Model_T_B$mean,sd.pred = Pred_Model_T_B$sd.pred)
  Logscore_T_C <- getLogRank(z_values,mean.pred = Pred_Model_T_C$trend + Pred_Model_T_C$mean,sd.pred = Pred_Model_T_C$sd.pred)
  
  z_std_A <- (Pred_Model_T_A$trend + Pred_Model_T_A$mean - z_values) / Pred_Model_T_A$sd.pred
  z_std_B <- (Pred_Model_T_B$trend + Pred_Model_T_B$mean - z_values) / Pred_Model_T_B$sd.pred
  z_std_C <- (Pred_Model_T_C$trend + Pred_Model_T_C$mean - z_values) / Pred_Model_T_C$sd.pred
  
  save(CRPS_T_A, CRPS_T_B, CRPS_T_C,
       Logscore_T_A, Logscore_T_B, Logscore_T_C,
       z_std_A,z_std_B,z_std_C,
       Pred_Model_T_A,Pred_Model_T_B, Pred_Model_T_C, 
       file = 'RData/Pred_Taper.RData')
  
}
  
# Metrics

load('RData/Pred_Taper.RData')
  
if(T){
  
  # RMSPE
  
  if(T){
    
    RMSPE_T_A <- numeric(100)
    
    for(ii in 1:100){
      RMSPE_T_A[ii] <- sqrt(mean((Pred_Model_T_A$mean[groups == levels(groups)[ii]] + 
                                    Pred_Model_T_A$trend[groups == levels(groups)[ii]] - 
                                    z_values[groups == levels(groups)[ii]])^2))
    }
  
    RMSPE_T_B <- numeric(100)
    
    for(ii in 1:100){
      RMSPE_T_B[ii] <- sqrt(mean((Pred_Model_T_B$mean[groups == levels(groups)[ii]] + 
                                    Pred_Model_T_B$trend[groups == levels(groups)[ii]] - 
                                    z_values[groups == levels(groups)[ii]])^2))
    }

    RMSPE_T_C <- numeric(100)
    
    for(ii in 1:100){
      RMSPE_T_C[ii] <- sqrt(mean((Pred_Model_T_C$mean[groups == levels(groups)[ii]] + 
                                    Pred_Model_T_C$trend[groups == levels(groups)[ii]] - 
                                    z_values[groups == levels(groups)[ii]])^2))
    }

    matrix_to_show <- t(round(cbind(c(mean(RMSPE_T_B), mean(RMSPE_T_A), mean(RMSPE_T_C)),
                                    c(sd(RMSPE_T_B), sd(RMSPE_T_A), sd(RMSPE_T_C))), 3))
    
    rownames(matrix_to_show) <- c('mean', 'se')
    colnames(matrix_to_show) <- c('M-STAT-T', 'T-NS-T', 'M-NS-T-G')
    
    cat('RMSPE\n')
    print(matrix_to_show)
    cat('------\n')
    
  }
  
  # CRPS
  
  if(T){
    
    CRPS_A_holdout <- numeric(100)
    
    for(ii in 1:100){
      CRPS_A_holdout[ii] <- mean(CRPS_T_A[groups == levels(groups)[ii]])
    }

    CRPS_B_holdout <- numeric(100)
    
    for(ii in 1:100){
      CRPS_B_holdout[ii] <- mean(CRPS_T_B[groups == levels(groups)[ii]])
    }

    CRPS_C_holdout <- numeric(100)
    
    for(ii in 1:100){
      CRPS_C_holdout[ii] <- mean(CRPS_T_C[groups == levels(groups)[ii]])
    }

    matrix_to_show <- t(round(cbind(c(mean(CRPS_B_holdout),mean(CRPS_A_holdout),mean(CRPS_C_holdout)),
                                    c(sd(CRPS_B_holdout),sd(CRPS_A_holdout),sd(CRPS_C_holdout))),3))
    
    rownames(matrix_to_show) <- c('mean', 'se')
    colnames(matrix_to_show) <- c('M-STAT-T', 'T-NS-T', 'M-NS-T-G')
    
    cat('CRPS\n')
    print(matrix_to_show)
    cat('------\n')
    
  }
  
  # Q_0.9 CRPS
  
  if(T){
    
    Q_CRPS_A <- numeric(100)
    
    for(ii in 1:100){
      Q_CRPS_A[ii] <- quantile(CRPS_T_A[groups == levels(groups)[ii]],probs = 0.95)
    }

    Q_CRPS_B <- numeric(100)
    
    for(ii in 1:100){
      Q_CRPS_B[ii] <- quantile(CRPS_T_B[groups == levels(groups)[ii]],probs = 0.95)
    }

    Q_CRPS_C <- numeric(100)
    
    for(ii in 1:100){
      Q_CRPS_C[ii] <- quantile(CRPS_T_C[groups == levels(groups)[ii]],probs=0.95)
    }

    matrix_to_show <- t(round(cbind(c(mean(Q_CRPS_B), mean(Q_CRPS_A), mean(Q_CRPS_C)),
                                    c(sd(Q_CRPS_B), sd(Q_CRPS_A), sd(Q_CRPS_C))), 3))
    
    rownames(matrix_to_show) <- c('mean', 'se')
    colnames(matrix_to_show) <- c('M-STAT-T', 'T-NS-T', 'M-NS-T-G')
    
    cat('0.95 Quantile CRPS\n')
    print(matrix_to_show)
    cat('---------\n')
    
  }
  
  # D_n
  
  if(T){
    
    KS_A <- numeric(100)
    
    for(ii in 1:100){
      KS_A[ii] <- ks.test(z_std_A[groups == levels(groups)[ii]],y = pnorm)$statistic
    }
   
    KS_B <- numeric(100)
    
    for(ii in 1:100){
      KS_B[ii] <- ks.test(z_std_B[groups == levels(groups)[ii]],
                          y = pnorm)$statistic
    }
    
    KS_C <- numeric(100)
    
    for(ii in 1:100){
      KS_C[ii] <- ks.test(z_std_C[groups == levels(groups)[ii]],
                          y = pnorm)$statistic
    }
    
    matrix_to_show <- t(round(cbind(c(mean(KS_B), mean(KS_A), mean(KS_C)),
                                    c(sd(KS_B), sd(KS_A), sd(KS_C))),3))
    
    rownames(matrix_to_show) <- c('mean', 'se')
    colnames(matrix_to_show) <- c('M-STAT-T', 'T-NS-T', 'M-NS-T-G')
    
    cat('KS statistic \n')
    print(matrix_to_show)
    cat('--------\n')
    
  }
  
  # Coverage_prob
  
  if(T){
    
    Cov_prob_A <- numeric(100)
    
    for(ii in 1:100){
      
      upper_bound <- Pred_Model_T_A$trend[groups == levels(groups)[ii]] + 
        Pred_Model_T_A$mean[groups == levels(groups)[ii]] +
        qnorm(1 - 0.025) * Pred_Model_T_A$sd.pred[groups == levels(groups)[ii]]
      
      lower_bound <- Pred_Model_T_A$trend[groups == levels(groups)[ii]] + 
        Pred_Model_T_A$mean[groups == levels(groups)[ii]] -
        qnorm(1 - 0.025) * Pred_Model_T_A$sd.pred[groups == levels(groups)[ii]]
      
      Cov_prob_A[ii] <-  1 - length(which(z_values[groups == levels(groups)[ii]] < lower_bound |
                                            z_values[groups == levels(groups)[ii]] > upper_bound  )) /
        length(z_values[groups == levels(groups)[ii]])
      
    }

    Cov_prob_B <- numeric(100)
    
    for(ii in 1:100){
      
      upper_bound <- Pred_Model_T_B$trend[groups == levels(groups)[ii]] + 
        Pred_Model_T_B$mean[groups == levels(groups)[ii]] +
        qnorm(1 - 0.025) * Pred_Model_T_B$sd.pred[groups == levels(groups)[ii]]
      
      lower_bound <- Pred_Model_T_B$trend[groups == levels(groups)[ii]] + 
        Pred_Model_T_B$mean[groups == levels(groups)[ii]] -
        qnorm(1 - 0.025) * Pred_Model_T_B$sd.pred[groups == levels(groups)[ii]]
      
      Cov_prob_B[ii] <-  1 - length(which(z_values[groups == levels(groups)[ii]] < lower_bound |
                                            z_values[groups == levels(groups)[ii]] > upper_bound  )) /
        length(z_values[groups == levels(groups)[ii]])
      
      
    }
    
    Cov_prob_C <- numeric(100)
    
    for(ii in 1:100){
      
      upper_bound <- Pred_Model_T_C$trend[groups == levels(groups)[ii]] + 
        Pred_Model_T_C$mean[groups == levels(groups)[ii]] +
        qnorm(1 - 0.025) * Pred_Model_T_C$sd.pred[groups == levels(groups)[ii]]
      
      lower_bound <- Pred_Model_T_C$trend[groups == levels(groups)[ii]] + 
        Pred_Model_T_C$mean[groups == levels(groups)[ii]] -
        qnorm(1 - 0.025) * Pred_Model_T_C$sd.pred[groups == levels(groups)[ii]]
      
      Cov_prob_C[ii] <-  1 - length(which(z_values[groups == levels(groups)[ii]] < lower_bound |
                                            z_values[groups == levels(groups)[ii]] > upper_bound  )) /
        length(z_values[groups == levels(groups)[ii]])
      
    }
    
    matrix_to_show <- t(round(cbind(c(mean(Cov_prob_B),mean(Cov_prob_A),mean(Cov_prob_C)),
                                    c(sd(Cov_prob_B),sd(Cov_prob_A),sd(Cov_prob_C))),3))
    
    rownames(matrix_to_show) <- c('mean', 'se')
    colnames(matrix_to_show) <- c('M-STAT-T', 'T-NS-T', 'M-NS-T-G')
    
    cat('coverage probability\n')
    print(matrix_to_show)
    cat('--------------\n')
    
  }
  
  # log-score
  
  if(T){
    
    LS_A <- numeric(100)
    
    for(ii in 1:100){
      LS_A[ii] <- mean(Logscore_T_A[groups == levels(groups)[ii]])
    }

    LS_B <- numeric(100)
    
    for(ii in 1:100){
      LS_B[ii] <-  mean(Logscore_T_B[groups == levels(groups)[ii]])
    }
    
    LS_C <- numeric(100)
    
    for(ii in 1:100){
      LS_C[ii] <-  mean(Logscore_T_C[groups == levels(groups)[ii]])
    }
    
    matrix_to_show <- t(round(cbind(c(mean(LS_B), mean(LS_A), mean(LS_C)),
                                    c(sd(LS_B), sd(LS_A), sd(LS_C))), 3))
    
    rownames(matrix_to_show) <- c('mean', 'se')
    colnames(matrix_to_show) <- c('M-STAT-T', 'T-NS-T', 'M-NS-T-G')
    
    cat('logscore\n')
    print(matrix_to_show)
    cat('--------\n')
    
  }
  
  # logliks
  
  vector_to_show <- c(getLoglik(Model_T_B),getLoglik(Model_T_A),getLoglik(Model_T_C))
  names(vector_to_show) <- c('M-STAT-T', 'T-NS-T', 'M-NS-T-G')
  
  cat('Logliks\n')
  print(vector_to_show)
  cat('--------\n')
  
  # times
  
  vector_to_show <- round(c(Time_T_B[3]/60,time_T_A[3]/60,Time_T_C[3]/60),2)
  names(vector_to_show) <- c('M-STAT-T', 'T-NS-T', 'M-NS-T-G')
  
  cat('system times\n')
  print(vector_to_show)
  cat('--------\n')
  
}
