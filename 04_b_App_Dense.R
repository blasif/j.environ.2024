source('00_system.R')

load('Datasets/app_dataset.RData')

# Initial Model (Model A)

if(T){
  
  list_formulas <- list("mean" = as.formula(" ~ 1 + lati + long + wind + elevation + cloud_c"),
                        "std.dev" = as.formula("  ~ 1  + log_cloud + log_elevation + log_wind + log_long + log_lati"),
                        "scale" = as.formula("  ~ 1  + log_cloud + log_elevation + log_wind + log_long + log_lati"),
                        "aniso" = as.formula("  ~ 1  + log_cloud + log_elevation + log_wind + log_long + log_lati"),
                        "tilt" = as.formula("  ~ 1  + log_cloud + log_elevation + log_wind + log_long + log_lati"),
                        "smooth" = as.formula("  ~ 1  + log_cloud + log_elevation + log_wind + log_long + log_lati"),
                        "nugget" = -Inf)
  
  test_coco <- coco(type = 'dense',
                    data = all_dfs[[2]],
                    locs = as.matrix(all_dfs[[2]][, 1:2]),
                    z = all_dfs[[2]]$prec,
                    model.list = list_formulas,
                    info = list('lambda' = 0.025 * sqrt(8), 
                                'smooth_limits' = c(0.5, 2.0)
                                )
                    )
  
  boundaries_A <- getBoundariesV2(coco.object = test_coco,
                                  mean.limits = c(-1.5, 0, 1.5),
                                  std.dev.limits = c(-2.5, 0, 2.5),
                                  scale.limits = c(-2.5, 0, 2.5),
                                  aniso.limits =  c(-2, 0, 2),
                                  tilt.limits =  c(-2, 0, 2),
                                  smooth.limits = c(-2, 0, 1.5),
                                  nugget.limits = c(-5, -3, 1))
  
  # Tailored boundaries
  if(T){
    first_var <- which(names(boundaries_A$theta_init) == "std.dev.limits")[1]
    first_range <- which(names(boundaries_A$theta_init) == "scale.limits")[1]
    first_smooth <- which(names(boundaries_A$theta_init) == "smooth.limits")[1]
    
    boundaries_A$theta_upper[c(first_var, first_range)] <- c(5, 5)
    boundaries_A$theta_lower[c(first_var, first_range)] <- c(-5, -1.75)
    
    boundaries_A$theta_init[first_range] <- log(sd(c(dist(test_coco@locs))))
    
    boundaries_A$theta_upper[first_smooth] <- 2
    boundaries_A$theta_lower[first_smooth] <- -3
    boundaries_A$theta_init[first_smooth] <- -3
    
    boundaries_A$theta_init[1] <- mean(test_coco@z)
    boundaries_A$theta_upper[1] <- boundaries_A$theta_init[1] + 3
    boundaries_A$theta_lower[1] <- boundaries_A$theta_init[1] - 3
    
  }
  
  
  
  Time_A <- system.time({Model_A <- cocoOptim(coco.object = test_coco,
                                              ncores = ncores,
                                              boundaries = boundaries_A,
                                              optim.type = "mle")})
  
  plot(Model_A)
  
  save(Time_A,Model_A, file = 'Rdata/Model_A.RData')  

}

# Model B - A (better) specification of the trend

if(T){
  
  list_formulas <- list("mean" = as.formula("  ~ 1 + wind + new_wind + BIO04 + BIO15 + cloud_c + elevation + lati + long"),
                        "std.dev" = as.formula("  ~ 1 + log_cloud + log_wind + log_elevation + log_lati + log_long"),
                        "scale" = as.formula("  ~ 1 + log_cloud + log_wind + log_elevation + log_lati + log_long"),
                        "aniso" = as.formula("  ~ 1 + log_cloud + log_wind + log_elevation + log_lati + log_long"),
                        "tilt" = as.formula("  ~ 1 + log_cloud + log_wind + log_elevation + log_lati + log_long"),
                        "smooth" = as.formula(" ~ 1 + log_cloud + log_wind + log_elevation + log_lati + log_long"),
                        "nugget" = -Inf)
  
  test_coco_classic <- coco(type = 'dense',
                            data = all_dfs[[2]],
                            locs = as.matrix(all_dfs[[2]][, 1:2]),
                            z = all_dfs[[2]]$prec,
                            model.list = list_formulas,
                            info = list('lambda' = 0.025 * sqrt(8), # before was 0.075
                                        'smooth_limits' = c(0.5, 2.0)
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
                                              ncores = ncores,
                                              boundaries = boundaries_B,
                                              optim.type = 'mle')})
  
  HESS_B <- getHessian(Model_B, ncores = ncores)
  
  aaa <- getCIs(coco.object = Model_B, inv.hess = solve(HESS_B), alpha = 0.2)
  
  which(aaa[,1] < 0 & aaa[,2] > 0)
  
  save(Time_B, Model_B, HESS_B, file = 'RData/Model_B.RData')
  
}

# Model C - (informally) Reduced

if(F){
  
  list_formulas <- list("mean" = as.formula("  ~ 1 + wind + new_wind + BIO04 + BIO15 + cloud_c + elevation + lati + long"),
                        "std.dev" = as.formula("  ~ 1 + log_cloud + log_wind + log_elevation"),
                        "scale" = as.formula("  ~ 1 + log_wind + log_elevation"),
                        "aniso" = as.formula("  ~ 1 + log_lati + log_long"),
                        "tilt" = as.formula("  ~ 1 + log_cloud"),
                        "smooth" = as.formula(" ~ 1 + log_wind + log_long"),
                        "nugget" = -Inf)
  
  test_coco_red <- coco(type = 'dense',
                        data = all_dfs[[2]],
                        locs = as.matrix(all_dfs[[2]][, 1:2]),
                        z = all_dfs[[2]]$prec,
                        model.list = list_formulas,
                        info = list('lambda' = 0.025 * sqrt(8),
                                    'smooth_limits' = c(0.5, 2.0)
                                    )
                        )
  
  # Tailored boundaries
  if(T){
    boundaries_C <- getBoundariesV2(coco.object = test_coco_red,
                                    mean.limits = c(-1.5, 0, 1.5),
                                    std.dev.limits = c(-2.5, 0, 2.5),
                                    scale.limits = c(-2.5, 0, 2.5),
                                    aniso.limits =  c(-2, 0, 2),
                                    tilt.limits =  c(-2, 0, 2),
                                    smooth.limits = c(-2, 0, 1.5),
                                    nugget.limits = c(-5, -3, 1))
    
    first_var <- which(names(boundaries_C$theta_init) == "std.dev.limits")[1]
    first_range <- which(names(boundaries_C$theta_init) == "scale.limits")[1]
    first_smooth <- which(names(boundaries_C$theta_init) == "smooth.limits")[1]
    
    boundaries_C$theta_upper[c(first_var, first_range)] <- c(5, 5)
    boundaries_C$theta_lower[c(first_var, first_range)] <- c(-5, -1.75)
    
    boundaries_C$theta_init[first_range] <- log(sd(c(dist(test_coco_red@locs))))
    
    boundaries_C$theta_upper[first_smooth] <- 1.5
    boundaries_C$theta_lower[first_smooth] <- -3
    boundaries_C$theta_init[first_smooth] <- -3
    
    boundaries_C$theta_init[1] <- mean(test_coco_red@z)
    boundaries_C$theta_upper[1] <- boundaries_C$theta_init[1] + 3
    boundaries_C$theta_lower[1] <- boundaries_C$theta_init[1] - 1  
    
  }
  
  time_C <- system.time({Model_C <- cocoOptim(coco.object = test_coco_red,
                                              ncores = ncores,
                                              boundaries = boundaries_C,
                                              optim.type = 'mle')})
  
  HESS_C <- getHessian(Model_C, ncores = ncores)
  
  aaa <- getCIs(Model_C, inv.hess = solve(HESS_C), alpha = 0.2)
  which(aaa[,1] < 0 & aaa[,2] > 0)
  
  save(Model_C, time_C, HESS_C, file = 'RData/Model_C.RData')
  
}

# Model D - Classic model

if(T){
  
  list_formulas <- list("mean" = as.formula("   ~ 1 + wind + new_wind + BIO04 + BIO15 + cloud_c + elevation + lati + long"),
                        "std.dev" = as.formula("  ~ 1 "),
                        "scale" = as.formula("  ~ 1 "),
                        "aniso" = as.formula("  ~ 1 "),
                        "tilt" = as.formula("  ~ 1 "),
                        "smooth" = as.formula(" ~ 1 "),
                        "nugget" = -Inf)
  
  test_coco_classic <- coco(type = 'dense',
                            data = all_dfs[[2]],
                            locs = as.matrix(all_dfs[[2]][, 1:2]),
                            z = all_dfs[[2]]$prec,
                            model.list = list_formulas,
                            info = list('smooth_limits' = c(0.5, 2.0)
                                        ))
  
  # Tailored boundaries
  if(T){
    boundaries_D <- getBoundariesV2(coco.object = test_coco_classic,
                                    mean.limits = c(-3, 0, 3),
                                    std.dev.limits = c(-2.5, 0, 2.5),
                                    scale.limits = c(-2.5, 0, 2.5),
                                    aniso.limits =  c(-2, 0, 2),
                                    tilt.limits =  c(-2, 0, 2),
                                    smooth.limits = c(-3, 0, 1.5),
                                    nugget.limits = c(-5, -3, 1))
    
    boundaries_D$theta_upper[1] <- 6
    
    first_var <- which(names(boundaries_D$theta_init) == "std.dev.limits")[1]
    first_range <- which(names(boundaries_D$theta_init) == "scale.limits")[1]
    first_smooth <- which(names(boundaries_D$theta_init) == "smooth.limits")[1]
    
    boundaries_D$theta_init[1] <- c(mean(test_coco_classic@z))
    
    boundaries_D$theta_upper[c(first_var, first_range)] <- c(4, 5)
    boundaries_D$theta_lower[c(first_var, first_range)] <- c(-5, -3)
    
    boundaries_D$theta_init[first_range] <- log(sd(c(dist(test_coco_classic@locs))))
    
    boundaries_D$theta_upper[first_smooth] <- 2
    
  }
  
  Time_D <- system.time({Model_D <- cocoOptim(coco.object = test_coco_classic,
                                              ncores = ncores,
                                              boundaries = boundaries_D,
                                              optim.type = 'mle')})
  
  getLoglik(Model_D)
  getBIC(Model_D)

  save(Time_D, Model_D, file = 'RData/Model_D.RData')
  
}

# Model E - Global smoothness

if(T){
  
  list_formulas <- list("mean" = as.formula("  ~ 1 + wind + new_wind + BIO04 + BIO15 + cloud_c + elevation + lati + long"),
                        "std.dev" = as.formula("  ~ 1 + log_cloud + log_wind + log_elevation"),
                        "scale" = as.formula("  ~ 1 + log_wind + log_elevation"),
                        "aniso" = as.formula("  ~ 1 + log_lati + log_long"),
                        "tilt" = as.formula("  ~ 1 + log_cloud"),
                        "smooth" = as.formula(" ~ 1"),
                        "nugget" = -Inf)
  
  test_coco_red <- coco(type = 'dense',
                        data = all_dfs[[2]],
                        locs = as.matrix(all_dfs[[2]][, 1:2]),
                        z = all_dfs[[2]]$prec,
                        model.list = list_formulas,
                        info = list('lambda' = 0.025 * sqrt(8),
                                    'smooth_limits' = c(0.5, 2.0)
                        ))
  
  # Tailored boundaries
  if(T){
    boundaries_E <- getBoundariesV2(coco.object = test_coco_red,
                                    mean.limits = c(-1.5, 0, 1.5),
                                    std.dev.limits = c(-2.5, 0, 2.5),
                                    scale.limits = c(-2.5, 0, 2.5),
                                    aniso.limits =  c(-2, 0, 2),
                                    tilt.limits =  c(-2, 0, 2),
                                    smooth.limits = c(-2, 0, 1.5),
                                    nugget.limits = c(-5, -3, 1))
    
    first_var <- which(names(boundaries_E$theta_init) == "std.dev.limits")[1]
    first_range <- which(names(boundaries_E$theta_init) == "scale.limits")[1]
    first_smooth <- which(names(boundaries_E$theta_init) == "smooth.limits")[1]
    
    boundaries_E$theta_upper[c(first_var, first_range)] <- c(5, 5)
    boundaries_E$theta_lower[c(first_var, first_range)] <- c(-5, -1.75)
    
    boundaries_E$theta_init[first_range] <- log(sd(c(dist(test_coco_red@locs))))
    
    boundaries_E$theta_upper[first_smooth] <- 1.5
    boundaries_E$theta_lower[first_smooth] <- -3
    boundaries_E$theta_init[first_smooth] <- -3
    
    boundaries_E$theta_init[1] <- mean(test_coco_red@z)
    boundaries_E$theta_upper[1] <- boundaries_E$theta_init[1] + 3
    boundaries_E$theta_lower[1] <- boundaries_E$theta_init[1] - 3  
    
  }
  
  time_E <- system.time({Model_E <- cocoOptim(coco.object = test_coco_red,
                                              ncores = ncores,
                                              boundaries = boundaries_E,
                                              optim.type = 'mle')})
  
  save(Model_E, time_E, file = 'RData/Model_E.RData')
  
}

# Predictions

if(T){
  
  newdataset <- all_dfs[[1]]
  
  newlocs <- all_dfs[[1]][ ,1:2]
  
  z_values <- all_dfs[[1]]$prec
  
  # Create holdouts
  if(T){
    set.seed(100621)
    hetero_holdouts <- kmeans(as.data.frame(scale(newdataset[,c(1,2,4:9)])),centers = 100,iter.max = 100)
    groups <- as.factor(hetero_holdouts$cluster)
    quilt.plot(newlocs, hetero_holdouts$cluster, nx = 150, ny = 150)
  }
  
  Pred_B <- cocoPredict(coco.object = Model_B, 
                         newdataset = newdataset, 
                         newlocs = newlocs,
                         type = 'pred')
  
  Pred_C <- cocons::cocoPredict(coco.object = Model_C, 
                                newdataset = newdataset, 
                                newlocs = newlocs,
                                type = 'pred')
  
  Pred_D <- cocons::cocoPredict(coco.object = Model_D, 
                                newdataset = newdataset, 
                                newlocs = newlocs,
                                type = 'pred')
  
  Pred_E <- cocons::cocoPredict(coco.object = Model_E, 
                                newdataset = newdataset, 
                                newlocs = newlocs,
                                type = 'pred')
  
  CRPS_B <- getCRPS(z_values,mean.pred = Pred_B$trend + Pred_B$mean,sd.pred = Pred_B$sd.pred)
  CRPS_C <- getCRPS(z_values,mean.pred = Pred_C$trend + Pred_C$mean,sd.pred = Pred_C$sd.pred)
  CRPS_D <- getCRPS(z_values,mean.pred = Pred_D$trend + Pred_D$mean,sd.pred = Pred_D$sd.pred)
  CRPS_E <- getCRPS(z_values,mean.pred = Pred_E$trend + Pred_E$mean,sd.pred = Pred_E$sd.pred)
  
  Logscore_B <- getLogRank(z_values,mean.pred = Pred_B$trend + Pred_B$mean,sd.pred = Pred_B$sd.pred)
  Logscore_C <- getLogRank(z_values,mean.pred = Pred_C$trend + Pred_C$mean,sd.pred = Pred_C$sd.pred)
  Logscore_D <- getLogRank(z_values,mean.pred = Pred_D$trend + Pred_D$mean,sd.pred = Pred_D$sd.pred)
  Logscore_E <- getLogRank(z_values,mean.pred = Pred_E$trend + Pred_E$mean,sd.pred = Pred_E$sd.pred)
  
  z_std_B <- (Pred_B$trend + Pred_B$mean - z_values) / Pred_B$sd.pred
  z_std_C <- (Pred_C$trend + Pred_C$mean - z_values) / Pred_C$sd.pred
  z_std_D <- (Pred_D$trend + Pred_D$mean - z_values) / Pred_D$sd.pred
  z_std_E <- (Pred_E$trend + Pred_E$mean - z_values) / Pred_E$sd.pred
  
  save(Pred_B, Pred_C, Pred_D, Pred_E, 
       CRPS_B, CRPS_C, CRPS_D, CRPS_E, 
       Logscore_B, Logscore_C, Logscore_D, Logscore_E,
       z_std_B, z_std_C, z_std_D, z_std_E,
       file = 'RData/Pred_dense.RData')
  
}

# Metrics

load('RData/Pred_dense.RData')

if(T){
  
  # RMSPE
  
  if(T){
    
    RMSPE_C <- numeric(100)
    
    for(ii in 1:100){
      RMSPE_C[ii] <- sqrt(mean((Pred_C$mean[groups == levels(groups)[ii]] + Pred_C$trend[groups == levels(groups)[ii]] - z_values[groups == levels(groups)[ii]])^2))
    }
    
    mean(RMSPE_C)
    sd(RMSPE_C)
    
    RMSPE_D <- numeric(100)
    
    for(ii in 1:100){
      RMSPE_D[ii] <- sqrt(mean((Pred_D$mean[groups == levels(groups)[ii]] + Pred_D$trend[groups == levels(groups)[ii]] - z_values[groups == levels(groups)[ii]])^2))
    }
    
    mean(RMSPE_D)
    sd(RMSPE_D)
    
    RMSPE_E <- numeric(100)
    
    for(ii in 1:100){
      RMSPE_E[ii] <- sqrt(mean((Pred_E$mean[groups == levels(groups)[ii]] + 
                                  Pred_E$trend[groups == levels(groups)[ii]] - 
                                  z_values[groups == levels(groups)[ii]])^2))
    }
    
    mean(RMSPE_E)
    sd(RMSPE_E)
    
    matrix_to_show <- t(round(cbind(c(mean(RMSPE_D),mean(RMSPE_C),mean(RMSPE_E)),
                                    c(sd(RMSPE_D),sd(RMSPE_C),sd(RMSPE_E))),3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS','M-NS-G')
    
    cat('RMSPE\n')
    print(matrix_to_show)
    cat('------\n')
    
  }
  
  # CRPS
  
  if(T){
    
    CRPS_C_holdout <- numeric(100)
    
    for(ii in 1:100){
      CRPS_C_holdout[ii] <- mean(CRPS_C[groups == levels(groups)[ii]])
    }
    
    mean(CRPS_C_holdout)
    sd(CRPS_C_holdout)
    
    CRPS_D_holdout <- numeric(100)
    
    for(ii in 1:100){
      CRPS_D_holdout[ii] <- mean(CRPS_D[groups == levels(groups)[ii]])
    }
    
    mean(CRPS_D_holdout)
    sd(CRPS_D_holdout)
    
    CRPS_E_holdout <- numeric(100)
    
    for(ii in 1:100){
      CRPS_E_holdout[ii] <- mean(CRPS_E[groups == levels(groups)[ii]])
    }
    
    mean(CRPS_E_holdout)
    sd(CRPS_E_holdout)
    
    matrix_to_show <- t(round(cbind(c(mean(CRPS_D_holdout),mean(CRPS_C_holdout),mean(CRPS_E_holdout)),
                                    c(sd(CRPS_D_holdout),sd(CRPS_C_holdout),sd(CRPS_E_holdout))),3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS','M-NS-G')
    
    cat('CRPS\n')
    print(matrix_to_show)
    cat('------\n')
    
  }
  
  # Q_0.9 CRPS
  
  if(T){
    
    Q_CRPS_C <- numeric(100)
    
    for(ii in 1:100){
      Q_CRPS_C[ii] <- quantile(CRPS_C[groups == levels(groups)[ii]],probs = 0.95)
    }
    
    mean(Q_CRPS_C)
    sd(Q_CRPS_C)
    
    Q_CRPS_D <- numeric(100)
    
    for(ii in 1:100){
      Q_CRPS_D[ii] <- quantile(CRPS_D[groups == levels(groups)[ii]],probs = 0.95)
    }
    
    mean(Q_CRPS_D)
    sd(Q_CRPS_D)
    
    Q_CRPS_E <- numeric(100)
    
    for(ii in 1:100){
      Q_CRPS_E[ii] <- quantile(CRPS_E[groups == levels(groups)[ii]],probs=0.95)
    }
    
    mean(Q_CRPS_E)
    sd(Q_CRPS_E)
    
    matrix_to_show <- t(round(cbind(c(mean(Q_CRPS_D),mean(Q_CRPS_C),mean(Q_CRPS_E)),
                                    c(sd(Q_CRPS_D),sd(Q_CRPS_C),sd(Q_CRPS_E))),3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS','M-NS-G')
    
    cat('0.95 Quantile CRPS\n')
    print(matrix_to_show)
    cat('---------\n')
    
  }
  
  # D_n
  
  if(T){
    
    KS_C <- numeric(100)
    
    for(ii in 1:100){
      KS_C[ii] <- ks.test(z_std_C[groups == levels(groups)[ii]],y = pnorm)$statistic
    }
    
    mean(KS_C)
    sd(KS_C)
    
    KS_D <- numeric(100)
    
    for(ii in 1:100){
      KS_D[ii] <- ks.test(z_std_D[groups == levels(groups)[ii]],
                          y = pnorm)$statistic
    }
    
    mean(KS_D)
    sd(KS_D)  
    
    KS_E <- numeric(100)
    
    for(ii in 1:100){
      KS_E[ii] <- ks.test(z_std_E[groups == levels(groups)[ii]],
                          y = pnorm)$statistic
    }
    
    mean(KS_E)
    sd(KS_E)  
    
    matrix_to_show <- t(round(cbind(c(mean(KS_D),mean(KS_C),mean(KS_E)),
                                    c(sd(KS_D),sd(KS_C),sd(KS_E))),3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS','M-NS-G')
    
    cat('KS statistic \n')
    print(matrix_to_show)
    cat('--------\n')
    
  }
  
  # Coverage_prob
  
  if(T){
    
    Cov_prob_C <- numeric(100)
    
    for(ii in 1:100){
      
      upper_bound <- Pred_C$trend[groups == levels(groups)[ii]] + 
        Pred_C$mean[groups == levels(groups)[ii]] +
        qnorm(1 - 0.025) * Pred_C$sd.pred[groups == levels(groups)[ii]]
      
      lower_bound <- Pred_C$trend[groups == levels(groups)[ii]] + 
        Pred_C$mean[groups == levels(groups)[ii]] -
        qnorm(1 - 0.025) * Pred_C$sd.pred[groups == levels(groups)[ii]]
      
      Cov_prob_C[ii] <-  1 - length(which(z_values[groups == levels(groups)[ii]] < lower_bound |
                                            z_values[groups == levels(groups)[ii]] > upper_bound  )) /
        length(z_values[groups == levels(groups)[ii]])
      
    }
    
    mean(Cov_prob_C)
    sd(Cov_prob_C)
    
    Cov_prob_D <- numeric(100)
    
    for(ii in 1:100){
      
      upper_bound <- Pred_D$trend[groups == levels(groups)[ii]] + 
        Pred_D$mean[groups == levels(groups)[ii]] +
        qnorm(1 - 0.025) * Pred_D$sd.pred[groups == levels(groups)[ii]]
      
      lower_bound <- Pred_D$trend[groups == levels(groups)[ii]] + 
        Pred_D$mean[groups == levels(groups)[ii]] -
        qnorm(1 - 0.025) * Pred_D$sd.pred[groups == levels(groups)[ii]]
      
      Cov_prob_D[ii] <-  1 - length(which(z_values[groups == levels(groups)[ii]] < lower_bound |
                                            z_values[groups == levels(groups)[ii]] > upper_bound  )) /
        length(z_values[groups == levels(groups)[ii]])
      
      
    }
    
    Cov_prob_E <- numeric(100)
    
    for(ii in 1:100){
      
      upper_bound <- Pred_E$trend[groups == levels(groups)[ii]] + 
        Pred_E$mean[groups == levels(groups)[ii]] +
        qnorm(1 - 0.025) * Pred_E$sd.pred[groups == levels(groups)[ii]]
      
      lower_bound <- Pred_E$trend[groups == levels(groups)[ii]] + 
        Pred_E$mean[groups == levels(groups)[ii]] -
        qnorm(1 - 0.025) * Pred_E$sd.pred[groups == levels(groups)[ii]]
      
      Cov_prob_E[ii] <-  1 - length(which(z_values[groups == levels(groups)[ii]] < lower_bound |
                                            z_values[groups == levels(groups)[ii]] > upper_bound  )) /
        length(z_values[groups == levels(groups)[ii]])
      
    }
    
    mean(Cov_prob_E)
    sd(Cov_prob_E)
    
    matrix_to_show <- t(round(cbind(c(mean(Cov_prob_D),mean(Cov_prob_C),mean(Cov_prob_E)),
                                    c(sd(Cov_prob_D),sd(Cov_prob_C),sd(Cov_prob_E))),3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS','M-NS-G')
    
    
    cat('coverage probability\n')
    print(matrix_to_show)
    cat('--------------\n')
    
  }
  
  # log-score
  
  if(T){
    
    LS_C <- numeric(100)
    
    for(ii in 1:100){
      LS_C[ii] <- mean(Logscore_C[groups == levels(groups)[ii]])
    }
    
    mean(LS_C)
    sd(LS_C)
    
    LS_D <- numeric(100)
    
    for(ii in 1:100){
      LS_D[ii] <-  mean(Logscore_D[groups == levels(groups)[ii]])
    }
    
    LS_E <- numeric(100)
    
    for(ii in 1:100){
      LS_E[ii] <-  mean(Logscore_E[groups == levels(groups)[ii]])
    }
    
    mean(LS_E)
    sd(LS_E)
    
    matrix_to_show <- t(round(cbind(c(mean(LS_D),mean(LS_C),mean(LS_E)),
                                    c(sd(LS_D),sd(LS_C),sd(LS_E))),3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS','M-NS-G')
    
    cat('logscore\n')
    print(matrix_to_show)
    cat('--------\n')
    
  }
  
  # logliks
  
  vector_to_show <- c(getLoglik(Model_D),getLoglik(Model_C),getLoglik(Model_E))
  names(vector_to_show) <- c('M-STAT','M-NS','M-NS-G')
  
  cat('Logliks\n')
  print(vector_to_show)
  cat('--------\n')
  
  # times
  
  vector_to_show <- round(c(Time_D[3]/60,time_C[3]/60,time_E[3]/60),2)
  names(vector_to_show) <- c('M-STAT','M-NS','M-NS-G')
  
  cat('system times\n')
  print(vector_to_show)
  cat('--------\n')
  
}
