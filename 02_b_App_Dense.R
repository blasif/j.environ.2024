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
                            z = all_dfs[[2]]$prec,
                            model.list = list_formulas,
                            info = list('smooth.limits' = c(0.5, 2.0)
                                        )
                            )
  
  # Tailored boundaries
  if(T){
    
    boundaries_STAT <- getBoundariesV4(test_coco_classic)
    
    boundaries_STAT$theta_upper['std.dev.limits'] <- 4
    boundaries_STAT$theta_lower['std.dev.limits'] <- -4
    
    boundaries_STAT$theta_upper['scale.limits'] <- 4
    boundaries_STAT$theta_lower['scale.limits'] <- -4
    
  }
  
  Time_STAT <- system.time({Model_STAT <- cocoOptim(coco.object = test_coco_classic,
                                              ncores = 7,
                                              boundaries = boundaries_STAT,
                                              optim.type = 'ml')})
  
  HESS_STAT <- getHessian(Model_STAT, ncores = ncores)

  save(Time_STAT, Model_STAT, HESS_STAT, file = 'RData/Model_STAT.RData')
  
}

# Predictions

if(T){
  
  # Create holdouts
  if(T){
    
    newdataset <- all_dfs[[1]]
    newlocs <- all_dfs[[1]][ , c("long","lati")]
    z_values <- all_dfs[[1]]$prec
    
    set.seed(100621)
    hetero_holdouts <- kmeans(as.data.frame(scale(newdataset[,c(1,2,4:9)])),centers = 100,iter.max = 100)
    groups <- as.factor(hetero_holdouts$cluster)
    quilt.plot(newlocs, hetero_holdouts$cluster, nx = 150, ny = 150)
    
    sample_to_tune_hyperparameters <- sample(1:100, 30)
    
    newdataset_final <- newdataset[!(groups %in% sample_to_tune_hyperparameters),]
    newlocs_final <- as.matrix(newdataset[!(groups %in% sample_to_tune_hyperparameters), c("long","lati")])
    z_values_final <- all_dfs[[1]]$prec[!(groups %in% sample_to_tune_hyperparameters)]
    
    final_groups <- c(1:100)[!(c(1:100) %in% sample_to_tune_hyperparameters)]
    
    obs_groups_final <- groups[groups %in% final_groups]
    
  }  
  
  load('RData/Model_NS.RData')

  Pred_NS <- cocoPredict(coco.object = Model_NS, 
                         newdataset = newdataset_final, 
                         newlocs = newlocs_final,
                         type = 'pred')
  
  Pred_STAT <- cocons::cocoPredict(coco.object = Model_STAT, 
                                newdataset = newdataset_final, 
                                newlocs = newlocs_final,
                                type = 'pred')
  
  CRPS_NS <- getCRPS(z_values_final, mean.pred = Pred_NS$systematic + Pred_NS$stochastic, sd.pred = Pred_NS$sd.pred)
  CRPS_STAT <- getCRPS(z_values_final, mean.pred = Pred_STAT$systematic + Pred_STAT$stochastic, sd.pred = Pred_STAT$sd.pred)
  
  Logscore_NS <- getLogScore(z_values_final, mean.pred = Pred_NS$systematic + Pred_NS$stochastic, sd.pred = Pred_NS$sd.pred)
  Logscore_STAT <- getLogScore(z_values_final, mean.pred = Pred_STAT$systematic + Pred_STAT$stochastic, sd.pred = Pred_STAT$sd.pred)

  z_std_NS <- (Pred_NS$systematic + Pred_NS$stochastic - z_values_final) / Pred_NS$sd.pred
  z_std_STAT <- (Pred_STAT$systematic + Pred_STAT$stochastic - z_values_final) / Pred_STAT$sd.pred

  save(Pred_NS, Pred_STAT,
       CRPS_NS, CRPS_STAT, 
       Logscore_NS, Logscore_STAT,
       z_std_NS, z_std_STAT,
       file = 'RData/Pred_dense.RData')
  
}

# Metrics

load('RData/Pred_dense.RData')

if(T){
  
  # RMSPE
  
  if(T){
    
    RMSPE_NS <- numeric(70)
    
    for(ii in 1:70){
      index_groups <- obs_groups_final == final_groups[ii]
      RMSPE_NS[ii] <- sqrt(mean((Pred_NS$systematic[index_groups] + Pred_NS$stochastic[index_groups] - z_values_final[index_groups])^2))
    }

    RMSPE_STAT <- numeric(70)
    
    for(ii in 1:70){
      index_groups <- obs_groups_final == final_groups[ii]
      RMSPE_STAT[ii] <- sqrt(mean((Pred_STAT$systematic[index_groups] + Pred_STAT$stochastic[index_groups] - z_values_final[index_groups])^2))
    }

    matrix_to_show <- t(round(cbind(c(mean(RMSPE_STAT),mean(RMSPE_NS)),
                                    c(sd(RMSPE_STAT),sd(RMSPE_NS))),3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS')
    
    cat('RMSPE\n')
    print(matrix_to_show)
    cat('------\n')
    
  }
  
  # CRPS
  
  if(T){
    
    CRPS_NS_holdout <- numeric(70)
    
    for(ii in 1:70){
      CRPS_NS_holdout[ii] <- mean(CRPS_NS[obs_groups_final == final_groups[ii]])
    }

    CRPS_STAT_holdout <- numeric(70)
    
    for(ii in 1:70){
      CRPS_STAT_holdout[ii] <- mean(CRPS_STAT[obs_groups_final == final_groups[ii]])
    }

    matrix_to_show <- t(round(cbind(c(mean(CRPS_STAT_holdout),mean(CRPS_NS_holdout)),
                                    c(sd(CRPS_STAT_holdout),sd(CRPS_NS_holdout))), 3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS')
    
    cat('CRPS\n')
    print(matrix_to_show)
    cat('------\n')
    
  }
  
  # Q_0.9 CRPS
  
  if(T){
    
    Q_CRPS_NS <- numeric(70)
    
    for(ii in 1:70){
      Q_CRPS_NS[ii] <- quantile(CRPS_NS[obs_groups_final == final_groups[ii]],probs = 0.95)
    }

    Q_CRPS_STAT <- numeric(70)
    
    for(ii in 1:70){
      Q_CRPS_STAT[ii] <- quantile(CRPS_STAT[obs_groups_final == final_groups[ii]],probs = 0.95)
    }

    matrix_to_show <- t(round(cbind(c(mean(Q_CRPS_STAT),mean(Q_CRPS_NS)),
                                    c(sd(Q_CRPS_STAT),sd(Q_CRPS_NS))), 3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS')
    
    cat('0.95 Quantile CRPS\n')
    print(matrix_to_show)
    cat('---------\n')
    
  }
  
  # D_n
  
  if(T){
    
    KS_NS <- numeric(70)
    
    for(ii in 1:70){
      KS_NS[ii] <- ks.test(z_std_NS[obs_groups_final == final_groups[ii]],y = pnorm)$statistic
    }

    KS_STAT <- numeric(70)
    
    for(ii in 1:70){
      KS_STAT[ii] <- ks.test(z_std_STAT[obs_groups_final == final_groups[ii]],
                          y = pnorm)$statistic
    }

    matrix_to_show <- t(round(cbind(c(mean(KS_STAT),mean(KS_NS)),
                                    c(sd(KS_STAT),sd(KS_NS))),3))
    
    rownames(matrix_to_show) <- c('mean', 'se')
    colnames(matrix_to_show) <- c('M-STAT', 'M-NS')
    
    cat('KS statistic \n')
    print(matrix_to_show)
    cat('--------\n')
    
  }
  
  # Coverage_prob
  
  if(T){
    
    Cov_prob_NS <- numeric(70)
    
    for(ii in 1:70){
      
      index_groups <- obs_groups_final == final_groups[ii]
      
      upper_bound <- Pred_NS$systematic[index_groups] + 
        Pred_NS$stochastic[index_groups] +
        qnorm(1 - 0.025) * Pred_NS$sd.pred[index_groups]
      
      lower_bound <- Pred_NS$systematic[index_groups] + 
        Pred_NS$stochastic[index_groups] -
        qnorm(1 - 0.025) * Pred_NS$sd.pred[index_groups]
      
      Cov_prob_NS[ii] <-  1 - length(which(z_values_final[index_groups] < lower_bound |
                                            z_values_final[index_groups] > upper_bound  )) /
        length(z_values_final[index_groups])
      
    }

    Cov_prob_STAT <- numeric(70)
    
    for(ii in 1:70){
      
      index_groups <- obs_groups_final == final_groups[ii]
      
      upper_bound <- Pred_STAT$systematic[index_groups] + 
        Pred_STAT$stochastic[index_groups] +
        qnorm(1 - 0.025) * Pred_STAT$sd.pred[index_groups]
      
      lower_bound <- Pred_STAT$systematic[index_groups] + 
        Pred_STAT$stochastic[index_groups] -
        qnorm(1 - 0.025) * Pred_STAT$sd.pred[index_groups]
      
      Cov_prob_STAT[ii] <-  1 - length(which(z_values_final[index_groups] < lower_bound |
                                            z_values_final[index_groups] > upper_bound  )) /
        length(z_values_final[index_groups])
      
      
    }

    matrix_to_show <- t(round(cbind(c(mean(Cov_prob_STAT),mean(Cov_prob_NS)),
                                    c(sd(Cov_prob_STAT),sd(Cov_prob_NS))),3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS')
    
    
    cat('coverage probability\n')
    print(matrix_to_show)
    cat('--------------\n')
    
  }
  
  # log-score
  
  if(T){
    
    LS_NS <- numeric(70)
    
    for(ii in 1:70){
      index_groups <- obs_groups_final == final_groups[ii]
      
      LS_NS[ii] <- mean(Logscore_NS[index_groups])
    }

    LS_STAT <- numeric(70)
    
    for(ii in 1:70){
      index_groups <- obs_groups_final == final_groups[ii]
      LS_STAT[ii] <-  mean(Logscore_STAT[index_groups])
    }
    
    matrix_to_show <- t(round(cbind(c(mean(LS_STAT), mean(LS_NS)),
                                    c(sd(LS_STAT), sd(LS_NS))), 3))
    
    rownames(matrix_to_show) <- c('mean', 'se')
    colnames(matrix_to_show) <- c('M-STAT', 'M-NS')
    
    cat('logscore\n')
    print(matrix_to_show)
    cat('--------\n')
    
  }
  
  # logliks
  
  vector_to_show <- c(getLoglik(Model_STAT),getLoglik(Model_NS))
  names(vector_to_show) <- c('M-STAT','M-NS')
  
  cat('Logliks\n')
  print(round(vector_to_show,0))
  cat('--------\n')
  
  # times
  
  vector_to_show <- round(c(Time_STAT[3]/60,Time_NS/60),2)
  names(vector_to_show) <- c('M-STAT','M-NS')
  
  cat('system times\n')
  print(vector_to_show)
  cat('--------\n')
  
}
