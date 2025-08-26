
source('00_system.R')
load('Datasets/app_dataset.RData')

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
                      locs = as.matrix(all_dfs[[3]][, c("long","lati")]),
                      z = all_dfs[[3]]$prec,
                      model.list = list_formulas,
                      info = list('lambda.reg' = 0.0125,
                                  'taper' = spam::cov.wend1,
                                  'delta' = 0.20, 
                                  'smooth.limits' = c(0.25, 1.5))
  )
  
  # Tailored boundaries
  
  if(T){
    
    boundaries_B <- getBoundariesV4(OT_classic)
    
    boundaries_B$theta_upper['std.dev.limits'] <- 5
    boundaries_B$theta_lower['std.dev.limits'] <- -5
    
    boundaries_B$theta_upper['scale.limits'] <- 5
    boundaries_B$theta_lower['scale.limits'] <- -5
    
  }
  
  Time_T_STAT <- system.time({ Model_T_STAT <- cocoOptim(coco.object = OT_classic,
                                                    ncores = 7,
                                                    optim.type = 'ml',
                                                    boundaries = boundaries_B,
                                                    optim.control = list(control = list(factr = 1e-7/.Machine$double.eps)))})
  
  save(Model_T_STAT, Time_T_STAT, file = 'RData/Model_T_STAT.RData')
  
}

# Predictions
if(T){
  
  load('RData/Model_T_NS.Rdata')
  
  newdataset <- all_dfs[[1]]
  newlocs <- all_dfs[[1]][, 1:2]
  
  z_values <- all_dfs[[1]]$prec
  
  # Create holdouts
  if(T){
    
    newdataset <- all_dfs[[1]]
    newlocs <- all_dfs[[1]][ ,c("long","lati")]
    z_values <- all_dfs[[1]]$prec
    
    set.seed(100621)
    hetero_holdouts <- kmeans(as.data.frame(scale(newdataset[,c(1,2,4:9)])),centers = 100,iter.max = 100)
    groups <- as.factor(hetero_holdouts$cluster)
    quilt.plot(newlocs, hetero_holdouts$cluster, nx = 150, ny = 150)
    
    sample_to_tune_hyperparameters <- sample(1:100,30)
    
    newdataset_final <- newdataset[!(groups %in% sample_to_tune_hyperparameters),]
    newlocs_final <- as.matrix(newdataset[!(groups %in% sample_to_tune_hyperparameters), c("long","lati")])
    z_values_final <- all_dfs[[1]]$prec[!(groups %in% sample_to_tune_hyperparameters)]
    
    final_groups <- c(1:100)[!(c(1:100) %in% sample_to_tune_hyperparameters)]
    
    obs_groups_final <- groups[groups %in% final_groups]
    
  }
  
  load('RData/Model_T_NS.Rdata')
  
  Pred_T_NS <- cocoPredict(coco.object = Model_T_NS, 
                                newdataset = newdataset_final, 
                                newlocs = newlocs_final,
                                type = 'pred')
  
  Pred_T_STAT <- cocoPredict(coco.object = Model_T_STAT, 
                                newdataset = newdataset_final, 
                                newlocs = newlocs_final,
                                type = 'pred')

  CRPS_T_NS <- getCRPS(z_values_final, mean.pred = Pred_T_NS$systematic + Pred_T_NS$stochastic, sd.pred = Pred_T_NS$sd.pred)
  CRPS_T_STAT <- getCRPS(z_values_final, mean.pred = Pred_T_STAT$systematic + Pred_T_STAT$stochastic, sd.pred = Pred_T_STAT$sd.pred)
  
  Logscore_T_NS <- getLogScore(z_values_final, mean.pred = Pred_T_NS$systematic + Pred_T_NS$stochastic, sd.pred = Pred_T_NS$sd.pred)
  Logscore_T_STAT <- getLogScore(z_values_final, mean.pred = Pred_T_STAT$systematic + Pred_T_STAT$stochastic, sd.pred = Pred_T_STAT$sd.pred)
  
  z_std_T_NS <- (Pred_T_NS$systematic + Pred_T_NS$stochastic - z_values_final) / Pred_T_NS$sd.pred
  z_std_T_STAT <- (Pred_T_STAT$systematic + Pred_T_STAT$stochastic - z_values_final) / Pred_T_STAT$sd.pred
  
  save(Pred_T_NS, Pred_T_STAT,
       CRPS_T_NS, CRPS_T_STAT, 
       Logscore_T_NS, Logscore_T_STAT,
       z_std_T_NS, z_std_T_STAT,
       file = 'RData/Pred_Taper.RData')
  
}
  
# Metrics

load('RData/Pred_Taper.RData')

if(T){
  
  # RMSPE
  
  if(T){
    
    RMSPE_T_NS <- numeric(70)
    
    for(ii in 1:70){
      index_groups <- obs_groups_final == final_groups[ii]
      RMSPE_T_NS[ii] <- sqrt(mean((Pred_T_NS$systematic[index_groups] + Pred_T_NS$stochastic[index_groups] - z_values_final[index_groups])^2))
    }

    RMSPE_T_STAT <- numeric(70)
    
    for(ii in 1:70){
      index_groups <- obs_groups_final == final_groups[ii]
      RMSPE_T_STAT[ii] <- sqrt(mean((Pred_T_STAT$systematic[index_groups] + Pred_T_STAT$stochastic[index_groups] - z_values_final[index_groups])^2))
    }

    matrix_to_show <- t(round(cbind(c(mean(RMSPE_T_STAT),mean(RMSPE_T_NS)),
                                    c(sd(RMSPE_T_STAT),sd(RMSPE_T_NS))),3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-T-STAT','M-T-NS')
    
    cat('RMSPE\n')
    print(matrix_to_show)
    cat('------\n')
    
  }
  
  # CRPS
  
  if(T){
    
    CRPS_T_NS_holdout <- numeric(70)
    
    for(ii in 1:70){
      CRPS_T_NS_holdout[ii] <- mean(CRPS_T_NS[obs_groups_final == final_groups[ii]])
    }

    CRPS_T_STAT_holdout <- numeric(70)
    
    for(ii in 1:70){
      CRPS_T_STAT_holdout[ii] <- mean(CRPS_T_STAT[obs_groups_final == final_groups[ii]])
    }

    matrix_to_show <- t(round(cbind(c(mean(CRPS_T_STAT_holdout),mean(CRPS_T_NS_holdout)),
                                    c(sd(CRPS_T_STAT_holdout),sd(CRPS_T_NS_holdout))), 3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-T-STAT','M-T-NS')
    
    cat('CRPS\n')
    print(matrix_to_show)
    cat('------\n')
    
  }
  
  # Q_0.9 CRPS
  
  if(T){
    
    Q_CRPS_T_NS <- numeric(70)
    
    for(ii in 1:70){
      Q_CRPS_T_NS[ii] <- quantile(CRPS_T_NS[obs_groups_final == final_groups[ii]], probs = 0.95)
    }
  
    Q_CRPS_T_STAT <- numeric(70)
    
    for(ii in 1:70){
      Q_CRPS_T_STAT[ii] <- quantile(CRPS_T_STAT[obs_groups_final == final_groups[ii]],probs = 0.95)
    }
    
    matrix_to_show <- t(round(cbind(c(mean(Q_CRPS_T_STAT),mean(Q_CRPS_T_NS)),
                                    c(sd(Q_CRPS_T_STAT),sd(Q_CRPS_T_NS))), 3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-T-STAT','M-T-NS')
    
    cat('0.95 Quantile CRPS\n')
    print(matrix_to_show)
    cat('---------\n')
    
  }
  
  # D_n
  
  if(T){
    
    KS_T_NS <- numeric(70)
    
    for(ii in 1:70){
      KS_T_NS[ii] <- ks.test(z_std_T_NS[obs_groups_final == final_groups[ii]],y = pnorm)$statistic
    }

    KS_T_STAT <- numeric(70)
    
    for(ii in 1:70){
      KS_T_STAT[ii] <- ks.test(z_std_T_STAT[obs_groups_final == final_groups[ii]],
                          y = pnorm)$statistic
    }

    matrix_to_show <- t(round(cbind(c(mean(KS_T_STAT),mean(KS_T_NS)),
                                    c(sd(KS_T_STAT),sd(KS_T_NS))),3))
    
    rownames(matrix_to_show) <- c('mean', 'se')
    colnames(matrix_to_show) <- c('M-T-STAT', 'M-T-NS')
    
    cat('KS statistic \n')
    print(matrix_to_show)
    cat('--------\n')
    
  }
  
  # Coverage_prob
  
  if(T){
    
    Cov_prob_T_NS <- numeric(70)
    
    for(ii in 1:70){
      
      index_groups <- obs_groups_final == final_groups[ii]
      
      upper_bound <- Pred_T_NS$systematic[index_groups] + 
        Pred_T_NS$stochastic[index_groups] +
        qnorm(1 - 0.025) * Pred_T_NS$sd.pred[index_groups]
      
      lower_bound <- Pred_T_NS$systematic[index_groups] + 
        Pred_T_NS$stochastic[index_groups] -
        qnorm(1 - 0.025) * Pred_T_NS$sd.pred[index_groups]
      
      Cov_prob_T_NS[ii] <-  1 - length(which(z_values_final[index_groups] < lower_bound |
                                            z_values_final[index_groups] > upper_bound  )) /
        length(z_values_final[index_groups])
      
    }

    Cov_prob_T_STAT <- numeric(70)
    
    for(ii in 1:70){
      
      index_groups <- obs_groups_final == final_groups[ii]
      
      upper_bound <- Pred_T_STAT$systematic[index_groups] + 
        Pred_T_STAT$stochastic[index_groups] +
        qnorm(1 - 0.025) * Pred_T_STAT$sd.pred[index_groups]
      
      lower_bound <- Pred_T_STAT$systematic[index_groups] + 
        Pred_T_STAT$stochastic[index_groups] -
        qnorm(1 - 0.025) * Pred_T_STAT$sd.pred[index_groups]
      
      Cov_prob_T_STAT[ii] <-  1 - length(which(z_values_final[index_groups] < lower_bound |
                                            z_values_final[index_groups] > upper_bound  )) /
        length(z_values_final[index_groups])
      
      
    }

    matrix_to_show <- t(round(cbind(c(mean(Cov_prob_T_STAT),mean(Cov_prob_T_NS)),
                                    c(sd(Cov_prob_T_STAT),sd(Cov_prob_T_NS))),3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-T-STAT','M-T-NS')
    
    
    cat('coverage probability\n')
    print(matrix_to_show)
    cat('--------------\n')
    
  }
  
  # log-score
  
  if(T){
    
    LS_T_NS <- numeric(70)
    
    for(ii in 1:70){
      index_groups <- obs_groups_final == final_groups[ii]
      
      LS_T_NS[ii] <- mean(Logscore_T_NS[index_groups])
    }

    LS_T_STAT <- numeric(70)
    
    for(ii in 1:70){
      index_groups <- obs_groups_final == final_groups[ii]
      LS_T_STAT[ii] <-  mean(Logscore_T_STAT[index_groups])
    }
    
    matrix_to_show <- t(round(cbind(c(mean(LS_T_STAT), mean(LS_T_NS)),
                                    c(sd(LS_T_STAT), sd(LS_T_NS))), 3))
    
    rownames(matrix_to_show) <- c('mean', 'se')
    colnames(matrix_to_show) <- c('M-T-STAT', 'M-T-NS')
    
    cat('logscore\n')
    print(matrix_to_show)
    cat('--------\n')
    
  }
  
  # logliks
  
  vector_to_show <- c(getLoglik(Model_T_STAT), getLoglik(Model_T_NS))
  names(vector_to_show) <- c('M-T-STAT', 'M-T-NS')
  
  cat('Logliks\n')
  print(vector_to_show)
  cat('--------\n')
  
  # times
  
  vector_to_show <- round(c(Time_T_STAT[3]/60, Time_T_NS/60), 2)
  names(vector_to_show) <- c('M-T-STAT','M-T-NS')
  
  cat('system times\n')
  print(vector_to_show)
  cat('--------\n')
  
}

