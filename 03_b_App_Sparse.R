
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
  
  Time_T_B <- system.time({ Model_T_B <- cocoOptim(coco.object = OT_classic,
                                                    ncores = "auto",
                                                    optim.type = 'ml',
                                                    boundaries = boundaries_B,
                                                    optim.control = list(control = list(trace = 5,
                                                                                       factr = 1e-7/.Machine$double.eps)))})
  
  save(Model_T_B, Time_T_B, file = 'RData/Model_T_B.RData')
  
}

# Predictions
if(T){
  
  load('RData/Model_T_A.Rdata')
  
  newdataset <- all_dfs[[1]]
  newlocs <- all_dfs[[1]][,1:2]
  
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
  
  load("RData/Model_T_A.RData")
  
  Pred_T_A <- cocoPredict(coco.object = Model_T_A, 
                                newdataset = newdataset_final, 
                                newlocs = newlocs_final,
                                type = 'pred')
  
  Pred_T_B <- cocoPredict(coco.object = Model_T_B, 
                                newdataset = newdataset_final, 
                                newlocs = newlocs_final,
                                type = 'pred')

  CRPS_T_A <- getCRPS(z_values_final, mean.pred = Pred_T_A$systematic + Pred_A$stochastic, sd.pred = Pred_T_A$sd.pred)
  CRPS_T_B <- getCRPS(z_values_final, mean.pred = Pred_T_B$systematic + Pred_T_B$stochastic, sd.pred = Pred_T_B$sd.pred)
  
  Logscore_T_A <- getLogScore(z_values_final, mean.pred = Pred_T_A$systematic + Pred_T_A$stochastic, sd.pred = Pred_T_A$sd.pred)
  Logscore_T_B <- getLogScore(z_values_final, mean.pred = Pred_T_B$systematic + Pred_T_B$stochastic, sd.pred = Pred_T_B$sd.pred)
  
  z_std_T_A <- (Pred_T_A$systematic + Pred_T_A$stochastic - z_values_final) / Pred_T_A$sd.pred
  z_std_T_B <- (Pred_T_B$systematic + Pred_T_B$stochastic - z_values_final) / Pred_T_B$sd.pred
  
  save(Pred_T_A, Pred_T_B,
       CRPS_T_A, CRPS_T_B, 
       Logscore_T_A, Logscore_T_B,
       z_std_T_A, z_std_T_B,
       file = 'RData/Pred_Taper.RData')
  
}
  
# Metrics

load('RData/Pred_Taper.RData')

if(T){
  
  # RMSPE
  
  if(T){
    
    RMSPE_T_A <- numeric(70)
    
    for(ii in 1:70){
      index_groups <- obs_groups_final == final_groups[ii]
      RMSPE_T_A[ii] <- sqrt(mean((Pred_T_A$systematic[index_groups] + Pred_T_A$stochastic[index_groups] - z_values_final[index_groups])^2))
    }
    
    mean(RMSPE_T_A)
    sd(RMSPE_T_A)
    
    RMSPE_T_B <- numeric(70)
    
    for(ii in 1:70){
      index_groups <- obs_groups_final == final_groups[ii]
      RMSPE_T_B[ii] <- sqrt(mean((Pred_T_B$systematic[index_groups] + Pred_T_B$stochastic[index_groups] - z_values_final[index_groups])^2))
    }
    
    mean(RMSPE_T_B)
    sd(RMSPE_T_B)
    
    matrix_to_show <- t(round(cbind(c(mean(RMSPE_T_B),mean(RMSPE_T_A)),
                                    c(sd(RMSPE_T_B),sd(RMSPE_T_A))),2))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS')
    
    cat('RMSPE\n')
    print(matrix_to_show)
    cat('------\n')
    
  }
  
  # CRPS
  
  if(T){
    
    CRPS_T_A_holdout <- numeric(70)
    
    for(ii in 1:70){
      CRPS_T_A_holdout[ii] <- mean(CRPS_T_A[obs_groups_final == final_groups[ii]])
    }
    
    mean(CRPS_T_A_holdout)
    sd(CRPS_T_A_holdout)
    
    CRPS_T_B_holdout <- numeric(70)
    
    for(ii in 1:70){
      CRPS_T_B_holdout[ii] <- mean(CRPS_T_B[obs_groups_final == final_groups[ii]])
    }
    
    mean(CRPS_T_B_holdout)
    sd(CRPS_T_B_holdout)
    
    matrix_to_show <- t(round(cbind(c(mean(CRPS_T_B_holdout),mean(CRPS_T_A_holdout)),
                                    c(sd(CRPS_T_B_holdout),sd(CRPS_T_A_holdout))), 3))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS')
    
    cat('CRPS\n')
    print(matrix_to_show)
    cat('------\n')
    
  }
  
  # Q_0.9 CRPS
  
  if(T){
    
    Q_CRPS_T_A <- numeric(70)
    
    for(ii in 1:70){
      Q_CRPS_T_A[ii] <- quantile(CRPS_T_A[obs_groups_final == final_groups[ii]], probs = 0.95)
    }
    
    mean(Q_CRPS_T_A)
    sd(Q_CRPS_T_A)
    
    Q_CRPS_T_B <- numeric(70)
    
    for(ii in 1:70){
      Q_CRPS_T_B[ii] <- quantile(CRPS_T_B[obs_groups_final == final_groups[ii]],probs = 0.95)
    }
    
    mean(Q_CRPS_T_B)
    sd(Q_CRPS_T_B)
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS')
    
    cat('0.95 Quantile CRPS\n')
    print(matrix_to_show)
    cat('---------\n')
    
  }
  
  # D_n
  
  if(T){
    
    KS_T_A <- numeric(70)
    
    for(ii in 1:70){
      KS_T_A[ii] <- ks.test(z_std_T_A[obs_groups_final == final_groups[ii]],y = pnorm)$statistic
    }
    
    mean(KS_T_A)
    sd(KS_T_A)
    
    KS_T_B <- numeric(70)
    
    for(ii in 1:70){
      KS_T_B[ii] <- ks.test(z_std_T_B[obs_groups_final == final_groups[ii]],
                          y = pnorm)$statistic
    }
    
    mean(KS_T_B)
    sd(KS_T_B)  
    
    matrix_to_show <- t(round(cbind(c(mean(KS_T_B),mean(KS_T_A)),
                                    c(sd(KS_T_B),sd(KS_T_A))),2))
    
    rownames(matrix_to_show) <- c('mean', 'se')
    colnames(matrix_to_show) <- c('M-STAT', 'M-NS')
    
    cat('KS statistic \n')
    print(matrix_to_show)
    cat('--------\n')
    
  }
  
  # Coverage_prob
  
  if(T){
    
    Cov_prob_T_A <- numeric(70)
    
    for(ii in 1:70){
      
      index_groups <- obs_groups_final == final_groups[ii]
      
      upper_bound <- Pred_T_A$systematic[index_groups] + 
        Pred_T_A$stochastic[index_groups] +
        qnorm(1 - 0.025) * Pred_T_A$sd.pred[index_groups]
      
      lower_bound <- Pred_T_A$systematic[index_groups] + 
        Pred_T_A$stochastic[index_groups] -
        qnorm(1 - 0.025) * Pred_T_A$sd.pred[index_groups]
      
      Cov_prob_T_A[ii] <-  1 - length(which(z_values_final[index_groups] < lower_bound |
                                            z_values_final[index_groups] > upper_bound  )) /
        length(z_values_final[index_groups])
      
    }
    
    mean(Cov_prob_T_A)
    sd(Cov_prob_T_A)
    
    Cov_prob_T_B <- numeric(70)
    
    for(ii in 1:70){
      
      index_groups <- obs_groups_final == final_groups[ii]
      
      upper_bound <- Pred_T_B$systematic[index_groups] + 
        Pred_T_B$stochastic[index_groups] +
        qnorm(1 - 0.025) * Pred_T_B$sd.pred[index_groups]
      
      lower_bound <- Pred_T_B$systematic[index_groups] + 
        Pred_T_B$stochastic[index_groups] -
        qnorm(1 - 0.025) * Pred_T_B$sd.pred[index_groups]
      
      Cov_prob_T_B[ii] <-  1 - length(which(z_values_final[index_groups] < lower_bound |
                                            z_values_final[index_groups] > upper_bound  )) /
        length(z_values_final[index_groups])
      
      
    }
    
    mean(Cov_prob_T_B)
    sd(Cov_prob_T_B)
    
    matrix_to_show <- t(round(cbind(c(mean(Cov_prob_T_B),mean(Cov_prob_T_A)),
                                    c(sd(Cov_prob_T_B),sd(Cov_prob_T_A))),2))
    
    rownames(matrix_to_show) <- c('mean','se')
    colnames(matrix_to_show) <- c('M-STAT','M-NS')
    
    
    cat('coverage probability\n')
    print(matrix_to_show)
    cat('--------------\n')
    
  }
  
  # log-score
  
  if(T){
    
    LS_T_A <- numeric(70)
    
    for(ii in 1:70){
      index_groups <- obs_groups_final == final_groups[ii]
      
      LS_T_A[ii] <- mean(Logscore_T_A[index_groups])
    }
    
    mean(LS_T_A)
    sd(LS_T_A)
    
    LS_T_B <- numeric(70)
    
    for(ii in 1:70){
      index_groups <- obs_groups_final == final_groups[ii]
      LS_T_B[ii] <-  mean(Logscore_T_B[index_groups])
    }
    
    matrix_to_show <- t(round(cbind(c(mean(LS_T_B), mean(LS_T_A)),
                                    c(sd(LS_T_B), sd(LS_T_A))), 2))
    
    rownames(matrix_to_show) <- c('mean', 'se')
    colnames(matrix_to_show) <- c('M-STAT', 'M-NS')
    
    cat('logscore\n')
    print(matrix_to_show)
    cat('--------\n')
    
  }
  
  # logliks
  
  vector_to_show <- c(getLoglik(Model_T_A),getLoglik(Model_T_B))
  names(vector_to_show) <- c('M-STAT','M-NS')
  
  cat('Logliks\n')
  print(vector_to_show)
  cat('--------\n')
  
  # times
  
  vector_to_show <- round(c(Time_B[3]/60,Time_A/60),2)
  names(vector_to_show) <- c('M-STAT','M-NS')
  
  cat('system times\n')
  print(vector_to_show)
  cat('--------\n')
  
}