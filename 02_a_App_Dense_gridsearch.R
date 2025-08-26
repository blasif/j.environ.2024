source('00_system.R')

load('Datasets/app_dataset.RData')

# Model B - A (better) specification of the trend

if(T){
  
  list_formulas <- list("mean" = as.formula("  ~ 1 + wind + new_wind + BIO04 + BIO15 + cloud_c + elevation + lati + long"),
                        "std.dev" = as.formula("  ~ 1 + BIO04 + BIO15 + new_wind + log_elevation + log_wind + log_cloud +  log_lati + log_long"),
                        "scale" = as.formula("  ~ 1 + BIO04 + BIO15 + new_wind +  log_elevation + log_wind + log_cloud +  log_lati + log_long"),
                        "aniso" = as.formula("  ~ 1 + BIO04 + BIO15 + new_wind +  log_elevation + log_wind + log_cloud +  log_lati + log_long"),
                        "tilt" = as.formula("  ~ 1 + BIO04 + BIO15 + new_wind + log_elevation + log_wind + log_cloud +  log_lati + log_long"),
                        "smooth" = as.formula(" ~ 1 + BIO04 + BIO15 + new_wind +  log_elevation + log_wind + log_cloud +  log_lati + log_long"),
                        "nugget" = -Inf)
  
  test_coco_classic <- coco(type = 'dense',
                            data = all_dfs[[2]],
                            locs = as.matrix(all_dfs[[2]][, c("long","lati")]),
                            z = all_dfs[[2]]$prec,
                            model.list = list_formulas,
                            info = list('smooth.limits' = c(0.5, 2.0)
                            )
  )
  
  if(T){
    
    tmp_DM <- getDesignMatrix(test_coco_classic@model.list,data = test_coco_classic@data)
    
    std_stuff <- getScale(tmp_DM$model.matrix)$std.covs
    
    tmp_lm <- lm(z ~ 1 + wind + new_wind + BIO04 + BIO15 + cloud_c + 
                   elevation + lati + long, data = cbind.data.frame('z' = test_coco_classic@z[,1],
                                                                    as.data.frame(std_stuff)))
    
    coefs_lm <- coef(tmp_lm)
    
    boundaries_B <- getBoundariesV2(coco.object = test_coco_classic,
                                    mean.limits = c(-Inf, 0, Inf),
                                    std.dev.limits = c(-2, 0, 2),
                                    scale.limits = c(-2, 0, 2),
                                    aniso.limits =  c(-2, 0, 2),
                                    tilt.limits =  c(-2, 0, 2),
                                    smooth.limits = c(-2, 0, 2),
                                    nugget.limits = c(-2, 0, 2))
    
    boundaries_B$theta_init[1:length(coefs_lm)] <- coefs_lm
    
    first_var <- which(names(boundaries_B$theta_init) == "std.dev.limits")[1]
    n_var <- length(which(names(boundaries_B$theta_init) == "std.dev.limits")) - 1
    
    first_range <- which(names(boundaries_B$theta_init) == "scale.limits")[1]
    n_range <- length(which(names(boundaries_B$theta_init) == "scale.limits")) - 1
    
    first_aniso <- which(names(boundaries_B$theta_init) == "aniso.limits")[1]
    n_aniso <- length(which(names(boundaries_B$theta_init) == "aniso.limits")) - 1
    
    first_tilt <- which(names(boundaries_B$theta_init) == "tilt.limits")[1]
    n_tilt <- length(which(names(boundaries_B$theta_init) == "tilt.limits")) - 1
    
    first_smooth <- which(names(boundaries_B$theta_init) == "smooth.limits")[1]
    n_smooth <- length(which(names(boundaries_B$theta_init) == "smooth.limits")) - 1
    
    boundaries_B$theta_upper[c(first_var, first_range)] <- c(3, 3)
    boundaries_B$theta_lower[c(first_var, first_range)] <- c(-3, -3)
    
    boundaries_B$theta_init[first_range] <- (log(sd(tmp_lm$residuals)) - log(sd(c(dist(test_coco_classic@locs)))))/2
    boundaries_B$theta_init[first_var] <- (log(sd(tmp_lm$residuals))  + log(sd(c(dist(test_coco_classic@locs)))))/2
    
    boundaries_B$theta_upper[first_smooth] <- 2
    boundaries_B$theta_lower[first_smooth] <- -3.5
    boundaries_B$theta_init[first_smooth] <- 0
    
    boundaries_B$theta_init[1] <- coefs_lm[1]
    boundaries_B$theta_upper[1] <- boundaries_B$theta_init[1] + 5
    boundaries_B$theta_lower[1] <- boundaries_B$theta_init[1] - 5
    
  }
  
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
    
    newdataset_hyper <- newdataset[which(groups %in% sample_to_tune_hyperparameters),]
    newlocs_hyper <- newdataset[which(groups %in% sample_to_tune_hyperparameters),c("long","lati")]
    z_values_hyper <- all_dfs[[1]]$prec[which(groups %in% sample_to_tune_hyperparameters)]
    
  }
  
  lambda_Sigma <- c(0, 0.2, 0.4)
  lambda_betas <- c(0, 0.05, 0.1)
  lambda_reg <- c(0.01, 0.02, 0.03)
  
  save_results <- list()
  
  for(ii in 3:3){
    
    save_results[[ii]] <- list()
    
    for(jj in 1:3){
      
      save_results[[ii]][[jj]] <- list()
      
      for(zz in 1:3){
        
        cat("ii:",ii," jj:",jj," zz:",zz, " ")
        
        test_coco_classic@info$lambda.Sigma <- lambda_Sigma[ii]
        test_coco_classic@info$lambda.betas <- lambda_betas[jj]
        test_coco_classic@info$lambda.reg <- lambda_reg[zz]
        
        Model_coco <- cocoOptim(coco.object = test_coco_classic,
                                boundaries = boundaries_B,
                                optim.type = 'ml')
        
        Pred_B <- cocoPredict(coco.object = Model_coco, 
                              newdataset = newdataset_hyper, 
                              newlocs = as.matrix(newlocs_hyper),
                              type = 'pred',
                              index.pred = 1)
        
        pred_median_CRPS <- mean(getCRPS(z.pred = z_values_hyper, 
                                           mean.pred = Pred_B$systematic + Pred_B$stochastic, 
                                           sd.pred = Pred_B$sd.pred))

        save_results[[ii]][[jj]][[zz]] <- list(Model_coco, 
                                               pred_median_CRPS,
                                               ii,
                                               jj,
                                               zz, 
                                               test_coco_classic@info$lambda.Sigma,
                                               test_coco_classic@info$lambda.betas ,
                                               test_coco_classic@info$lambda.reg)
        
      }
    }
  }
  
  to_stack <- data.frame('ii' = 1,'jj' = 1,'zz' =1,
                         'lambda_Sigma' = lambda_Sigma[1], 
                         'lambda_betas' = lambda_betas[1] , 
                         'lambda_reg' = lambda_reg[1], 
                         'CRPS' =  save_results[[1]][[1]][[1]][[2]])
  
  for(ii in 1:3){
    print(ii)
    for(jj in 1:3){
      for(zz in 1:3){
        
        to_stack <- rbind(to_stack, c(ii,jj,zz,lambda_Sigma[ii],lambda_betas[jj],lambda_reg[zz],save_results[[ii]][[jj]][[zz]][[2]]))
        
      }
    }
  }
  
  to_stack <- to_stack[-1,]
  
  index_best <- which.min(to_stack$CRPS)
  
  Model_NS <- save_results[[ to_stack$ii[index_best] ]][[to_stack$jj[index_best] ]][[ to_stack$zz[index_best] ]][[1]]
  
  # Compute time simplified model
  
  Model_NS@info$lambda.betas <- 0
  Model_NS@info$lambda.Sigma <- 0
  
  Time_NS <- system.time(cocoOptim(coco.object = Model_NS,
                          boundaries = Model_NS@info$boundaries,
                          optim.type = 'ml', ncores = 7))[3]
  
  Model_NS@info$lambda.betas <- to_stack$lambda_betas[index_best]
  Model_NS@info$lambda.Sigma <- to_stack$lambda_Sigma[index_best]
  
  HESS_NS <- getHessian(Model_NS)
  
  save(Time_NS, Model_NS, HESS_NS, file = 'RData/Model_NS.RData')
  
}
