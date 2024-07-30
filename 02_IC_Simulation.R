
source('00_system.R')

sample_sizes <- c(100, 200, 400, 800)

create_covariate <- function(x){
  
  return(sin(x*10))
  if(x < 0){
    0
  }else{
    1
  }
  
}

for(zz in 1:4){
  
  set.seed(1231231212)
  
  locs <- cbind(runif(sample_sizes[zz],min = -0.5, max = 0.5), 
                runif(sample_sizes[zz], min = -0.5, max = 0.5))
  
  df_training_mc <- data.frame('x' = locs[, 1],
                               'y' = locs[, 2],
                               'cov_x' = Vectorize(create_covariate)(locs[,1]))
  
  list_formulas <- list("mean" = 0,
                        "std.dev" = as.formula(" ~ 1 + cov_x"),
                        "scale" = as.formula(" ~ 1 + cov_x"),
                        "aniso" = 0,
                        "tilt" = 0,
                        "smooth" = 1.5,
                        "nugget" = -Inf)

  CoCo_object_nn <- coco(type = 'dense',
                          model.list = list_formulas,
                          locs = as.matrix(locs),
                          z = numeric(length = dim(df_training_mc)[1]),
                          data = df_training_mc)
  
  true_pars <-  c(0, 0.5, log(0.25), -0.25)
  
  coco_realiz <- cocoSim(coco.object = CoCo_object_nn, 
                               n = 1000, 
                               type = 'classic',
                               pars = true_pars, seed = 181222,
                               standardize = TRUE)
  
  IC_matrix <- matrix(ncol = 3, nrow = 1000)
  
  IC_list <- list()
  
  for(ii in 1:1000){
    
    list_formulas <- list("mean" = 0,
                          "std.dev" = as.formula(" ~ 1 + cov_x"),
                          "scale" = as.formula(" ~ 1 + cov_x"),
                          "aniso" = 0,
                          "tilt" = 0,
                          "smooth" = 1.5,
                          "nugget" = -Inf)
    
    suppressWarnings(CoCo_object_ms <- coco(type = 'dense',
                            model.list = list_formulas,
                            locs = as.matrix(locs),
                            z =  t(coco_realiz[ii, , drop = FALSE]),
                            data = df_training_mc))
    
    coco_boundaries <- getBoundaries(CoCo_object_ms, -2.5, 2.5)
    
    coco_boundaries$theta_upper[c(2, 4)] <- c(1.25, 1.25)
    coco_boundaries$theta_lower[c(2, 4)] <- c(-1.25, -1.25)
    
    optim_ms_mle <- cocoOptim(coco.object = CoCo_object_ms, 
                              ncores = ncores, 
                              boundaries = coco_boundaries,
                              optim.type = 'mle',
                              optim.control = list('control' = list('trace' = 1)))
    
    sd_cov_x <- solve(getHessian(optim_ms_mle,ncores = ncores))
    
    IC_list[[ii]] <- list('pars' = optim_ms_mle@output$par,
                          'hess' = sd_cov_x)
    
  }
  
  save(IC_list, file = paste0('RData/ic_sim_results.', sample_sizes[zz], '.RData'))
  
}

# results

alpha_values <- c(0.1, 0.05, 0.01)

vector_ns <- c(100, 200, 400, 800)

matrix_results_sd <- matrix(ncol = 4, nrow = 3)
matrix_results_scale <- matrix(ncol = 4, nrow = 3)

matrix_results_sum <- matrix(ncol = 4, nrow = 3)
matrix_results_diff <- matrix(ncol = 4, nrow = 3)

matrix_chi <- matrix(ncol=4,nrow=3)
matrix_chi_range <- matrix(ncol=4,nrow=3)

list_bias <- list()

for(ww in 1:4){
  
  load(file = paste0('RData/ic_sim_results.', sample_sizes[ww],'.RData'))
  
  cis_sd <- matrix(ncol = 2, nrow = 1000)
  cis_scale <- matrix(ncol = 2, nrow = 1000)
  cis_range_ns <- matrix(ncol = 1, nrow = 1000)
  cis_raw_sum <- matrix(ncol = 2, nrow = 1000)
  cis_raw_diff <- matrix(ncol = 2, nrow = 1000)
  
  n <- vector_ns[ww]
  
  list_bias[[ww]] <-  matrix(ncol=2,nrow=1000)
  
  for(jj in 1:3){
    
    for(ii in 1:length(IC_list)){
      
      # Bias
      if(jj == 1){
        list_bias[[ww]][ii,1] <- sum(IC_list[[ii]]$pars[c(2,4)])/2
        list_bias[[ww]][ii,2] <- diff(rev(IC_list[[ii]]$pars[c(2,4)]))/2
        
      }
      
      modHESS <- getModHess(optim_ms_mle,inv.hess = IC_list[[ii]]$hess)
      
      cis_sd[ii, ] <- mean(IC_list[[ii]]$pars[c(2,4)]) + c(-1,1) * qnorm(1-(alpha_values[jj]/2)) * sqrt(modHESS[2,2])
      cis_scale[ii, ] <- diff(IC_list[[ii]]$pars[c(4,2)])/2 + c(-1,1) * qnorm(1-(alpha_values[jj]/2)) * sqrt(modHESS[4,4])
      
    }
    
    matrix_results_sd[jj, ww] <- length(which(cis_sd[, 1] < true_pars[2] & cis_sd[, 2] > true_pars[2])) / length(IC_list)
    matrix_results_scale[jj, ww] <- length(which(cis_scale[, 1] < true_pars[4] & cis_scale[, 2] > true_pars[4])) / length(IC_list)

  }
  
}

# Bias

round(apply(list_bias[[1]], 2, mean), 3)
round(apply(list_bias[[2]], 2, mean), 3)
round(apply(list_bias[[3]], 2, mean), 3)
round(apply(list_bias[[4]], 2, mean), 3)

# Cov.Prob.

t(matrix_results_sd)
t(matrix_results_scale)
