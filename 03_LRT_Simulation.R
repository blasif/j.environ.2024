
# loglikelihood ratio test

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

for(ww in 1:4){
  
  set.seed(31416)
  locs <- cbind(runif(sample_sizes[ww]) - 0.5, runif(sample_sizes[ww]) - 0.5)
  
  df_training_mc <- data.frame('x' = locs[,1],
                               'y' = locs[,2],
                               'cov_x' = Vectorize(create_covariate)(locs[,1]))
  
  list_formulas <- list("mean" = 0,
                        "std.dev" = as.formula(" ~ 1"),
                        "scale" = as.formula(" ~ 1"),
                        "aniso" = 0,
                        "tilt" = 0,
                        "smooth" = 1.5,
                        "nugget" = -Inf)
  
  CoCo_object_nn <- coco(type = 'dense',
                          model.list = list_formulas,
                          locs = as.matrix(locs),
                          z = numeric(length = dim(df_training_mc)[1]),
                          data = df_training_mc)
  
  true_pars <-  c(0, log(0.15))
  
  a_ver_la_repe_nn <- cocoSim(coco.object = CoCo_object_nn, n = 1000, type = 'classic',
                               pars = true_pars, seed = 181222 * 2,
                               standardize = TRUE)
  
  neg2loglikelihood_matrix <- matrix(ncol = 3, nrow = 1000)
  
  for(ii in 1:1000){

    list_formulas <- list("mean" = as.formula(" ~ 1"),
                          "std.dev" = as.formula(" ~ 1 + cov_x"),
                          "scale" = as.formula(" ~ 1  + cov_x"),
                          "aniso" = 0,
                          "tilt" = 0,
                          "smooth" = 1.5,
                          "nugget" = -Inf)
    
    suppressWarnings(CoCo_object_ms <- coco(type = 'dense',
                            model.list = list_formulas,
                            locs = as.matrix(locs),
                            z =  c(a_ver_la_repe_nn[ii, ]),
                            data = df_training_mc))
    
    tmp_boundaries <- getBoundaries(CoCo_object_ms, lower.value = -5, 2.5)
    
    tmp_boundaries$theta_upper[c(3, 5)] <- c(1, 1)
    tmp_boundaries$theta_lower[c(3, 5)] <- c(-1, -1)
    
    optim_ms_mle_full <- cocoOptim(coco.object = CoCo_object_ms, 
                                   ncores = ncores,
                                   boundaries = tmp_boundaries,
                                   optim.type = 'mle')
    
    list_formulas <- list("mean" = as.formula(" ~ 1"),
                          "std.dev" = as.formula(" ~ 1"),
                          "scale" = as.formula(" ~ 1  + cov_x"),
                          "aniso" = 0,
                          "tilt" = 0,
                          "smooth" = 1.5,
                          "nugget" = -Inf)
    
    suppressWarnings(CoCo_object_ms <- coco(type = 'dense',
                            model.list = list_formulas,
                            locs = as.matrix(locs),
                            z =  c(a_ver_la_repe_nn[ii, ]),
                            data = df_training_mc))
    
    tmp_boundaries <- getBoundaries(CoCo_object_ms,lower.value = -5, 2.5)
    tmp_boundaries$theta_upper[c(4)] <- c(1)
    tmp_boundaries$theta_lower[c(4)] <- c(-1)
    
    optim_ms_mle <- cocoOptim(coco.object = CoCo_object_ms, 
                               ncores = ncores, 
                               boundaries = tmp_boundaries,
                               optim.type = 'mle')
    
    list_formulas <- list("mean" = as.formula(" ~ 1"),
                          "std.dev" = as.formula(" ~ 1"),
                          "scale" = as.formula(" ~ 1"),
                          "aniso" = 0,
                          "tilt" = 0,
                          "smooth" = 1.5,
                          "nugget" = -Inf)
    
    suppressWarnings(CoCo_object_ms <- coco(type = 'dense',
                            model.list = list_formulas,
                            locs = as.matrix(locs),
                            z =  c(a_ver_la_repe_nn[ii, ]),
                            data = df_training_mc))
    
    tmp_boundaries <- getBoundaries(CoCo_object_ms, lower.value = -3, 3)
    
    optim_ms_mle_r <- cocoOptim(coco.object = CoCo_object_ms, 
                                 ncores = ncores, 
                                 boundaries = tmp_boundaries,
                                 optim.type = 'mle')
    
    neg2loglikelihood_matrix[ii,] <- c(optim_ms_mle_full@output$value, 
                                       optim_ms_mle@output$value, 
                                       optim_ms_mle_r@output$value)
    
  }
  
  save(neg2loglikelihood_matrix, 
       file = paste0('RData/lrt_sim_results.', sample_sizes[ww],'.RData'))
  
}

matrix_one <- matrix(nrow=4,ncol=3)
matrix_two <- matrix(nrow=4,ncol=3)

alpha_vector <- c(0.9,0.95,0.99)

ks_values <- matrix(nrow=4,ncol=2)

for(ii in 1:4){
  
  load(paste0('RData/lrt_sim_results.', sample_sizes[ii], '.RData'))
  
  ks_values[ii,1] <- ks.test(neg2loglikelihood_matrix[, 3] - neg2loglikelihood_matrix[, 2], y = function(x){pchisq(x, df = 1)})$statistic
  ks_values[ii,2] <- ks.test(neg2loglikelihood_matrix[, 3] - neg2loglikelihood_matrix[, 1], y = function(x){pchisq(x, df = 2)})$statistic
    
  for(jj in 1:3){
    
    matrix_one[ii,jj] <- sum((neg2loglikelihood_matrix[, 3] - neg2loglikelihood_matrix[, 2])[!is.na(neg2loglikelihood_matrix[,1])] > 
                              qchisq(alpha_vector[jj], df = 1)) / sum(!is.na(neg2loglikelihood_matrix[, 1]))
    
    matrix_two[ii,jj] <- sum((neg2loglikelihood_matrix[, 3] - neg2loglikelihood_matrix[, 1])[!is.na(neg2loglikelihood_matrix[,1])] > 
                               qchisq(alpha_vector[jj], df = 2)) / sum(!is.na(neg2loglikelihood_matrix[, 1]))
      
  }

}

# q = 1
matrix_one

# q = 2
matrix_two

# ks
round(ks_values,3)
