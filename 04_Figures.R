rm(list = ls())
source("00_system.R")
load("Datasets/app_dataset.RData")

# Kernels

if (T) {
  
  createEllipse <- function(center, a, b, angle = 0, steps = 100) {
    theta <- seq(0, 2 * pi, length.out = steps)
    ellipse_points <- cbind(center[1] + a * cos(theta), center[2] + b * sin(theta))
    ellipse_points <- ellipse_points %*% matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2)
    return(ellipse_points)
  }
  
  if (T) {
    pdf("Figures/univariate_r.pdf", width = master_width, height = master_height)
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    plot(
      x = 0, y = 0, xlim = c(-2, 2), ylim = c(-2, 2), cex = 1,
      xlab = "x", ylab = "y", pch = 20, cex.lab = 1.5, cex.axis = 1.5
    )
    abline(h = 0, lty = 2)
    abline(v = 0, lty = 2)
    
    r_vector_sq <- c(0.25, 0.5, 1, 2)
    
    alpha_i <- pi / 2
    
    for (ii in 1:length(r_vector_sq)) {
      r_sq <- r_vector_sq[ii]
      e_1 <- 1
      e_2 <- r_sq
      
      arrows(
        x0 = 0, x1 = 1,
        y0 = 0, y1 = 0, lwd = 2, lty = 1,
        cex = 0.5, angle = 15, length = 0.1
      )
      
      arrows(
        x0 = 0, x1 = 0,
        y0 = 0, y1 = e_2, lwd = 2, lty = 1,
        cex = 0.5, angle = 15, length = 0.1
      )
      
      ellipse_points <- createEllipse(c(0, 0),
                                      a = 1,
                                      b = r_sq, angle = 0, steps = 2000
      )
      
      lines(ellipse_points, col = viridis::viridis(n = length(r_vector_sq))[ii], lwd = 2.0)
    }
    
    legend(bty = "n","bottomright",
           legend = c(expression(r^2 ~ "= 0.25", r^2 ~ "= 0.50", r^2 ~ "= 1.00", r^2 ~ "= 2.00")),
           col = viridis::viridis(length(r_vector_sq)), lty = rep(1, length(r_vector_sq)), cex = 1.25
    )
    
    dev.off()
  }
  
  if (T) {
    # univariate alpha
    
    pdf("Figures/univariate_alpha.pdf", width = master_width, height = master_height)
    
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    plot(
      x = 0, y = 0, xlim = c(-1.5, 1.5), ylim = c(-1, 1), asp = 1, cex = 1,
      xlab = "x", ylab = "y", pch = 20, cex.lab = 1.5, cex.axis = 1.5
    )
    abline(h = 0, lty = 2)
    abline(v = 0, lty = 2)
    
    alpha_vector <- c(pi / 10, pi / 6, pi / 4, pi / 3)
    
    for (ii in 1:length(alpha_vector)) {
      r <- 1
      
      alpha_i <- alpha_vector[ii]
      
      A_part <- sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2)
      
      e_1 <- 0.5 * ((r + 1) + sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2))
      e_2 <- 0.5 * ((r + 1) - sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2))
      
      a_a_1 <- 2 * r^(0.5) * cospi(alpha_i / pi)
      a_b_1 <- r - 1 - A_part
      
      e_b_1 <- 2 * r^(0.5) * cos(alpha_i)
      d_a_1_new <- r - 1 + A_part
      
      magnitude_vector_one <- sqrt((a_a_1)^2 + a_b_1^2)
      
      y_length <- sqrt(e_2) / magnitude_vector_one
      
      arrows(
        x0 = 0, x1 = (-1) * y_length * a_a_1,
        y0 = 0, y1 = (-1) * y_length * a_b_1, lwd = 2, lty = 1,
        cex = 0.5, angle = 15, length = 0.1
      )
      
      magnitude_vector_two <- sqrt((e_b_1)^2 + (d_a_1_new)^2)
      
      y_length <- sqrt(e_1) / magnitude_vector_two
      
      arrows(
        x0 = 0, x1 = y_length * e_b_1,
        y0 = 0, y1 = y_length * d_a_1_new, lwd = 2, lty = 1,
        cex = 0.5, angle = 15, length = 0.1
      )
      
      # Plot ellipse
      
      ellipse_points <- createEllipse(c(0, 0),
                                      a = sqrt((max(e_1, e_2))),
                                      b = sqrt((min(e_1, e_2))),
                                      atan2(x = e_b_1, y = d_a_1_new),
                                      steps = 2000
      )
      lines(ellipse_points, col = viridis::viridis(n = length(alpha_vector))[ii], lwd = 2.0)
    }
    
    legend(bty = "n","bottomright",
           legend = expression(omega ~ "=" ~ pi / 10, omega ~ "=" ~ pi / 6, omega ~ "=" ~ pi / 4, omega ~ "=" ~ pi / 3), col = viridis::viridis(length(alpha_vector)),
           lwd = rep(2, length(alpha_vector)), lty = rep(1, length(alpha_vector)), cex = 1.25
    )
    dev.off()
  }
  
  # How does r and omega affects the tilt?
  
  if (T) {
    alpha_vector <- log(seq(exp(0 + 1e-2), exp(pi / 2 - 1e-2), length.out = 1000))
    r_vector <- seq(0.05, 0.95, length.out = 10)

    pdf("Figures/angle_analysis.pdf", width = master_width, height = master_height)
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    vector_change_alpha <- numeric(length = length(alpha_vector))
    for (jj in 1:length(r_vector)) {
      r <- r_vector[jj]
      for (ii in 1:length(alpha_vector)) {
        alpha_i <- alpha_vector[ii]
        
        A_part <- sqrt((r + 1)^2 - 4 * r * sinpi(alpha_i / pi)^2)
        e_b_1 <- 2 * r^(0.5) * cos(alpha_i)
        d_a_1_new <- r - 1 + A_part
        
        vector_change_alpha[ii] <- atan2(x = e_b_1, y = d_a_1_new) # * 180 / pi
      }
      
      if (jj == 1) {
        plot(
          x = alpha_vector, y = vector_change_alpha, type = "l", ylim = c(0, pi / 4 + 0.1),
          col = viridis::viridis(length(r_vector))[jj], cex = 1, lwd = 2,
          xlab = expression(omega), ylab = "rotation angle", pch = 20,
          cex.lab = 1.5, cex.axis = 1.5,
          xaxt = "n", yaxt = "n"
        )
        abline(h = pi / 4, col = "blue", lty = 2)
        abline(v = pi / 2, col = "blue", lty = 2)
      } else {
        lines(
          x = alpha_vector,
          y = vector_change_alpha,
          col = viridis::viridis(length(r_vector))[jj], lwd = 2
        )
      }
    }
    
    axis(1, at = c(0, pi / 8, pi / 4, pi / 3, pi / 2), labels = c(0, expression(pi / 8), expression(pi / 4), expression(pi / 3), expression(pi / 2)), cex.lab = 1.5, cex.axis = 1.5)
    axis(2, at = c(0, pi / 12, pi / 6, pi / 4), labels = c(0, expression(pi / 12), expression(pi / 6), expression(pi / 4)), cex.lab = 1.5, cex.axis = 1.5)
    
    
    par(mar = c(4.5, 0, 0.5, 2.5))
    tmp_z <- matrix(1:10, nrow = 1)
    tmp_x <- 1
    tmp_y <- seq(min(r_vector), max(r_vector), len = 10) 
    tmp_y <- r_vector 
    image(tmp_x, tmp_y, tmp_z, col = viridis::viridis(11), axes = FALSE, xlab = "", ylab = "")
    axis(4, lwd = 2, cex.axis = 1.5, at = seq(0.05, 0.95, length.out = 10))
    
    dev.off()
  }
  
  if (T) {
    alpha_vector <- c(pi / 8, pi / 3, pi / 2.1)
    
    r_vector <- c(0.3, 0.5, 0.9)^2
    
    vec_letters <- c("a", "b", "c")
    
    for (jj in 1:length(r_vector)) {
      pdf(paste0("Figures/bivariate_alpha_r_simp_", vec_letters[jj], ".pdf"), width = master_width, height = master_height)
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      plot(
        x = 0, y = 0, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), asp = 1, cex = 1,
        xlab = "x", ylab = "y", pch = 20, cex.lab = 1.5, cex.axis = 1.5
      )
      
      abline(h = 0, lty = 2)
      abline(v = 0, lty = 2)
      
      r <- r_vector[jj]
      
      for (ii in 1:length(alpha_vector)) {
        alpha_i <- alpha_vector[ii]
        
        if (alpha_i == pi / 2) {
          print("linear independent eigenvectors")
          r <- r_vector[ii]
          
          A_stuff <- sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2)
          
          e_1 <- ((r + 1) + sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2))
          e_2 <- ((r + 1) - sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2))
          
          magnitude_vector_one <- sqrt((0)^2 + 1^2)
          y_length <- sqrt(max(e_2, e_1)) / magnitude_vector_one
          
          arrows(
            x0 = 0, x1 = ifelse(r > 1, 0, y_length),
            y0 = 0, y1 = ifelse(r > 1, y_length, 0), lwd = 2, lty = 1,
            cex = 0.5, angle = 5, length = 0.1
          )
          
          magnitude_vector_two <- sqrt((1)^2 + (0)^2)
          
          y_length <- sqrt(min(e_1, e_2)) / magnitude_vector_two
          
          arrows(
            x0 = 0, x1 = ifelse(r > 1, y_length, 0),
            y0 = 0, y1 = ifelse(r > 1, 0, y_length), lwd = 2, lty = 1,
            cex = 0.5, angle = 5, length = 0.1
          )
          
          ellipse_points <- createEllipse(c(0, 0),
                                          a = sqrt(ifelse(r > 1, (min(e_1, e_2)), max(e_1, e_2))),
                                          b = sqrt(ifelse(r > 1, (max(e_1, e_2)), min(e_1, e_2))), 0, steps = 2000
          )
          
          lines(ellipse_points, col = viridis::viridis(n = length(alpha_vector))[ii], lwd = 2)
          
          next
        }
        
        if (F) {
          A_stuff <- sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2)
          
          e_1 <- c(
            2 * sqrt(r) * cos(alpha_i),
            (r - 1 + A_stuff)
          )
          
          e_2 <- c(
            2 * sqrt(r) * cos(alpha_i),
            (r - 1 - A_stuff)
          )
          
          magnitude_1 <- sqrt(e_1[1]^2 + e_1[2]^2)
          magnitude_2 <- sqrt(e_2[1]^2 + e_2[2]^2)
          
          arrows(
            x0 = 0, x1 = e_1[1],
            y0 = 0, y1 = e_1[2], lwd = 2, lty = 1,
            cex = 0.5, angle = 15, length = 0.1
          )
          
          arrows(
            x0 = 0, x1 = e_2[1],
            y0 = 0, y1 = e_2[2], lwd = 2, lty = 1,
            cex = 0.5, angle = 15, length = 0.1
          )
          
          ellipse_points <- createEllipse(c(0, 0),
                                          a = magnitude_1,
                                          b = magnitude_2, atan2(x = e_1[2], y = e_1[1]), steps = 2000
          )
          
          lines(ellipse_points, col = viridis::viridis(n = length(alpha_vector))[ii], lwd = 2)
          
          e_1 <- 0.5 * ((r + 1) + sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2))
          e_2 <- 0.5 * ((r + 1) - sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2))
          
          a_a_1 <- 2 * r^(0.5) * cospi(alpha_i / pi)
          a_b_1 <- r - 1 - A_stuff
          
          e_b_1 <- 2 * r^(0.5) * cos(alpha_i)
          d_a_1_new <- r - 1 + A_stuff
          
          magnitude_vector_one <- sqrt((a_a_1)^2 + a_b_1^2)
          
          y_length <- sqrt(e_2) / magnitude_vector_one
          
          magnitude_vector_two <- sqrt((e_b_1)^2 + (d_a_1_new)^2)
          
          y_length <- sqrt(e_1) / magnitude_vector_two
          
          arrows(
            x0 = 0, x1 = y_length * e_b_1,
            y0 = 0, y1 = y_length * d_a_1_new, lwd = 2, lty = 1,
            cex = 0.5, angle = 15, length = 0.1
          )
          
          # Plot ellipse

          ellipse_points <- createEllipse(c(0, 0),
                                          a = sqrt((max(e_1, e_2))),
                                          b = sqrt((min(e_1, e_2))), atan2(x = e_b_1, y = d_a_1_new), steps = 2000
          )
          
          lines(ellipse_points, col = viridis::viridis(n = length(alpha_vector))[ii], lwd = 2)
        }
        
        A_one <- sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2)
        
        e_1 <- 0.5 * ((r + 1) + sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2))
        e_2 <- 0.5 * ((r + 1) - sqrt((r + 1)^2 - 4 * r * sin(alpha_i)^2))
        
        a_a_1 <- 2 * r^(0.5) * cospi(alpha_i / pi)
        a_b_1 <- r - 1 - A_one
        
        e_b_1 <- 2 * r^(0.5) * cos(alpha_i)
        d_a_1_new <- r - 1 + A_one
        
        magnitude_vector_one <- sqrt((a_a_1)^2 + a_b_1^2)
        
        y_length <- sqrt(e_2) / magnitude_vector_one
        
        arrows(
          x0 = 0, x1 = (-1) * y_length * a_a_1,
          y0 = 0, y1 = (-1) * y_length * a_b_1, lwd = 2, lty = 1,
          cex = 0.5, angle = 15, length = 0.1
        )
        
        magnitude_vector_two <- sqrt((e_b_1)^2 + (d_a_1_new)^2)
        
        y_length <- sqrt(e_1) / magnitude_vector_two
        
        arrows(
          x0 = 0, x1 = y_length * e_b_1,
          y0 = 0, y1 = y_length * d_a_1_new, lwd = 2, lty = 1,
          cex = 0.5, angle = 15, length = 0.1
        )
        
        # Plot ellipse
        
        ellipse_points <- createEllipse(c(0, 0),
                                        a = sqrt((max(e_1, e_2))),
                                        b = sqrt((min(e_1, e_2))), atan2(x = e_b_1, y = d_a_1_new), steps = 2000
        )
        
        lines(ellipse_points, col = viridis::viridis(n = length(alpha_vector))[ii], lwd = 2)
      }
      
      legend(bty = "n","bottomright",
             legend = c(expression(omega ~ "=" ~ pi / 8), expression(omega ~ "=" ~ pi / 3), expression(omega ~ "=" ~ 2 ~ pi / 10)), col = viridis::viridis(length(alpha_vector)),
             lwd = rep(2, length(alpha_vector)), lty = rep(1, length(alpha_vector)), cex = 1.25
      )
      dev.off()
    }
  }
}

# Noisy / Categorical covariates plots

if (T) {
  
  ### Noisy how affect
  
  locs_x <- seq(-10, 10, by = 0.05)
  locs <- cbind(locs_x, rep(0, length(locs_x)))
  
  create_covariate <- function(x) {
    if (x < 0) {
      tmp_x <- 1 / (1 + exp(-x - 5)) * 2
      tmp_x <- tmp_x + abs(rnorm(1, 0, sd = tmp_x * 0.25))
    } else {
      tmp_x <- exp(-x + 5) / (1 + exp(-x + 5)) * 2
      tmp_x <- tmp_x + abs(rnorm(1, 0, sd = tmp_x * 0.25))
    }
  }
  
  v_cv <- Vectorize(create_covariate)
  
  create_covariate_no_noise <- function(x) {
    if (x < 0) {
      tmp_x <- 1 / (1 + exp(-x - 5)) * 2
    } else {
      tmp_x <- exp(-x + 5) / (1 + exp(-x + 5)) * 2
    }
  }
  
  v_cv_nn <- Vectorize(create_covariate_no_noise)
  
  df_training <- data.frame(
    "x" = locs[, 1],
    "y" = locs[, 2],
    "cov_x" = v_cv(locs[, 1])
  )
  
  list_formulas <- list(
    "mean" = 0,
    "std.dev" = as.formula(" ~ 1"),
    "scale" = as.formula(" ~ 1 + cov_x"),
    "aniso" = 0,
    "tilt" = 0,
    "smooth" = 2.5,
    "nugget" = -Inf
  )
  
  CoCo_object <- coco(
    type = "dense",
    model.list = list_formulas,
    locs = as.matrix(locs),
    z = matrix(numeric(length = dim(df_training)[1]), ncol = 1),
    data = df_training
    )
  
  ####
  
  df_training_nn <- data.frame(
    "x" = locs[, 1],
    "y" = locs[, 2],
    "cov_x" = v_cv_nn(locs[, 1])
  )
  
  list_formulas <- list(
    "mean" = 0,
    "std.dev" = as.formula(" ~ 1"),
    "scale" = as.formula(" ~ 1 + cov_x"),
    "aniso" = 0,
    "tilt" = 0,
    "smooth" = 2.5,
    "nugget" = -Inf
  )
  
  CoCo_object_nn <- coco(
    type = "dense",
    model.list = list_formulas,
    locs = as.matrix(locs),
    z = numeric(length = dim(df_training)[1]),
    data = df_training_nn
  )
  
  no_noise_sim <- cocoSim(CoCo_object_nn,
                              pars = c(0, log(0.25), 1), seed = 181222 * 8,
                              standardize = FALSE
  )
  
  noise_sim <- cocoSim(CoCo_object,
                           pars = c(0, log(0.25), 1),
                           seed = 181222 * 8,
                           standardize = FALSE
  )
  
  pdf("Figures/noisy_covs_a.pdf", width = master_width, height = master_height)
  par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
  layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
  plot(df_training_nn$x, df_training_nn$cov_x,
       type = "l", col = "blue", ylim = c(-3.3, 3.3),
       cex = 0.3, xlab = expression(s), ylab = expression(y), pch = 20, cex.lab = 1.5, cex.axis = 1.45, lwd = 2
  )
  lines(df_training_nn$x, no_noise_sim, col = "red", lwd = 2)
  axis(4, col = "blue", cex.axis = 1.45)
  dev.off()
  
  ## another one
  
  locs_x <- seq(-10, 10, by = 0.05)
  locs <- cbind(locs_x, rep(0, length(locs_x)))
  
  create_covariate <- function(x) {
    if (x < 0) {
      0
    } else {
      1
    }
  }
  
  v_cv <- Vectorize(create_covariate)
  
  df_training <- data.frame(
    "x" = locs[, 1],
    "y" = locs[, 2],
    "cov_x" = v_cv(locs[, 1])
  )
  
  list_formulas <- list(
    "mean" = 0,
    "std.dev" = as.formula(" ~ 1"),
    "scale" = as.formula(" ~ 1 + cov_x"),
    "aniso" = 0,
    "tilt" = 0,
    "smooth" = 2.5,
    "nugget" = -Inf
  )
  
  CoCo_object <- coco(
    type = "dense",
    model.list = list_formulas,
    locs = as.matrix(locs),
    z = matrix(numeric(length = dim(df_training)[1]), ncol = 1),
    data = df_training
  )
  
  one_jump <- cocoSim(CoCo_object,
                           pars = c(0, log(0.25), 3),
                           seed = 181222 * 2,
                           standardize = FALSE
  )
  
  create_covariate_many_jumps <- function(x) {
    if (round(x, 0) %% 2 > 0) {
      1
    } else {
      0
    }
  }
  
  ####
  
  df_training_mc <- data.frame(
    "x" = locs[, 1],
    "y" = locs[, 2],
    "cov_x" = Vectorize(create_covariate_many_jumps)(locs[, 1])
  )
  
  list_formulas <- list(
    "mean" = 0,
    "std.dev" = as.formula(" ~ 1"),
    "scale" = as.formula(" ~ 1 + cov_x"),
    "aniso" = 0,
    "tilt" = 0,
    "smooth" = 2.5,
    "nugget" = -Inf
  )
  
  CoCo_object_nn <- coco(
    type = "dense",
    model.list = list_formulas,
    locs = as.matrix(locs),
    z = numeric(length = dim(df_training)[1]),
    data = df_training_mc
  )
  
  many_jumps_nn <- cocoSim(CoCo_object_nn,
                              pars = c(0, log(0.25), 3), seed = 181222 * 2,
                              standardize = FALSE
  )
  
  pdf("Figures/categorical_covs_a.pdf", width = master_width, height = master_height)
  par(mfrow = c(1, 2), oma = master_oma, mar = master_mars)
  layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
  plot(df_training$x, df_training$cov_x,
       type = "l", col = "blue", ylim = c(-3.3, 3.3),
       cex = 0.3, xlab = expression(s), ylab = expression(y), pch = 20, cex.lab = 1.5, cex.axis = 1.45, lwd = 2
  )
  axis(4, cex.lab = 1.5, cex.axis = 1.45, col = "blue")
  lines(df_training$x, one_jump, col = "red", lwd = 2)
  dev.off()
  
  pdf("Figures/categorical_covs_b.pdf", width = master_width, height = master_height)
  par(mfrow = c(1, 2), oma = master_oma, mar = master_mars)
  layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
  plot(df_training_mc$x, df_training_mc$cov_x,
       type = "l", col = "blue", ylim = c(-3.3, 3.3),
       cex = 0.3, xlab = expression(s), ylab = expression(y), pch = 20, cex.lab = 1.5, cex.axis = 1.45, lwd = 2
  )
  lines(df_training_mc$x, many_jumps_nn, col = "red", lwd = 2)
  axis(4, cex.lab = 1.5, cex.axis = 1.45, col = "blue")
  dev.off()
}

# Plot regularization example plots

if (T) {
  
  sample_sizes <- 250 + 400
  
  create_covariate_one <- Vectorize(function(x) {
    return(sin(x * 10))
    if (x < 0) {
      0
    } else {
      1
    }
  })
  
  seq_x <- seq(0, 1, length.out = floor(sqrt(750)))
  
  locs <- expand.grid(seq_x, seq_x)
  
  df_training_mc <- data.frame(
    "x" = locs[, 1],
    "y" = locs[, 2],
    "cov_x" = create_covariate_one(locs[, 1]),
    "cov_y" = create_covariate_one(locs[, 2]),
    "cov_xy" = create_covariate_one(locs[, 1] * locs[, 2])
  )
  
  list_formulas <- list(
    "mean" = as.formula(" ~ 1 + cov_x + cov_y"),
    "std.dev" = as.formula(" ~ 1 + cov_x + cov_y"),
    "scale" = as.formula(" ~ 1  + cov_x + cov_y"),
    "aniso" = 0,
    "tilt" = 0,
    "smooth" = 1.5,
    "nugget" = -Inf
  )
  
  CoCo_object_nn <- coco(
    type = "dense",
    model.list = list_formulas,
    locs = as.matrix(locs),
    z = numeric(length = dim(df_training_mc)[1]),
    data = df_training_mc
  )
  
  intensity_effect <- seq(0, 5, length.out = 30)
  
  list_results <- list()
  
  true_pars <- c(
    c(0, 0.25, -0.25) * 2,
    0, 0.5, 0.25,
    log(0.25), -0.25, 0.2
  )
  
  regula_sim <- cocoSim(
    coco.object = CoCo_object_nn,
    n = 1,
    type = "classic",
    pars = true_pars, seed = 181222,
    standardize = TRUE
  )
  
  lambda_seq <- seq(0,1,length.out = 50)

  # Correctly specified
  
  if(T){
    
    vector_lambda_cs <- numeric(length = length(lambda_seq))
    cond_numbers_cs <- numeric(length = length(lambda_seq))
    log_liks_cs <- numeric(length = length(lambda_seq))
    RMSPE_vector_cs <- numeric(length = length(lambda_seq))
    CRPS_vector_cs <- numeric(length = length(lambda_seq))
    matrix_parameters_cs <- matrix(ncol = 10, nrow = length(lambda_seq))
    
    for (ii in 1:length(lambda_seq)) {
      
      list_formulas <- list(
        "mean" = as.formula(" ~ 1 + cov_x + cov_y"),
        "std.dev" = as.formula(" ~ 1 + cov_x + cov_y"),
        "scale" = as.formula(" ~ 1  + cov_x + cov_y"),
        "aniso" = 0,
        "tilt" = 0,
        "smooth" = as.formula(" ~ 1"),
        "nugget" = -Inf
      )
      
      set.seed(12222)
      index_training <- sort(sample(650, 250))
      to_use <- (c(regula_sim[index_training,1]) - mean(c(regula_sim[index_training,1]))) / sd(c(regula_sim[index_training,1]))
      
      CoCo_object_ms <- coco(
        type = "dense",
        model.list = list_formulas,
        locs = as.matrix(locs[index_training, ]),
        z = to_use,
        data = df_training_mc[index_training, ],
        info = list(
          "lambda.reg" = lambda_seq[ii],
          "smooth.limits" = c(0.5, 3)
        )
      )
      
      tmp_boundaries <- getBoundaries(CoCo_object_ms, -20, 20)
      
      Optim_new <- cocoOptim(
        coco.object = CoCo_object_ms,
        ncores = ncores,
        boundaries = tmp_boundaries,
        optim.type = "ml"
      )
      
      cond_numbers_cs[ii] <- kappa(cov2cor(getCovMatrix(Optim_new)),exact = TRUE)
      log_liks_cs[ii] <- getLoglik(Optim_new)
      
      preds_aa <- cocoPredict(Optim_new,
                              newdataset = df_training_mc[-index_training, ],
                              newlocs = as.matrix(locs[-index_training, ]),
                              type = "pred"
      )
      
      to_use_pred <- (c(regula_sim[-index_training]) - mean(c(regula_sim[index_training]))) / sd(c(regula_sim[index_training]))
      
      RMSPE_vector_cs[ii] <- sqrt(mean((to_use_pred - (preds_aa$systematic + preds_aa$stochastic))^2))
      
      CRPS_vector_cs[ii] <- mean(getCRPS(to_use_pred, mean.pred = preds_aa$systematic + preds_aa$stochastic, sd.pred = preds_aa$sd.pred))
      
      matrix_parameters_cs[ii, ] <- unlist(getEstims(Optim_new))[c(1:3,4:9,16)]
      
    }
    
    # Prediction metrics under well specified model
    
    if(T){
      pdf("Figures/regu_pred.pdf", width = master_width_large, height = master_height)
      master_mars_2 <- master_mars
      master_mars_2[2] <- master_mars_2[2] + 1.9
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars_2)
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      plot(lambda_seq, CRPS_vector_cs / CRPS_vector_cs[1], type = 'l', 
           xlab = expression(lambda), ylab = expression("Ratio wrt. " * hat(bold(vartheta))[ML]), lwd = 2, cex.lab = 1.5, cex.axis = 1.5,
           xaxt = "n", yaxt = "n",xaxs = "i", yaxs = "i", ylim = c(0.94, 1.06))
      abline(h = 1, lty = 3, lwd = 2, col = 'gray90')
      lines(lambda_seq, RMSPE_vector_cs / RMSPE_vector_cs[1], col = 'red',lwd = 2)
      legend(bty = "n",'topleft',legend = c("RMSPE", "CRPS"),
             lty = c(1, 1), 
             col = c('black', 'red'), ncol = 1, cex = 1.25
      )
      axis(side = 1, cex.axis = 1.5)
      axis(side = 2, cex.axis = 1.5)
      dev.off()
    }
    
    # Condition Number under well specified model
    if(T){
      
      pdf("Figures/regu_cn.pdf", width = master_width_large, height = master_height)
      master_mars_2 <- master_mars
      master_mars_2[2] <- master_mars_2[2] + 1.9
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars_2)
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      plot(lambda_seq, cond_numbers_cs / cond_numbers_cs[1], type = 'l', 
           xlab = expression(lambda), ylab = expression("Ratio wrt. " * hat(bold(vartheta))[ML]), lwd = 2, cex.lab = 1.5, cex.axis = 1.5,
           xaxt = "n", yaxt = "n",xaxs = "i",yaxs = "i",ylim = c(0,1.05))
      axis(side = 1, cex.axis = 1.5)
      axis(side = 2, cex.axis = 1.5)
      abline(h = 1, lty = 3, lwd = 2, col = 'gray90')
      dev.off()
    }
    
    # Relative bias under well specified model
    if(T){
      
      pdf("Figures/regu_bias.pdf", width = master_width_large, height = master_height)
      master_mars_2 <- master_mars
      master_mars_2[2] <- master_mars_2[2] + 1.9
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars_2)
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      matplot_plot <- apply(matrix_parameters_cs, MARGIN = 1, function(x) x - matrix_parameters_cs[1, ])
      matplot(x = lambda_seq, y = t(matplot_plot), type = 'l', lty = c(1,1,1,1,1,1,1,1,1,1),
              col = c(1,'gray90', 'gray90', 2, 'gray90', 'gray90', 3, 'gray90', 'gray90', 4), xaxs = "i", yaxs = "i",
              ylab =expression("Change wrt. " * hat(bold(vartheta))[ML]), 
              xlab = expression(lambda), lwd = 2, ylim = c(-1.5,1.2), cex.lab = 1.5, cex.axis = 1.5)
      legend(bty = "n",'topleft', legend = c(expression(beta[0]), 
                                   expression(alpha[1]), 
                                   expression(theta[ms * "," * 1]), 
                                   expression(xi[1]), 
                                   's.v. coefs.'), lty = c(1, 1, 1, 1, 1),
             col = c(1, 2, 3, 4, 'gray90'), cex = 1.25, lwd = 2, ncol = 3)
      dev.off()
      master_mars[2] <- master_mars[2] - 0.6
    }
    
}
  
}

# Plot Covariates

if (T) {
  
  df_app <- all_dfs[[4]]
  
  if (T) {
    if (T) {
      
      png("Figures/illu_prec.png",
          width = master_width * 3, height = master_height * 2, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      
      tp_prec <- df_app$prec
      
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(tp_prec,
                                     breaks = seq(min(tp_prec), max(tp_prec), length.out = 128), labels = F
           )],
           cex = 0.35, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1
      )
      map(add = TRUE, resolution = 0,lwd=2)
      rect(xleft = 7.029167 - 0.35, ybottom = 46.179167 - 0.35, xright = 7.029167 + 0.35, ytop = 46.179167 + 0.35, border = "blue", lwd = 2, lty = 1)
      text(x = 7.029167 - 0.3, y = 46.179167 - 0.275, labels = "A", cex = 2)
      rect(xleft = 9.029167 - 0.35, ybottom = 46.920833 - 0.35, xright = 9.029167 + 0.35, ytop = 46.920833 + 0.35, border = "blue", lwd = 2, lty = 1)
      text(x = 9.029167 - 0.3, y = 46.920833 - 0.275, labels = "B", cex = 2)
      
      rect(xleft = 7.354167 - 0.35, ybottom = 47.10417 - 0.35, xright = 7.354167 + 0.35, ytop = 47.10417 + 0.35, border = "blue", lwd = 2, lty = 1)
      text(x = 7.354167 - 0.3, y = 47.10417 - 0.275, labels = "C", cex = 2)
      
      # ref locations
      
      points(x = 7.029167, y = 46.179167, pch=23,col='black',cex=2,bg='white')
      points(x = 9.029167, y = 46.920833, pch=23,col='black',cex=2,bg='white')
      points(x = 7.354167, y = 47.10417, pch=23,col='black',cex=2,bg='white')
      
      master_mars_second_2 <- master_mars_second
      master_mars_second_2[4] <- 7
      par(mar = master_mars_second_2)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(min(tp_prec), max(tp_prec), len = 100)
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 2)
      dev.off()
    }
    
    if (T) {
      png("Figures/illu_wind.png",
          width = master_width_large,
          height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      tp_wind <- df_app$wind
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(tp_wind,
                                     breaks = seq(min(tp_wind), max(tp_wind), length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.lab = 1.5, cex.axis = 1.5
      )
      map(add = TRUE, resolution = 0)
      par(mar = master_mars_second)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(min(tp_wind), max(tp_wind), len = 100) # supposing 3 and 2345 are the range of your data
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 1.5)
      dev.off()
    }
    
    if (T) {
      png("Figures/illu_cloud.png",
          width = master_width_large,
          height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      tp_wind <- df_app$cloud_c
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(tp_wind,
                                     breaks = seq(min(tp_wind), max(tp_wind), length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.lab = 1.5, cex.axis = 1.5
      )
      map(add = TRUE, resolution = 0)
      par(mar = master_mars_second)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(min(tp_wind), max(tp_wind), len = 100) # supposing 3 and 2345 are the range of your data
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 1.5)
      dev.off()
    }
    
    if (T) {
      png("Figures/illu_elev.png",
          width = master_width_large,
          height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      tp_wind <- df_app$elevation
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(tp_wind,
                                     breaks = seq(min(tp_wind), max(tp_wind), length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.lab = 1.5, cex.axis = 1.5
      )
      map(add = TRUE, resolution = 0)
      par(mar = master_mars_second)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(min(tp_wind), max(tp_wind), len = 100) # supposing 3 and 2345 are the range of your data
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 1.5)
      dev.off()
    }
    
    if (T) {
      png("Figures/illu_new_wind.png",
          width = master_width_large,
          height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      tp_wind <- df_app$new_wind
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(tp_wind,
                                     breaks = seq(min(tp_wind), max(tp_wind), length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.lab = 1.5, cex.axis = 1.5
      )
      map(add = TRUE, resolution = 0)
      par(mar = master_mars_second)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(min(tp_wind), max(tp_wind), len = 100) # supposing 3 and 2345 are the range of your data
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 1.5)
      dev.off()
    }
    
    if (T) {
      png("Figures/illu_BIO04.png",
          width = master_width_large,
          height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      tp_wind <- df_app$BIO04
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(tp_wind,
                                     breaks = seq(min(tp_wind), max(tp_wind), length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.lab = 1.5, cex.axis = 1.5
      )
      map(add = TRUE, resolution = 0)
      par(mar = master_mars_second)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(min(tp_wind), max(tp_wind), len = 100) # supposing 3 and 2345 are the range of your data
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 1.5)
      dev.off()
    }
    
    if (T) {
      png("Figures/illu_BIO15.png",
          width = master_width_large,
          height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      tp_wind <- df_app$BIO15
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(tp_wind,
                                     breaks = seq(min(tp_wind), max(tp_wind), length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.lab = 1.5, cex.axis = 1.5
      )
      map(add = TRUE, resolution = 0)
      par(mar = master_mars_second)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(min(tp_wind), max(tp_wind), len = 100) # supposing 3 and 2345 are the range of your data
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 1.5)
      dev.off()
    }
  }
}

# Plot Models

# Reduced dense
if (T) {
  
  load("RData/Model_B.RData")
  
  x <- Model_A
  
  x@data <- all_dfs[[4]]
  x@locs <- as.matrix(all_dfs[[4]][, 1:2])
  
  spat_effects_c <- getSpatEffects(x)
  
  # Sd
  if (T) {
    png(
      file = "Figures/reduced_model_se.png",
      width = master_width_large,
      height = master_height_large,
      units = "in",
      res = master_res
    )
    
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
     
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    
    plot(x@locs[, 1], x@locs[, 2],
         col = tim.colors(128)[cut(spat_effects_c$sd,
                                   breaks = seq(min(spat_effects_c$sd), max(spat_effects_c$sd), length.out = 128), labels = F
         )],
         cex = 0.2, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
         asp = 1.5
    )
    map(add = TRUE, resolution = 0)
    par(mar = master_mars_second)
    
    tmp_z <- matrix(1:100, nrow = 1)
    tmp_x <- 1
    tmp_y <- seq(min(spat_effects_c$sd), max(spat_effects_c$sd), len = 100)
    image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
    axis(4, lwd = 2, cex.axis = 1.5)
    dev.off()
  }
  
  # Scale
  if (T) {
    png(
      file = "Figures/reduced_model_scale.png",
      width = master_width_large,
      height = master_height_large,
      units = "in",
      res = master_res
    )
    
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
     
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    
    plot(x@locs[, 1], x@locs[, 2],
         col = tim.colors(128)[cut(spat_effects_c$scale_x,
                                   breaks = seq(min(spat_effects_c$scale_x), max(spat_effects_c$scale_x), length.out = 128), labels = F
         )],
         cex = 0.2, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
         asp = 1.5
    )
    map(add = TRUE, resolution = 0)
    par(mar = master_mars_second)
    
    tmp_z <- matrix(1:100, nrow = 1)
    tmp_x <- 1
    tmp_y <- seq(min(spat_effects_c$scale_x), max(spat_effects_c$scale_x), len = 100)
    image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
    axis(4, lwd = 2, cex.axis = 1.5)
    dev.off()
  }
  
  # Mean
  
  if (T) {
    
    png(
      file = "Figures/reduced_model_trend.png",
      width = master_width_large,
      height = master_height_large,
      units = "in",
      res = master_res
    )
    
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
    
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    
    spat_mean <- getSpatMean(x)
    
    plot(x@locs[, 1], x@locs[, 2],
         col = tim.colors(128)[cut(spat_mean,
                                   breaks = seq(min(spat_mean), max(spat_mean), length.out = 128), labels = F
         )],
         cex = 0.2, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
         asp = 1.5
    )
    map(add = TRUE, resolution = 0)
    par(mar = master_mars_second)
    
    tmp_z <- matrix(1:100, nrow = 1)
    tmp_x <- 1
    tmp_y <- seq(min(spat_mean), max(spat_mean), len = 100)
    image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
    axis(4, lwd = 2, cex.axis = 1.5)
    dev.off()
    
  }
  
}

# Taper
if (T) {
  
  load("RData/Model_T_A.RData")
  
  x <- Model_T_A
  
  x@data <- all_dfs[[4]]
  x@locs <- as.matrix(all_dfs[[4]][, 1:2])
  x@type <- "sparse"
  names(x@model.list)[2] <- "std.dev"
  
  spat_effects_c <- getSpatEffects(x)
  
  # se
  if (T) {
    png(
      file = "Figures/taper_se.png",
      width = master_width_large,
      height = master_height_large,
      units = "in",
      res = master_res
    )
    
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
     
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    
    plot(x@locs[, 1], x@locs[, 2],
         col = tim.colors(128)[cut(spat_effects_c$sd,
                                   breaks = seq(min(spat_effects_c$sd), max(spat_effects_c$sd), length.out = 128), labels = F
         )],
         cex = 0.2, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
         asp = 1.5
    )
    map(add = TRUE, resolution = 0)
    par(mar = master_mars_second)
    
    tmp_z <- matrix(1:100, nrow = 1)
    tmp_x <- 1
    tmp_y <- seq(min(spat_effects_c$sd), max(spat_effects_c$sd), len = 100)
    image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
    axis(4, lwd = 2, cex.axis = 1.5)
    dev.off()
  }
  
  # scale
  if (T) {
    png(
      file = "Figures/taper_scale.png",
      width = master_width_large,
      height = master_height_large,
      units = "in",
      res = master_res
    )
    
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
     
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    
    plot(x@locs[, 1], x@locs[, 2],
         col = tim.colors(128)[cut(spat_effects_c$scale_x,
                                   breaks = seq(min(spat_effects_c$scale_x), max(spat_effects_c$scale_x), length.out = 128), labels = F
         )],
         cex = 0.2, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
         asp = 1.5
    )
    map(add = TRUE, resolution = 0)
    par(mar = master_mars_second)
    
    tmp_z <- matrix(1:100, nrow = 1)
    tmp_x <- 1
    tmp_y <- seq(min(spat_effects_c$scale_x), max(spat_effects_c$scale_x), len = 100)
    image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
    axis(4, lwd = 2, cex.axis = 1.5)
    dev.off()
  }
  
  if (T) {
    
    png(
      file = "Figures/taper_trend.png",
      width = master_width_large,
      height = master_height_large,
      units = "in",
      res = master_res
    )
    
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
    
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    
    spat_mean <- getSpatMean(x)
    
    plot(x@locs[, 1], x@locs[, 2],
         col = tim.colors(128)[cut(spat_mean,
                                   breaks = seq(min(spat_mean), max(spat_mean), length.out = 128), labels = F
         )],
         cex = 0.2, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
         asp = 1.5
    )
    map(add = TRUE, resolution = 0)
    par(mar = master_mars_second)
    
    tmp_z <- matrix(1:100, nrow = 1)
    tmp_x <- 1
    tmp_y <- seq(min(spat_mean), max(spat_mean), len = 100)
    image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
    axis(4, lwd = 2, cex.axis = 1.5)
    dev.off()
    
  }
  
}

# Correlation plots - Dense model

if (T) {
  
  # A
  if (T) {
    load("RData/Model_A.RData")
    
    load("RData/Model_B.RData")
    
    x <- Model_A
    
    ref_loc <- c(7.029167, 46.179167)
    
    subset_plot <- which(all_dfs[[4]][, "long"] < ref_loc[1] + 0.35 & all_dfs[[4]][, "long"] > ref_loc[1] - 0.35 &
                           all_dfs[[4]][, "lati"] < ref_loc[2] + 0.3 & all_dfs[[4]][, "lati"] > ref_loc[2] - 0.3)
    
    x@data <- all_dfs[[4]][subset_plot, ]
    x@locs <- as.matrix(all_dfs[[4]][subset_plot, c('long','lati')])
    
    index_loc <- which.min(as.matrix(nearest.dist(x@locs, matrix(ref_loc, ncol = 2))))
    
    tmp_info <- cocons::getDesignMatrix(
      model.list = x@model.list,
      data = x@data
    )
    
    testt <- getCovMatrix(x)
    
    cor_testt <- cov2cor(testt)
    
    y <- Model_B
    
    subset_plot <- which(all_dfs[[4]][, "long"] < ref_loc[1] + 0.35 & all_dfs[[4]][, "long"] > ref_loc[1] - 0.35 &
                           all_dfs[[4]][, "lati"] < ref_loc[2] + 0.3 & all_dfs[[4]][, "lati"] > ref_loc[2] - 0.3)
    
    y@data <- all_dfs[[4]][subset_plot, ]
    
    y@locs <- as.matrix(all_dfs[[4]][subset_plot, c("long","lati")])
    
    tmp_info <- cocons::getDesignMatrix(
      model.list = y@model.list,
      data = y@data
    )
    
    testt_y <- getCovMatrix(y)
    
    cor_testt_y <- cov2cor(testt_y)
    
    rm(testt_y)
    
    #iso_info <- quilt.plot(x@locs, cor_testt_y[index_loc, ], nx = 84, ny = 72, plot = FALSE)
    #contour(iso_info$x, iso_info$y, iso_info$z, add = TRUE, col = "white", levels = c(0.99, 0.95, 0.9, 0.85, 0.8), lty = 1)
    
    png(
      file = "Figures/corr_dense.png",
      width = master_width_large,
      height = master_height_large,
      units = "in",
      res = master_res * 2
    )
    
    par(mfrow = c(1, 1), oma = master_oma, mar = c(4.05, 4, 0.5, 0.5))
    quilt.plot(x@locs, cor_testt[index_loc, ], nx = 84, ny = 72, zlim = c(0.8, 1), xlab = "Longitude (°)", ylab = "Latitude (°)")
    text(x = min(x@locs[, 1]) + 0.05, y = min(x@locs[, 2]) + 0.05, labels = "A", cex = 2)
    cor_info <- quilt.plot(x@locs, cor_testt[index_loc, ], nx = 84, ny = 72, zlim = c(0.8, 1), plot = FALSE)
    elev_info <- quilt.plot(x@locs, x@data$elevation, nx = 84, ny = 72, plot = FALSE)
    contour(elev_info$x, elev_info$y, elev_info$z,nlevels = 6 , add = TRUE)
    points(x@locs[index_loc, 1], x@locs[index_loc, 2], col = "white", pch = 18,cex=1.0)
    map(add = TRUE, resolution = 0)
    iso_info <- quilt.plot(x@locs, cor_testt_y[index_loc, ], nx = 84, ny = 72, plot = FALSE)
    contour(iso_info$x, iso_info$y, iso_info$z, add = TRUE, col = "white", levels = c(0.99, 0.95, 0.9, 0.85, 0.8), lty = 1,cex=1.5)
    dev.off()
  }
  
  # B
  
  if (T) {
    load("RData/Model_A.RData")
    
    load("RData/Model_B.RData")
    
    x <- Model_A
    
    ref_loc <- c(9.029167, 46.920833)
    
    subset_plot <- which(all_dfs[[4]][, "long"] < ref_loc[1] + 0.35 & all_dfs[[4]][, "long"] > ref_loc[1] - 0.35 & all_dfs[[4]][, "lati"] < ref_loc[2] + 0.3 & all_dfs[[4]][, "lati"] > ref_loc[2] - 0.3)
    
    # subset_plot <- which(all_dfs[[4]][,1] < 9.2 & all_dfs[[4]][,1] > 8.6 & all_dfs[[4]][,2] < 47.5 & all_dfs[[4]][,2] > 46,75)
    
    x@data <- all_dfs[[4]][subset_plot, ]
    x@locs <- as.matrix(all_dfs[[4]][subset_plot, c("long","lati")])
    
    index_loc <- which.min(as.matrix(nearest.dist(x@locs, matrix(ref_loc, ncol = 2))))
    
    tmp_info <- cocons::getDesignMatrix(
      model.list = x@model.list,
      data = x@data
    )
    
    testt <- getCovMatrix(x)
    
    cor_testt <- cov2cor(testt)
    
    index_loc_2 <- 320
    
    # Classic
    
    y <- Model_B
    
    subset_plot <- which(all_dfs[[4]][, "long"] < ref_loc[1] + 0.35 & all_dfs[[4]][, "long"] > ref_loc[1] - 0.35 & all_dfs[[4]][, "lati"] < ref_loc[2] + 0.3 & all_dfs[[4]][, "lati"] > ref_loc[2] - 0.3)
    
    y@data <- all_dfs[[4]][subset_plot, ]
    y@locs <- as.matrix(all_dfs[[4]][subset_plot, c("long","lati")])
    
    tmp_info <- cocons::getDesignMatrix(
      model.list = y@model.list,
      data = y@data
    )
    
    testt_y <- getCovMatrix(y)
    
    cor_testt_y <- cov2cor(testt_y)
    
    rm(testt_y)
    
    iso_info <- quilt.plot(x@locs, cor_testt_y[index_loc, ], nx = 84, ny = 72, plot = FALSE)
    # contour(iso_info$x, iso_info$y, iso_info$z, add = TRUE, col = "white", levels = c(0.99, 0.95, 0.9, 0.85, 0.8), lty = 1)
    
    png(
      file = "Figures/corr_dense_two.png",
      width = master_width_large,
      height = master_height_large,
      units = "in",
      res = master_res * 2
    )
    
    par(mfrow = c(1, 1), oma = master_oma, mar = c(4.05, 4, 0.5, 0.5))
    quilt.plot(x@locs, cor_testt[index_loc, ], nx = 84, ny = 72, zlim = c(0.8, 1), xlab = "Longitude (°)", ylab = "Latitude (°)")
    text(x = min(x@locs[, 1]) + 0.05, y = min(x@locs[, 2]) + 0.05, labels = "B", cex = 2)
    cor_info <- quilt.plot(x@locs, cor_testt[index_loc, ], nx = 84, ny = 72, zlim = c(0.8, 1), plot = FALSE)
    elev_info <- quilt.plot(x@locs, x@data$cloud_c, nx = 84, ny = 72, plot = FALSE)
    contour(elev_info$x, elev_info$y, elev_info$z, add = TRUE)
    points(x@locs[index_loc, 1], x@locs[index_loc, 2], col = "white", pch = 18,cex=1.0)
    iso_info <- quilt.plot(x@locs, cor_testt_y[index_loc, ], nx = 84, ny = 72, plot = FALSE)
    contour(iso_info$x, iso_info$y, iso_info$z, add = TRUE, col = "white", levels = c(0.99, 0.95, 0.9, 0.85, 0.8, 0.75), lty = 1,cex=1.5)
    dev.off()
  }
  
  # C
  if (T) {
    
    load("RData/Model_A.RData")
    
    load("RData/Model_B.RData")
    
    x <- Model_A
    
    ref_loc <- c(7.354167, 47.10417)
    
    subset_plot <- which(all_dfs[[4]][, "long"] < ref_loc[1] + 0.35 & all_dfs[[4]][, "long"] > ref_loc[1] - 0.35 & all_dfs[[4]][, "lati"] < ref_loc[2] + 0.3 & all_dfs[[4]][, "lati"] > ref_loc[2] - 0.3)
    
    # subset_plot <- which(all_dfs[[4]][,1] < 9.2 & all_dfs[[4]][,1] > 8.6 & all_dfs[[4]][,2] < 47.5 & all_dfs[[4]][,2] > 46,75)
    
    x@data <- all_dfs[[4]][subset_plot, ]
    x@locs <- as.matrix(all_dfs[[4]][subset_plot, c("long","lati")])
    
    index_loc <- which.min(as.matrix(nearest.dist(x@locs, matrix(ref_loc, ncol = 2))))
    
    tmp_info <- cocons::getDesignMatrix(
      model.list = x@model.list,
      data = x@data
    )
    
    testt <- getCovMatrix(x)
    
    cor_testt <- cov2cor(testt)
    
    # Classic
    
    y <- Model_D
    
    subset_plot <- which(all_dfs[[4]][, "long"] < ref_loc[1] + 0.35 & all_dfs[[4]][, "long"] > ref_loc[1] - 0.35 & all_dfs[[4]][, "lati"] < ref_loc[2] + 0.3 & all_dfs[[4]][, "lati"] > ref_loc[2] - 0.3)
    
    y@data <- all_dfs[[4]][subset_plot, ]
    y@locs <- as.matrix(all_dfs[[4]][subset_plot, c("long","lati")])
    
    tmp_info <- cocons::getDesignMatrix(
      model.list = y@model.list,
      data = y@data
    )
    
    testt_y <- getCovMatrix(y)
    
    cor_testt_y <- cov2cor(testt_y)
    
    rm(testt_y)
    
    iso_info <- quilt.plot(y@locs, cor_testt_y[index_loc, ], nx = 84, ny = 72, plot = FALSE)
    
    # contour(iso_info$y,iso_info$y, iso_info$z,add=TRUE,col='white',levels = c(0.99,0.95,0.9,0.85,0.8),lty=1)
    
    png(
      file = "Figures/corr_dense_three.png",
      width = master_width_large,
      height = master_height_large,
      units = "in",
      res = master_res * 2
    )
    
    par(mfrow = c(1, 1), oma = master_oma, mar = c(4.05, 4, 0.5, 0.5))
    quilt.plot(x@locs, cor_testt[index_loc, ], nx = 84, ny = 72, zlim = c(0.8, 1), xlab = "Longitude (°)", ylab = "Latitude (°)")
    text(x = min(x@locs[, 1]) + 0.05, y = min(x@locs[, 2]) + 0.05, labels = "C", cex = 2)
    cor_info <- quilt.plot(x@locs, cor_testt[index_loc, ], nx = 84, ny = 72, zlim = c(0.8, 1), plot = FALSE)
    elev_info <- quilt.plot(x@locs, x@data$wind, nx = 84, ny = 72, plot = FALSE)
    contour(elev_info$x, elev_info$y, elev_info$z, add = TRUE)
    points(x@locs[index_loc, 1], x@locs[index_loc, 2], col = "white", pch = 18,cex=1.0)
    iso_info <- quilt.plot(y@locs, cor_testt_y[index_loc, ], nx = 84, ny = 72, plot = FALSE)
    contour(iso_info$x, iso_info$y, iso_info$z, add = TRUE, col = "white", levels = c(0.99, 0.95, 0.9, 0.85, 0.8, 0.75), lty = 1,cex=1.5)
    dev.off()
  }
}

# Correlation plots - Sparse model

  if(T){
    
    # A
    if (T) {
      
      load("RData/Model_T_A.RData")
      
      load("RData/Model_T_B.RData")
      
      x <- Model_T_A
      
      ref_loc <- c(7.029167, 46.179167)
      
      subset_plot <- which(all_dfs[[4]][, "long"] < ref_loc[1] + 0.35 & all_dfs[[4]][, "long"] > ref_loc[1] - 0.35 & all_dfs[[4]][, "lati"] < ref_loc[2] + 0.3 & all_dfs[[4]][, "lati"] > ref_loc[2] - 0.3)
      
      x@data <- all_dfs[[4]][subset_plot, ]
      x@locs <- as.matrix(all_dfs[[4]][subset_plot, c("long","lati")])
      
      index_loc <- which.min(as.matrix(nearest.dist(x@locs, matrix(ref_loc, ncol = 2))))
      
      tmp_info <- cocons::getDesignMatrix(
        model.list = x@model.list,
        data = x@data
      )
      
      testt <- getCovMatrix(x)
      
      cor_testt <- cov2cor(as.matrix(testt))
      
      # Classic
      
      y <- Model_T_B
      
      subset_plot <- which(all_dfs[[4]][, "long"] < ref_loc[1] + 0.35 & all_dfs[[4]][, "long"] > ref_loc[1] - 0.35 & all_dfs[[4]][, "lati"] < ref_loc[2] + 0.3 & all_dfs[[4]][, "lati"] > ref_loc[2] - 0.3)
      
      y@data <- all_dfs[[4]][subset_plot, ]
      y@locs <- as.matrix(all_dfs[[4]][subset_plot, c("long","lati")])
      
      tmp_info <- cocons::getDesignMatrix(
        model.list = y@model.list,
        data = y@data
      )
      
      testt_y <- getCovMatrix(y)
      
      cor_testt_y <- cov2cor(as.matrix(testt_y))
      
      rm(testt_y)
      
      iso_info <- quilt.plot(x@locs, cor_testt_y[index_loc, ], nx = 84, ny = 72, plot = FALSE)
      
      png(
        file = "Figures/T_corr_dense.png",
        width = master_width_large,
        height = master_height_large,
        units = "in",
        res = master_res * 2
      )
      
      par(mfrow = c(1, 1), oma = master_oma, mar = c(4.05, 4, 0.5, 0.5))
      quilt.plot(x@locs, cor_testt[index_loc, ], nx = 84, ny = 72, zlim = c(0.05, 1), xlab = "Longitude (°)", ylab = "Latitude (°)")
      text(x = min(x@locs[, 1]) + 0.05, y = min(x@locs[, 2]) + 0.05, labels = "A", cex = 2)
      cor_info <- quilt.plot(x@locs, cor_testt[index_loc, ], nx = 84, ny = 72, zlim = c(0.05, 1), plot = FALSE)
      elev_info <- quilt.plot(x@locs, x@data$elevation, nx = 84, ny = 72, plot = FALSE)
      contour(elev_info$x, elev_info$y, elev_info$z, add = TRUE)
      map(add = TRUE, resolution = 0)
      
      points(x@locs[index_loc, 1], x@locs[index_loc, 2], col = "white", pch = 18,cex=1.0)
      iso_info <- quilt.plot(x@locs, cor_testt_y[index_loc, ], nx = 84, ny = 72, plot = FALSE)
      contour(iso_info$x, iso_info$y, iso_info$z, add = TRUE, col = "white", levels = c(0.9, 0.7, 0.4, 0.2), lty = 1)
      dev.off()
    }
    
    # B
    if (T) {
      
      load("RData/Model_T_A.RData")
      
      load("RData/Model_T_B.RData")
      
      x <- Model_T_A
      
      ref_loc <- c(9.029167, 46.920833)
      
      subset_plot <- which(all_dfs[[4]][, "long"] < ref_loc[1] + 0.35 & all_dfs[[4]][, "long"] > ref_loc[1] - 0.35 & all_dfs[[4]][, "lati"] < ref_loc[2] + 0.3 & all_dfs[[4]][, "lati"] > ref_loc[2] - 0.3)
      
      x@data <- all_dfs[[4]][subset_plot, ]
      x@locs <- as.matrix(all_dfs[[4]][subset_plot, c("long","lati")])
      
      index_loc <- which.min(as.matrix(nearest.dist(x@locs, matrix(ref_loc, ncol = 2))))
      
      tmp_info <- cocons::getDesignMatrix(
        model.list = x@model.list,
        data = x@data
      )
      
      testt <- getCovMatrix(x)
      
      cor_testt <- cov2cor(as.matrix(testt))
      
      # Classic
      
      y <- Model_T_B
      
      subset_plot <- which(all_dfs[[4]][, "long"] < ref_loc[1] + 0.35 & all_dfs[[4]][, "long"] > ref_loc[1] - 0.35 & all_dfs[[4]][, "lati"] < ref_loc[2] + 0.3 & all_dfs[[4]][, "lati"] > ref_loc[2] - 0.3)
      
      y@data <- all_dfs[[4]][subset_plot, ]
      y@locs <- as.matrix(all_dfs[[4]][subset_plot, c("long","lati")])
      
      tmp_info <- cocons::getDesignMatrix(
        model.list = y@model.list,
        data = y@data
      )
      
      testt_y <- getCovMatrix(y)
      
      cor_testt_y <- cov2cor(as.matrix(testt_y))
      
      rm(testt_y)
      
      iso_info <- quilt.plot(x@locs, cor_testt_y[index_loc, ], nx = 84, ny = 72, plot = FALSE)
      
      png(
        file = "Figures/T_corr_dense_two.png",
        width = master_width_large,
        height = master_height_large,
        units = "in",
        res = master_res * 2
      )
      
      par(mfrow = c(1, 1), oma = master_oma, mar = c(4.05, 4, 0.5, 0.5))
      quilt.plot(x@locs, cor_testt[index_loc, ], nx = 84, ny = 72, zlim = c(0.05, 1), xlab = "Longitude (°)", ylab = "Latitude (°)")
      text(x = min(x@locs[, 1]) + 0.05, y = min(x@locs[, 2]) + 0.05, labels = "B", cex = 2)
      cor_info <- quilt.plot(x@locs, cor_testt[index_loc, ], nx = 84, ny = 72, zlim = c(0.05, 1), plot = FALSE)
      elev_info <- quilt.plot(x@locs, x@data$cloud_c, nx = 84, ny = 72, plot = FALSE)
      contour(elev_info$x, elev_info$y, elev_info$z, add = TRUE)
      points(x@locs[index_loc, 1], x@locs[index_loc, 2], col = "white", pch = 18,cex=1.0)
      iso_info <- quilt.plot(x@locs, cor_testt_y[index_loc, ], nx = 84, ny = 72, plot = FALSE)
      contour(iso_info$x, iso_info$y, iso_info$z, add = TRUE, col = "white", levels = c(0.9, 0.7, 0.4, 0.2), lty = 1)
      dev.off()
    }
    
    # C
    if (T) {
      
      load("RData/Model_T_A.RData")
      load("RData/Model_T_B.RData")
      
      x <- Model_T_A
      
      ref_loc <- c(7.354167, 47.10417)
      
      subset_plot <- which(all_dfs[[4]][, "long"] < ref_loc[1] + 0.35 & all_dfs[[4]][, "long"] > ref_loc[1] - 0.35 & all_dfs[[4]][, "lati"] < ref_loc[2] + 0.3 & all_dfs[[4]][, "lati"] > ref_loc[2] - 0.3)
      
      x@data <- all_dfs[[4]][subset_plot, ]
      x@locs <- as.matrix(all_dfs[[4]][subset_plot, c("long","lati")])
      
      index_loc <- which.min(as.matrix(nearest.dist(x@locs, matrix(ref_loc, ncol = 2))))
      
      tmp_info <- cocons::getDesignMatrix(
        model.list = x@model.list,
        data = x@data
      )
      
      testt <- getCovMatrix(x)
      
      cor_testt <- cov2cor(as.matrix(testt))
      
      # Classic
      
      y <- Model_T_B
      
      subset_plot <- which(all_dfs[[4]][, "long"] < ref_loc[1] + 0.35 & all_dfs[[4]][, "long"] > ref_loc[1] - 0.35 & all_dfs[[4]][, "lati"] < ref_loc[2] + 0.3 & all_dfs[[4]][, "lati"] > ref_loc[2] - 0.3)
      
      y@data <- all_dfs[[4]][subset_plot, ]
      y@locs <- as.matrix(all_dfs[[4]][subset_plot, c("long","lati")])
      
      tmp_info <- cocons::getDesignMatrix(
        model.list = y@model.list,
        data = y@data
      )
      
      testt_y <- getCovMatrix(y)
      
      cor_testt_y <- cov2cor(as.matrix(testt_y))
      
      rm(testt_y)
      
      iso_info <- quilt.plot(x@locs, cor_testt_y[index_loc, ], nx = 84, ny = 72, plot = FALSE)
      
      png(
        file = "Figures/T_corr_dense_three.png",
        width = master_width_large,
        height = master_height_large,
        units = "in",
        res = master_res * 2
      )
      
      par(mfrow = c(1, 1), oma = master_oma, mar = c(4.05, 4, 0.5, 0.5))
      quilt.plot(x@locs, cor_testt[index_loc, ], nx = 84, ny = 72, zlim = c(0.05, 1), xlab = "Longitude (°)", ylab = "Latitude (°)")
      text(x = min(x@locs[, 1]) + 0.05, y = min(x@locs[, 2]) + 0.05, labels = "C", cex = 2)
      cor_info <- quilt.plot(x@locs, cor_testt[index_loc, ], nx = 84, ny = 72, zlim = c(0.05, 1), plot = FALSE)
      elev_info <- quilt.plot(x@locs, x@data$wind, nx = 84, ny = 72, plot = FALSE)
      contour(elev_info$x, elev_info$y, elev_info$z, add = TRUE)
      points(x@locs[index_loc, 1], x@locs[index_loc, 2], col = "white", pch = 18,cex=1.0)
      iso_info <- quilt.plot(x@locs, cor_testt_y[index_loc, ], nx = 84, ny = 72, plot = FALSE)
      contour(iso_info$x, iso_info$y, iso_info$z, add = TRUE, col = "white", levels = c(0.9, 0.7, 0.4, 0.2), lty = 1)
      dev.off()
    }
  }
# Appendix 

if (T) {
  
  # Datasets
  
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
    
    final_groups <- c(1:100)[!(c(1:100) %in% sample_to_tune_hyperparameters)]
    
    obs_groups_final <- groups[groups %in% final_groups]
    
  }  
  
  # Training dense Precipitation
  if(T){
    
    df_app <- all_dfs[[2]]
    
    png("Figures/illu_prec_dense_training.png",
        width = master_width_large, height = master_height_large, res = master_res,
        units = "in"
    )
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
    
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    
    tp_prec <- df_app$prec
    
    plot(df_app[, "long"], df_app[, "lati"],
         col = tim.colors(128)[cut(tp_prec,
                                   breaks = seq(min(tp_prec), max(tp_prec), length.out = 128), labels = F
         )],
         cex = 1, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1.5
    )
    map(add = TRUE, resolution = 0,lwd=2)
    master_mars_second_2 <- master_mars_second
    master_mars_second_2[4] <- 2
    par(mar = master_mars_second_2)
    
    tmp_z <- matrix(1:100, nrow = 1)
    tmp_x <- 1
    tmp_y <- seq(min(tp_prec), max(tp_prec), len = 100)
    image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
    axis(4, lwd = 2, cex.axis = 2)
    dev.off()
    
  }
  
  # Training Sparse Precipitation
  if(T){
    
    df_app <- all_dfs[[3]]
    
    png("Figures/illu_prec_sparse_training.png",
        width = master_width_large, height = master_height_large, res = master_res,
        units = "in"
    )
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
    
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    
    tp_prec <- df_app$prec
    
    plot(df_app[, "long"], df_app[, "lati"],
         col = tim.colors(128)[cut(tp_prec,
                                   breaks = seq(min(tp_prec), max(tp_prec), length.out = 128), labels = F
         )],
         cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1.5
    )
    map(add = TRUE, resolution = 0,lwd=2)
    master_mars_second_2 <- master_mars_second
    master_mars_second_2[4] <- 2
    par(mar = master_mars_second_2)
    
    tmp_z <- matrix(1:100, nrow = 1)
    tmp_x <- 1
    tmp_y <- seq(min(tp_prec), max(tp_prec), len = 100)
    image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
    axis(4, lwd = 2, cex.axis = 2)
    dev.off()
    
  }
  
  # Test Precipitation
  if(T){
    
    df_app <- all_dfs[[1]]
    
    df_app <- newdataset_final
    
    png("Figures/illu_prec_test.png",
        width = master_width_large, height = master_height_large, res = master_res,
        units = "in"
    )
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
    
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    
    tp_prec <- df_app$prec
    
    plot(df_app[, "long"], df_app[, "lati"],
         col = tim.colors(128)[cut(tp_prec,
                                   breaks = seq(min(tp_prec), max(tp_prec), length.out = 128), labels = F
         )],
         cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1.5
    )
    map(add = TRUE, resolution = 0,lwd=2)
    master_mars_second_2 <- master_mars_second
    master_mars_second_2[4] <- 2
    par(mar = master_mars_second_2)
    
    tmp_z <- matrix(1:100, nrow = 1)
    tmp_x <- 1
    tmp_y <- seq(min(tp_prec), max(tp_prec), len = 100)
    image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
    axis(4, lwd = 2, cex.axis = 2)
    dev.off()
    
  }
  
  # Mean predictions - Dense scenario

  if(T){
    
    load('RData/Pred_dense.RData')
    
    unique_boundaries <- range(c(Pred_A$systematic + Pred_A$stochastic,
                                 Pred_B$systematic + Pred_B$stochastic)
    )
    
    # M_STAT
    if(T){
      
      png("Figures/point_pred_dense_M_STAT.png",
          width = master_width_large, height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      
      point_pred <- Pred_B$systematic + Pred_B$stochastic
      
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(point_pred,
                                     breaks = seq(unique_boundaries[1], unique_boundaries[2], length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1.5
      )
      map(add = TRUE, resolution = 0,lwd=2)
      master_mars_second_2 <- master_mars_second
      master_mars_second_2[4] <- 2
      par(mar = master_mars_second_2)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(unique_boundaries[1], unique_boundaries[2], len = 100)
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 2)
      dev.off()
      
      
    }
    
    # M_NS
    if(T){
      
      png("Figures/point_pred_dense_M_NS.png",
          width = master_width_large, height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      
      point_pred <- Pred_A$systematic + Pred_A$stochastic
      
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(point_pred,
                                     breaks = seq(unique_boundaries[1], unique_boundaries[2], length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1.5
      )
      map(add = TRUE, resolution = 0,lwd=2)
      master_mars_second_2 <- master_mars_second
      master_mars_second_2[4] <- 2
      par(mar = master_mars_second_2)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(unique_boundaries[1], unique_boundaries[2], len = 100)
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 2)
      dev.off()
      
      
    }
    
  }
  
  # SD predictions - Dense scenario
  
  if(T){
    
    load('RData/Pred_dense.RData')
    
    unique_boundaries <- range(c(Pred_A$sd.pred,
                                 Pred_B$sd.pred)
    )
    
    # M_STAT
    if(T){
      
      png("Figures/sd_pred_dense_M_STAT.png",
          width = master_width_large, height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      
      point_pred <- Pred_B$sd.pred
      
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(point_pred,
                                     breaks = seq(unique_boundaries[1], unique_boundaries[2], length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1.5
      )
      map(add = TRUE, resolution = 0,lwd=2)
      master_mars_second_2 <- master_mars_second
      master_mars_second_2[4] <- 2
      par(mar = master_mars_second_2)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(unique_boundaries[1], unique_boundaries[2], len = 100)
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 2)
      dev.off()
      
      
    }
    
    # M_NS
    if(T){
      
      png("Figures/sd_dense_M_NS.png",
          width = master_width_large, height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      
      point_pred <- Pred_A$sd.pred
      
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(point_pred,
                                     breaks = seq(unique_boundaries[1], unique_boundaries[2], length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1.5
      )
      map(add = TRUE, resolution = 0,lwd=2)
      master_mars_second_2 <- master_mars_second
      master_mars_second_2[4] <- 2
      par(mar = master_mars_second_2)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(unique_boundaries[1], unique_boundaries[2], len = 100)
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 2)
      dev.off()
      
      
    }
    
  }
  
  # Mean predictions - Sparse scenario
  
  if(T){
    
    load('RData/Pred_Taper.RData')
    
    unique_boundaries <- range(c(Pred_T_A$systematic + Pred_T_A$stochastic,
                                 Pred_T_B$systematic + Pred_T_B$stochastic)
    )
    
    # M_STAT_T
    if(T){
      
      png("Figures/point_pred_sparse_M_STAT_T.png",
          width = master_width_large, height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      
      point_pred <- Pred_T_B$systematic + Pred_T_B$stochastic
      
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(point_pred,
                                     breaks = seq(unique_boundaries[1], unique_boundaries[2], length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1.5
      )
      map(add = TRUE, resolution = 0,lwd=2)
      master_mars_second_2 <- master_mars_second
      master_mars_second_2[4] <- 2
      par(mar = master_mars_second_2)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(unique_boundaries[1], unique_boundaries[2], len = 100)
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 2)
      dev.off()
      
      
    }
    
    # M_NS_T
    if(T){
      
      png("Figures/point_pred_sparse_M_NS_T.png",
          width = master_width_large, height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      
      point_pred <- Pred_T_A$systematic + Pred_T_A$stochastic
      
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(point_pred,
                                     breaks = seq(unique_boundaries[1], unique_boundaries[2], length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1.5
      )
      map(add = TRUE, resolution = 0,lwd=2)
      master_mars_second_2 <- master_mars_second
      master_mars_second_2[4] <- 2
      par(mar = master_mars_second_2)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(unique_boundaries[1], unique_boundaries[2], len = 100)
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 2)
      dev.off()
      
    }
    
  }
  
  # SD predictions - Dense scenario
  
  if(T){
    
    load('RData/Pred_Taper.RData')
    
    unique_boundaries <- range(c(Pred_T_A$sd.pred,
                                 Pred_T_B$sd.pred)
    )
    
    # M_STAT_T
    if(T){
      
      png("Figures/sd_pred_sparse_M_STAT_T.png",
          width = master_width_large, height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      
      point_pred <- Pred_T_B$sd.pred
      
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(point_pred,
                                     breaks = seq(unique_boundaries[1], unique_boundaries[2], length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1.5
      )
      map(add = TRUE, resolution = 0,lwd=2)
      master_mars_second_2 <- master_mars_second
      master_mars_second_2[4] <- 2
      par(mar = master_mars_second_2)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(unique_boundaries[1], unique_boundaries[2], len = 100)
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 2)
      dev.off()
      
    }
    
    # M_NS_T
    if(T){
      
      png("Figures/sd_sparse_M_NS_T.png",
          width = master_width_large, height = master_height_large, res = master_res,
          units = "in"
      )
      par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
      
      layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
      
      point_pred <- Pred_T_A$sd.pred
      
      plot(df_app[, 1], df_app[, 2],
           col = tim.colors(128)[cut(point_pred,
                                     breaks = seq(unique_boundaries[1], unique_boundaries[2], length.out = 128), labels = F
           )],
           cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1.5
      )
      map(add = TRUE, resolution = 0,lwd=2)
      master_mars_second_2 <- master_mars_second
      master_mars_second_2[4] <- 2
      par(mar = master_mars_second_2)
      
      tmp_z <- matrix(1:100, nrow = 1)
      tmp_x <- 1
      tmp_y <- seq(unique_boundaries[1], unique_boundaries[2], len = 100)
      image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
      axis(4, lwd = 2, cex.axis = 2)
      dev.off()
      
      
    }
    
  }
  
  if(T){
    
    df_app <- all_dfs[[3]]
    
    png("Figures/illu_prec_sparse_training.png",
        width = master_width_large, height = master_height_large, res = master_res,
        units = "in"
    )
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
    
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    
    tp_prec <- df_app$prec
    
    plot(df_app[, 1], df_app[, 2],
         col = tim.colors(128)[cut(tp_prec,
                                   breaks = seq(min(tp_prec), max(tp_prec), length.out = 128), labels = F
         )],
         cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1.5
    )
    map(add = TRUE, resolution = 0,lwd=2)
    master_mars_second_2 <- master_mars_second
    master_mars_second_2[4] <- 2
    par(mar = master_mars_second_2)
    
    tmp_z <- matrix(1:100, nrow = 1)
    tmp_x <- 1
    tmp_y <- seq(min(tp_prec), max(tp_prec), len = 100)
    image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
    axis(4, lwd = 2, cex.axis = 2)
    dev.off()
    
  }
  
  if(T){
    
    df_app <- newdataset_final
    
    png("Figures/illu_prec_test.png",
        width = master_width_large, height = master_height_large, res = master_res,
        units = "in"
    )
    par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
    
    layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), nrow = 1))
    
    tp_prec <- df_app$prec
    
    plot(df_app[, 1], df_app[, 2],
         col = tim.colors(128)[cut(tp_prec,
                                   breaks = seq(min(tp_prec), max(tp_prec), length.out = 128), labels = F
         )],
         cex = 0.3, xlab = "Longitude (°)", ylab = "Latitude (°)", pch = 20, cex.axis = 2, cex.lab = 2, asp = 1.5
    )
    map(add = TRUE, resolution = 0,lwd=2)
    master_mars_second_2 <- master_mars_second
    master_mars_second_2[4] <- 2
    par(mar = master_mars_second_2)
    
    tmp_z <- matrix(1:100, nrow = 1)
    tmp_x <- 1
    tmp_y <- seq(min(tp_prec), max(tp_prec), len = 100)
    image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
    axis(4, lwd = 2, cex.axis = 2)
    dev.off()
    
  }
  
}
