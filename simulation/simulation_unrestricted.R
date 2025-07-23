# Simulation code for Model 1-3, no restrictions on weights coefficients

library(tidyverse)
library(quadprog)
library(Synth)

# library(foreach)

########### sub-functions ####################

# linear regression for detrending
lsq <- function(z){
  tt = length(z)
  y = z[(h+p):tt]
  X = embed(z[1:(tt-h)],p)
  y <- ifelse(is.na(y), 0, y)
  X <- ifelse(is.na(X), 0, X)
  df = data.frame(cbind(y,X))
  model = lm(y~.,data=df)
  return(model)
}

# function: Extract residuals and form a matrix
extract_residuals <- function(lm_list) {

  residuals_list <- lapply(lm_list, function(item) {
    return(residuals(item))
  })
  
  # Combine residuals into a matrix (each column corresponds to a model)
  residuals_matrix <- do.call(cbind, residuals_list)
  return(residuals_matrix)
}

trend_predict <- function(z, bhat) {
  tt = length(z)
  predictions <- numeric(h)
  for (i in 1:h){
    lags_h <- z[(tt-h+i):(tt-h-p+i+1)]
    x <- c(1,lags_h)
    predictions[i] <- sum(x*bhat)
  }
  return(predictions)
}

run_sc_unres <- function(y_matrix, T0) {
  X1 <- as.matrix(y_matrix[1:T0, 1])
  X0 <- cbind(rep(1, T0), as.matrix(y_matrix[1:T0, -1]))
  
  D <- t(X0) %*% X0
  d <- t(X0) %*% X1
  weights.sc <- solve(D) %*% d
  
  Y1hat.SC <- cbind(rep(1, T0 + Fh), as.matrix(y_matrix[1:(T0 + Fh), -1])) %*% weights.sc
  
  errorSC_vec <- Y1hat.SC - as.matrix(y_matrix[1:(T0 + Fh), 1])
  errorSC_vec_pre <- errorSC_vec[1:T0]
  errorSC_vec_post <- errorSC_vec[(T0 + 1):(T0 + Fh)]
  
  MSE.SC_pre <- sum(errorSC_vec_pre^2) / T0
  MSE.SC_post <- sum(errorSC_vec_post^2) / Fh
  
  list(
    MSE_pre = MSE.SC_pre,
    MSE_post = MSE.SC_post
  )
}


run_sbc_unres <- function(y_matrix, T0) {
  # Step 1: Detrend control units
  list_detrended_control <- apply(y_matrix[1:(T0 + Fh), -1], 2, function(col) {
    result <- lsq(col)
    return(result)
  })
  detrended_control <- extract_residuals(list_detrended_control)
  
  # Step 2: Detrend treated unit (pre-treatment)
  detrended_treated_pre <- lsq(y_matrix[1:T0, 1])$residuals
  detrended_control_pre <- detrended_control[1:length(detrended_treated_pre), ]
  
  # Step 3: Fit synthetic control on detrended series
  X1 <- as.matrix(detrended_treated_pre)
  X0 <- as.matrix(detrended_control_pre)
  
  D <- t(X0) %*% X0
  d <- t(X0) %*% X1
  weights.sbc <- solve(D) %*% d
  
  cychat_pre <- X0 %*% weights.sbc
  cychat_post <- tail(detrended_control, Fh) %*% weights.sbc
  
  # Step 4: Estimate and extrapolate trend for treated unit
  trend_pre <- y_matrix[(h + p):T0, 1] - detrended_treated_pre
  bhat_trend <- lsq(y_matrix[1:T0, 1])$coefficients
  trend_post <- trend_predict(y_matrix[1:T0, 1], bhat_trend)
  
  # Step 5: Compute residuals and MSE
  errorSBC_vec_pre <- trend_pre + cychat_pre - y_matrix[(h + p):T0, 1]
  errorSBC_vec_post <- trend_post + cychat_post - y_matrix[(T0 + 1):(T0 + Fh), 1]
  
  MSE.SBC_pre <- sum(errorSBC_vec_pre^2) / (T0 - h - p + 1)
  MSE.SBC_post <- sum(errorSBC_vec_post^2) / Fh
  
  # Return all results
  list(
    MSE_pre = MSE.SBC_pre,
    MSE_post = MSE.SBC_post
  )
}


########### Parameters ####################

N <- 12 # number of individuals
nsim <- 10000  # number of simulations

h = 2  # detrending horizon
p = 2 # number of lags
Fh = h # forecast periods after treatment, by default set it equal to h

T0_vec <- c(50,100,200) # vector of pre-treatment periods
phi_vec <- c(.2,.5,.8) # autoregressive coef in model2
start_time <- Sys.time() 



# initialize result matrices, col 1-3 are Model 1, col 4-6 are Model 2, col 7-9 are Model 3

result_mse_pre <- matrix(nrow = 3, ncol = 9) # mse ratio pre-treatment
result_mse_post <- matrix(nrow = 3, ncol = 9) # mse ratio post-treatment

for (iT0 in 1:3) {
  
  MSE_SC_pre_mat <- matrix(nrow = nsim,ncol = 9)
  MSE_SC_post_mat <- matrix(nrow =nsim,ncol =9)
  MSE_SBC_pre_mat <- matrix(nrow =nsim,ncol =9)
  MSE_SBC_post_mat <- matrix(nrow =nsim,ncol = 9)

  for (iphi in 1:3){
    for (j in 1:nsim){
      
    
      T0 <- T0_vec[iT0]
      phi_f <- phi_vec[iphi]
      TT = T0 + Fh 
      
      # ====================== Model 1   ====================== #
      # DGP 1.1
      y_matrix <- replicate(N, {
        mu <- 0
        cumsum(rnorm(TT) + mu)
      })
      MSE_sc <- run_sc_unres(y_matrix,T0 = T0)
      MSE_sbc <- run_sbc_unres(y_matrix,T0 = T0)
      MSE_SC_pre_mat[j,1] <- MSE_sc$MSE_pre
      MSE_SC_post_mat[j,1] <- MSE_sc$MSE_post
      MSE_SBC_pre_mat[j,1] <- MSE_sbc$MSE_pre
      MSE_SBC_post_mat[j,1] <- MSE_sbc$MSE_post
      
      
      # DGP 1.2
      y_matrix <- replicate(N, {
        mu <- 0.5
        cumsum(rnorm(TT) + mu)
      })
      MSE_sc <- run_sc_unres(y_matrix,T0 = T0)
      MSE_sbc <- run_sbc_unres(y_matrix,T0 = T0)
      MSE_SC_pre_mat[j,2] <- MSE_sc$MSE_pre
      MSE_SC_post_mat[j,2] <- MSE_sc$MSE_post
      MSE_SBC_pre_mat[j,2] <- MSE_sbc$MSE_pre
      MSE_SBC_post_mat[j,2] <- MSE_sbc$MSE_post
      
      # DGP 1.3
      y_matrix <- replicate(N, {
        mu <- rnorm(1, 0, sqrt(0.25))
        cumsum(rnorm(TT) + mu)
      })
      MSE_sc <- run_sc_unres(y_matrix,T0 = T0)
      MSE_sbc <- run_sbc_unres(y_matrix,T0 = T0)
      MSE_SC_pre_mat[j,3] <- MSE_sc$MSE_pre
      MSE_SC_post_mat[j,3] <- MSE_sc$MSE_post
      MSE_SBC_pre_mat[j,3] <- MSE_sbc$MSE_pre
      MSE_SBC_post_mat[j,3] <- MSE_sbc$MSE_post
      
      # ====================== Model 2   ====================== #
      r <- 2
  
      # 1) simulate factor loadings Λ (N × r)
      Lambda <- matrix(rnorm(N * r), nrow = N, ncol = r)
      
      # 2) simulate r stationary AR(1) factors f_t
      Fmat <- matrix(0, nrow = r, ncol = TT)
      Fmat[,1] <- rnorm(r)
      for (t in 2:TT) {
        Fmat[,t] <- phi_f * Fmat[,t-1] + rnorm(r,sd=1)
      }
      
      # 3) simulate idiosyncratic errors u_it ~ iid N(0,1)
      U <- matrix(rnorm(TT * N,sd=1), nrow = TT, ncol = N)
      
      # 4) build ε_it = [Λ f_t]_i + u_it
      eps <- t(Lambda %*% Fmat) + U
      
      # 5) form the random walks y_matrix_t^(i) = y_matrix_{t-1}^(i) + μ_i + ε_it
      y_matrix <- matrix(0, nrow = TT, ncol = N)
      for (t in 1:TT) {
        if (t == 1) {
          y_matrix[t, ] <- eps[t, ]
        } else {
          y_matrix[t, ] <- y_matrix[t-1, ] + eps[t, ]
        }
      }
      
      MSE_sc <- run_sc_unres(y_matrix,T0 = T0)
      MSE_sbc <- run_sbc_unres(y_matrix,T0 = T0)
      MSE_SC_pre_mat[j,3+iphi] <- MSE_sc$MSE_pre
      MSE_SC_post_mat[j,3+iphi] <- MSE_sc$MSE_post
      MSE_SBC_pre_mat[j,3+iphi] <- MSE_sbc$MSE_pre
      MSE_SBC_post_mat[j,3+iphi] <- MSE_sbc$MSE_post
      
      # ====================== Model 3   ====================== #
      
      s <- 2
      Lambda_trend <- matrix(rnorm((.5*N)*s,sd=T0^{-1/3}), nrow=.5*N, ncol=s)
      eta_trend    <- matrix(rnorm(TT*s,sd=1), nrow=TT, ncol=s)
      f_trend  <- apply(eta_trend, 2, cumsum)
      trend_mat <- t(Lambda_trend %*% t(f_trend))
      y_matrix1  <- trend_mat + eps[,1:(.5*N)]
      
      y_matrix2 <- matrix(0, nrow = TT, ncol = .5*N)
      for (t in 1:TT) {
        if (t == 1) {
          y_matrix2[t, ] <- eps[t, (.5*N+1):N]
        } else {
          y_matrix2[t, ] <- y_matrix2[t-1, ] + eps[t, (.5*N+1):N]
        }
      }
      
      y_matrix <- cbind(y_matrix1,y_matrix2)
      MSE_sc <- run_sc_unres(y_matrix,T0 = T0)
      MSE_sbc <- run_sbc_unres(y_matrix,T0 = T0)
      MSE_SC_pre_mat[j,6+iphi] <- MSE_sc$MSE_pre
      MSE_SC_post_mat[j,6+iphi] <- MSE_sc$MSE_post
      MSE_SBC_pre_mat[j,6+iphi] <- MSE_sbc$MSE_pre
      MSE_SBC_post_mat[j,6+iphi] <- MSE_sbc$MSE_post
      

    }
  }
  
  MSE_SC_pre_mean <- colMeans(MSE_SC_pre_mat)
  MSE_SC_post_mean <- colMeans(MSE_SC_post_mat)
  MSE_SBC_pre_mean <- colMeans(MSE_SBC_pre_mat)
  MSE_SBC_post_mean <- colMeans(MSE_SBC_post_mat)
  
  MSE_pre_ratio <- MSE_SBC_pre_mean / MSE_SC_pre_mean
  MSE_post_ratio <- MSE_SBC_post_mean / MSE_SC_post_mean
  result_mse_pre[iT0, ] <- MSE_pre_ratio
  result_mse_post[iT0, ] <- MSE_post_ratio

}

end_time <- Sys.time()


