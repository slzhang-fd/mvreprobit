#' @importFrom truncnorm rtruncnorm
#' @importFrom mvtnorm rmvnorm
#' @noRd
sample_Y_star_3re <- function(Y_pos_ind, Y_zero_ind, Y_star, Sigma_e, mu_v_u_w){
  K <- ncol(Y_star)
  if(K==1){
    Y_star[Y_pos_ind[,1],1] <- rtruncnorm(1, a=0, mean = mu_v_u_w[Y_pos_ind[,1]], sd = 1)
    Y_star[Y_zero_ind[,1],1] <- rtruncnorm(1, b=0, mean = mu_v_u_w[Y_zero_ind[,1]], sd = 1)
  } else{
    for(k in 1:K){
      sigma11 <- Sigma_e[k,k]
      sigma12 <- Sigma_e[k,-k]
      sigma22_inv <- solve(Sigma_e[-k,-k])
      sigma_Y_cond <- sigma11 - sigma12 %*% sigma22_inv %*% sigma12
      mu_Y_cond <- mu_v_u_w[,k] + (Y_star[,-k] - mu_v_u_w[,-k]) %*% sigma22_inv %*% sigma12

      Y_star[Y_pos_ind[,k],k] <- rtruncnorm(1, a=0, mean = mu_Y_cond[Y_pos_ind[,k]], sd = sqrt(sigma_Y_cond))
      Y_star[Y_zero_ind[,k],k] <- rtruncnorm(1, b=0, mean = mu_Y_cond[Y_zero_ind[,k]], sd = sqrt(sigma_Y_cond))
    }
  }

  Y_star
}
#' @noRd
sample_X_all1_3re <- function(X_all, Sigma_e_inv, Sigma_x_inv, x_len, temp){
  # temp <- rowsum(Y_star - Xbeta - U_all[i_ind,] - V_all[j_ind,], jt_ind, reorder = T)
  for(pp in unique(x_len)){
    sigma_X_cond <- solve(pp * Sigma_e_inv + Sigma_x_inv)
    x_loc <- which(x_len==pp)
    mu_X_cond <- temp[x_loc,] %*% Sigma_e_inv %*% sigma_X_cond
    X_all[x_loc,] <- mu_X_cond + rmvnorm(length(x_loc), sigma = sigma_X_cond)
  }
  X_all
}
#' @noRd
Optim_step_3re <- function(x_covs, Y_star, U_all, V_all, W_all, coeffs,
                            B_e, Sigma_e_inv, B_u, B_v, B_w, j_ind, i_ind, jt_ind){
  K <- ncol(Y_star)
  p <- nrow(coeffs)
  ## update coeffs
  params <- optim(par = c(coeffs), fn = neg_loglik_Y_star_coeff, gr = neg_loglik_Y_star_coeff_deriv,
                  x_covs=x_covs, temp0=Y_star - V_all[j_ind,] - U_all[i_ind,] - W_all[jt_ind,],
                  Sigma_e_inv=Sigma_e_inv,
                  method = "BFGS")$par
  coeffs <- matrix(params,p,K)
  ## update Sigma_e
  if(K>1){
    params <- optim(par = c(B_e[lower.tri(B_e)]), fn = neg_loglik_Y_star_sigma,
                    temp=Y_star - x_covs %*% coeffs - V_all[j_ind,] - U_all[i_ind,] - W_all[jt_ind,],
                    method = "BFGS")$par
    B_e[lower.tri(B_e)] <- params
  }
  # B_e <- diag(1/sqrt(rowSums(B_e^2))) %*% B_e

  ## update Sigma_u
  params <- optim(par = B_u[lower.tri(B_u,diag = T)], fn = neg_loglik_mvnorm,
                  x = U_all, method = "BFGS")$par
  B_u[lower.tri(B_u, diag = T)] <- params

  ## update Sigma_v
  params <- optim(par = B_v[lower.tri(B_v,diag = T)], fn = neg_loglik_mvnorm,
                  x = V_all, method = "BFGS")$par
  B_v[lower.tri(B_v, diag = T)] <- params

  ## update Sigma_w
  params <- optim(par = B_w[lower.tri(B_w,diag = T)], fn = neg_loglik_mvnorm,
                  x = W_all, method = "BFGS")$par
  B_w[lower.tri(B_w, diag = T)] <- params

  list('coeffs' = coeffs,
       'B_e' = B_e,
       'B_u' = B_u,
       'B_v' = B_v,
       'B_w' = B_w)
}

#' @importFrom MCMCpack riwish
#' @noRd
sample_sigma_direct <- function(X){
  ## inverse gamma
  # arma::mat sample_Sigma_new1(const arma::mat &e1234){
  #   arma::mat S = arma::eye(4,4) + e1234.t() * e1234;
  #   return riwish(4 + e1234.n_rows, S);
  # }
  K <- ncol(X)
  N <- nrow(X)
  S <- diag(K) + t(X) %*% X
  return( riwish(K + N, S) )
}
#' @importFrom matrixcalc is.positive.definite
#' @noRd
sample_sigmae_MH <- function(Sigma, Sigma_inv, X, step_size = 0.01){
  K <- ncol(X)
  N <- nrow(X)
  for(kk in 1:(K*(K-1)/2)){
    perturb <- step_size * rnorm(1)
    Sigma_new <- Sigma
    Sigma_new[lower.tri(Sigma_new)][kk] <- Sigma_new[lower.tri(Sigma_new)][kk] + perturb
    Sigma_new[upper.tri(Sigma_new)] <- t(Sigma_new)[upper.tri(Sigma_new)]
    if(is.positive.definite(Sigma_new)){
      Sigma_new_inv <- solve(Sigma_new)
      eet <- t(X) %*% X
      log_odds <- #sum(dnorm(Sigma_new[lower.tri(Sigma_new)], sd = 10, log = T)) -
        #sum(dnorm(Sigma[lower.tri(Sigma)], sd = 10, log = T)) +
        -0.5 * N * (as.numeric(determinant(Sigma_new)$modulus) - as.numeric(determinant(Sigma)$modulus)) -
        0.5 * (sum(diag(eet %*% Sigma_new_inv)) - sum(diag(eet %*% Sigma_inv)))
      if(runif(1) < exp(log_odds)){
        Sigma <- Sigma_new
      }
    }
  }
  return(Sigma)
}
#' @noRd
sample_sigmae_MH0 <- function(Sigma, Sigma_inv, X, step_size = 0.007){
  K <- ncol(X)
  N <- nrow(X)
  perturb <- step_size * rnorm((K*(K-1)/2))
  Sigma_new <- Sigma
  Sigma_new[lower.tri(Sigma_new)] <- Sigma_new[lower.tri(Sigma_new)] + perturb
  Sigma_new[upper.tri(Sigma_new)] <- t(Sigma_new)[upper.tri(Sigma_new)]
  if(is.positive.definite(Sigma_new)){
    Sigma_new_inv <- solve(Sigma_new)
    eet <- t(X) %*% X
    log_odds <- #sum(dnorm(Sigma_new[lower.tri(Sigma_new)], sd = 10, log = T)) -
      #sum(dnorm(Sigma[lower.tri(Sigma)], sd = 10, log = T)) +
      -0.5 * N * (as.numeric(determinant(Sigma_new)$modulus) - as.numeric(determinant(Sigma)$modulus)) -
      0.5 * (sum(diag(eet %*% Sigma_new_inv)) - sum(diag(eet %*% Sigma_inv)))
    if(runif(1) < exp(log_odds)){
      Sigma <- Sigma_new
    }
  }
  return(Sigma)
}
#' @noRd
sample_coeffs <- function(X, XtX, Z, Sigma_inv){
  K <- ncol(Z)
  N <- nrow(Z)
  sig2_beta0 <- 100
  V_beta <- kronecker(Sigma_inv, XtX)
  diag(V_beta) <- diag(V_beta) + 1.0/sig2_beta0
  V_beta <- solve(V_beta)
  mu_beta <- V_beta %*% kronecker(Sigma_inv, t(X)) %*% as.vector(Z)
  params <- mu_beta + t(rmvnorm(1, sigma = V_beta))
  return(matrix(params, ncol = K))
}
#' @noRd
sample_params_3re <- function(x_covs, XtX, Y_star, U_all, V_all, W_all, coeffs,
                               Sigma_e, Sigma_e_inv, Sigma_u, Sigma_v, Sigma_w, j_ind, i_ind, jt_ind, sample_params_3re){
  K <- ncol(Y_star)
  p <- nrow(coeffs)
  ## sample coeffs
  coeffs <- sample_coeffs(x_covs, XtX, Y_star - V_all[j_ind,] - U_all[i_ind,] - W_all[jt_ind,], Sigma_e_inv)
  ## sample Sigma_e
  if(K > 1){
    Sigma_e <- sample_sigmae_MH(Sigma_e, Sigma_e_inv,
                                Y_star - x_covs %*% coeffs - V_all[j_ind,] - U_all[i_ind,] - W_all[jt_ind,], sample_params_3re)
  }

  ## sample Sigma_u
  Sigma_u <- sample_sigma_direct(U_all)

  ## sample Sigma_v
  Sigma_v <- sample_sigma_direct(V_all)

  ## sample Sigma_w
  Sigma_w <- sample_sigma_direct(W_all)

  list('coeffs' = coeffs,
       'Sigma_e' = Sigma_e,
       'Sigma_u' = Sigma_u,
       'Sigma_v' = Sigma_v,
       'Sigma_w' = Sigma_w)
}
#' @noRd
my_seconds_to_period = function(x) {
  days = round(x %/% (60 * 60 * 24))
  hours = round((x - days*60*60*24) %/% (60 * 60))
  minutes = round((x - days*60*60*24 - hours*60*60) %/% 60)
  seconds = round(x - days*60*60*24 - hours*60*60 - minutes*60)
  days_str = ifelse(days == 0, "", paste0(days, "d "))
  hours_str = ifelse((hours == 0 & days == 0), "", paste0(hours, "hour "))
  minutes_str = ifelse((minutes == 0 & days == 0 & hours == 0), "", paste0(minutes, "min "))
  seconds_str = paste0(seconds, "s")
  final_str = paste0(days_str, hours_str, minutes_str, seconds_str)
  return(final_str)
}
#' @noRd
cov2corr <- function(x){
  K <- sqrt(length(c(x)))
  x <- matrix(c(x), K, K)
  scale <- diag(1 / sqrt(diag(x)), nrow = K, ncol = K)
  return(c(scale %*% x %*% scale))
}
#' @noRd
summary_res_basic <- function(x, name, rownum, colnum, dimnames){
  res <- list(matrix(colMeans(matrix(x, ncol = rownum * colnum)), rownum, colnum, dimnames=dimnames),
              matrix(apply(matrix(x, ncol = rownum * colnum), 2, sd), rownum, colnum, dimnames=dimnames),
              matrix(apply(matrix(x, ncol = rownum * colnum), 2, quantile, 0.025), rownum, colnum, dimnames=dimnames),
              matrix(apply(matrix(x, ncol = rownum * colnum), 2, quantile, 0.975), rownum, colnum, dimnames=dimnames))
  names(res) <- paste0(name, c('', '_sd', '_low2.5', '_high97.5'))
  return(res)
}
#' @export
res_compare_plot <- function(refit_res, true_val, burn_in) {
  refit_res <- refit_res[(burn_in + 1):nrow(refit_res), ]
  yl <- apply(refit_res, 2, quantile, 0.025)
  yu <- apply(refit_res, 2, quantile, 0.975)
  ## to remove warnings of arrows commmand
  yu <- pmax(yl+mean(yl)*0.001, yu)
  x <- seq_len(ncol(refit_res))
  plot(x, true_val, ylim = range(c(yl, yu)), pch = 19, cex = .3, xaxt = "n")
  axis(1, at = x)
  arrows(x, yl, x, yu, angle = 90, code = 3, length = .04)
}
