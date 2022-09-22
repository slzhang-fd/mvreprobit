#' @noRd
sample_params_1fac <- function(x_covs, XtX, Y_star, U_all, coeffs,
                               Sigma_e, Sigma_e_inv, Sigma_u, i_ind, cor_step_size){
  K <- ncol(Y_star)
  p <- nrow(coeffs)
  ## update coeffs
  coeffs <- sample_coeffs(x_covs, XtX, Y_star - U_all[i_ind,], Sigma_e_inv)
  ## update Sigma_e
  if(K>1){
    Sigma_e <- sample_sigmae_MH(Sigma_e, Sigma_e_inv,
                                Y_star - x_covs %*% coeffs - U_all[i_ind,], cor_step_size)
  }

  ## update Sigma_u
  Sigma_u <- sample_sigma_direct(U_all)

  list('coeffs' = coeffs,
       'Sigma_e' = Sigma_e,
       'Sigma_u' = Sigma_u)
}

#' @export
mvreprobit_1factor_Gibbs <- function(Y, x_covs, i_ind, max_steps, cor_step_size = 0.01){
  K <- ncol(Y)
  rcd_num <- nrow(Y)
  ind_num <- max(i_ind)
  cov_num <- ncol(x_covs)

  t_len <- aggregate(i_ind, by=list(i_ind), length)$x
  ## initialize random variables
  Y_star <- matrix(0, rcd_num, K)
  U_all <- matrix(rnorm(ind_num*K),ind_num,K)
  ## initialize parameters
  coeffs <- matrix(0, cov_num, K)
  Sigma_e <- Sigma_u <- diag(rep(1,K))

  coeffs_all <- matrix(0, max_steps, cov_num*K)
  Sigma_e_all <- Sigma_u_all <- matrix(0, max_steps, K*K)
  Y_pos_ind <- Y==1
  Y_zero_ind <- Y==0
  XtX <- t(x_covs)%*%x_covs

  init <- numeric(max_steps)
  end <- numeric(max_steps)
  extra <- 6
  width <- 30
  time <- remainining <- 0
  for(iter in 1:max_steps){
    # cat('\r iter: ', iter, 'mu=',coeffs[1,])
    init[iter] <- Sys.time()
    step <- round(iter / max_steps * (width - extra))
    text <- sprintf('|%s%s|% 3s%% | Execution time:%s | Estimated time remaining:%s | rejection rate:%.3f       ',
                    strrep('=', step), strrep(' ', width - step - extra),
                    round(iter / max_steps * 100), my_seconds_to_period(time),
                    my_seconds_to_period(remainining), rejection_rate)
    cat(text)

    Sigma_e_inv <- solve(Sigma_e)
    Sigma_u_inv <- solve(Sigma_u)
    Xbeta <- x_covs %*% coeffs
    ## stochastic E step
    # sample Y_star
    Y_star <- sample_Y_star_3fac(Y_pos_ind, Y_zero_ind, Y_star, Sigma_e,
                                 Xbeta + U_all[i_ind,])
    cat(' sample Y done |')
    # sample U
    U_all <- sample_X_all1_3fac(U_all, Sigma_e_inv, Sigma_u_inv, t_len,
                                rowsum(Y_star - Xbeta, i_ind, reorder = T))
    cat('sample U done |')
    # update parameters
    params <- sample_params_1fac(x_covs, XtX, Y_star, U_all, coeffs,
                                 Sigma_e, Sigma_e_inv, Sigma_u, i_ind, cor_step_size)
    coeffs <- params$coeffs
    Sigma_e <- params$Sigma_e
    Sigma_u <- params$Sigma_u
    # store results
    coeffs_all[iter,] <- coeffs
    Sigma_e_all[iter,] <- Sigma_e
    Sigma_u_all[iter,] <- Sigma_u
    # if(iter %% 100 == 0 && K > 1) cat(' rejection rate:', mean(diff(Sigma_e_all[(iter-99):iter,2])==0))
    ## progress bar
    end[iter] <- Sys.time()
    time <- round(sum(end - init), 0)
    est <- max_steps * (mean(end[end != 0] - init[init != 0])) - time
    remainining <- round(est, 0)
    cat(if (iter == max_steps) '\n' else '\r')
    if(iter %% 100 == 0 && K > 1)
      rejection_rate = mean(diff(Sigma_e_all[(iter-99):iter,2])==0)
  }
  return(list('coeffs_all'=coeffs_all,
              'Sigma_e_all'=Sigma_e_all,
              'Sigma_u_all'=Sigma_u_all,
              'Y_star' = Y_star,
              'U_all' = U_all))
}
#' @export
res_summary_1fac_Gibbs <- function(mvreprobit_res, xcov, burnin){
  K <- sqrt(ncol(mvreprobit_res$Sigma_e_all))
  p <- ncol(mvreprobit_res$coeffs_all) / K
  stem_len <- nrow(mvreprobit_res$coeffs_all)
  cov_names <- names(xcov)
  process_names <- NULL
  if( K == 4 ) process_names <- c('tpprac', 'tpfin', 'fpprac', 'fpfin')
  CORR_u_all <- matrix(0, stem_len - burnin, 16)
  for(iter in (burnin+1):stem_len){
    CORR_u_all[iter-burnin,] <- cov2corr(mvreprobit_res$Sigma_u_all[iter,])
  }
  return(c(summary_res_basic(mvreprobit_res$coeffs_all[-(1:burnin),], 'coeffs', p, K, list(cov_names,process_names)),
           summary_res_basic(mvreprobit_res$Sigma_e_all[-(1:burnin),], 'Sigma_e', K, K, list(process_names,process_names)),
           summary_res_basic(mvreprobit_res$Sigma_u_all[-(1:burnin),], 'Sigma_u', K, K, list(process_names,process_names)),
           summary_res_basic(CORR_u_all, 'CORR_u', K, K, list(process_names,process_names))))
}
#' @export
mvreprobit_1factor_Gibbs_addon <- function(Y, x_covs, i_ind, max_steps, res, cor_step_size = 0.02){
  K <- ncol(Y)
  rcd_num <- nrow(Y)
  ind_num <- max(i_ind)
  cov_num <- ncol(x_covs)
  current_steps <- nrow(res$coeffs_all)

  t_len <- aggregate(i_ind, by=list(i_ind), length)$x
  ## initialize random variables
  Y_star <- res$Y_star
  U_all <- res$U_all
  ## initialize parameters
  coeffs <- matrix(res$coeffs_all[current_steps,], cov_num, K)
  Sigma_e <- matrix(res$Sigma_e_all[current_steps,], K, K)
  Sigma_u <- matrix(res$Sigma_u_all[current_steps,], K, K)

  coeffs_all <- matrix(0, max_steps, cov_num*K)
  Sigma_e_all <- Sigma_u_all <- matrix(0, max_steps, K*K)
  Y_pos_ind <- Y==1
  Y_zero_ind <- Y==0
  XtX <- t(x_covs)%*%x_covs
  for(iter in 1:max_steps){
    cat('\r iter: ', iter, 'mu=',coeffs[1,])
    Sigma_e_inv <- solve(Sigma_e)
    Sigma_u_inv <- solve(Sigma_u)
    Xbeta <- x_covs %*% coeffs
    ## stochastic E step
    # sample Y_star
    Y_star <- sample_Y_star_3fac(Y_pos_ind, Y_zero_ind, Y_star, Sigma_e,
                                 Xbeta + U_all[i_ind,])
    cat(' sample Y done |')
    # sample U
    U_all <- sample_X_all1_3fac(U_all, Sigma_e_inv, Sigma_u_inv, t_len,
                                rowsum(Y_star - Xbeta, i_ind, reorder = T))
    cat('sample U done |')
    # update parameters
    params <- sample_params_1fac(x_covs, XtX, Y_star, U_all, coeffs,
                                 Sigma_e, Sigma_e_inv, Sigma_u, i_ind, cor_step_size)
    coeffs <- params$coeffs
    Sigma_e <- params$Sigma_e
    Sigma_u <- params$Sigma_u
    # store results
    coeffs_all[iter,] <- coeffs
    Sigma_e_all[iter,] <- Sigma_e
    Sigma_u_all[iter,] <- Sigma_u
    if(iter %% 100 == 0 && K > 1) cat(' rejection rate:', mean(diff(Sigma_e_all[1:iter,2])==0))
  }
  return(list('coeffs_all'=coeffs_all,
              'Sigma_e_all'=Sigma_e_all,
              'Sigma_u_all'=Sigma_u_all,
              'Y_star' = Y_star,
              'U_all' = U_all))
}
