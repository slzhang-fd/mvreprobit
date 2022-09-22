#' @export
mvreprobit_3factor_Gibbs <- function(Y, x_covs, j_ind, i_ind, jt_ind, max_steps, cor_step_size = 0.01){
  K <- ncol(Y)
  rcd_num <- nrow(Y)
  ind_num <- max(i_ind)
  grp_num <- max(j_ind)
  jt_num <- max(jt_ind)
  cov_num <- ncol(x_covs)

  t_len <- aggregate(i_ind, by=list(i_ind), length)$x
  it_len <- aggregate(j_ind, by=list(j_ind), length)$x
  i_len <- aggregate(jt_ind, by=list(jt_ind), length)$x
  ## initialize random variables
  Y_star <- matrix(0, rcd_num, K)
  U_all <- matrix(rnorm(ind_num*K),ind_num,K)
  V_all <- matrix(rnorm(grp_num*K),grp_num,K)
  W_all <- matrix(rnorm(jt_num*K),jt_num,K)
  ## initialize parameters
  coeffs <- matrix(0, cov_num, K)
  Sigma_e <- Sigma_v <- Sigma_u <- Sigma_w <- diag(rep(1,K))

  coeffs_all <- matrix(0, max_steps, cov_num*K)
  Sigma_e_all <- Sigma_v_all <- Sigma_u_all <- Sigma_w_all <- matrix(0, max_steps, K*K)
  Y_pos_ind <- Y==1
  Y_zero_ind <- Y==0
  XtX <- t(x_covs)%*%x_covs

  init <- numeric(max_steps)
  end <- numeric(max_steps)
  extra <- 6
  width <- 30
  time <- remainining <- 0
  rejection_rate <- 0
  for(iter in 1:max_steps){
    init[iter] <- Sys.time()
    step <- round(iter / max_steps * (width - extra))
    text <- sprintf('|%s%s|% 3s%% | Execution time:%s | Estimated time remaining:%s | rejection rate:%.3f       ',
                    strrep('=', step), strrep(' ', width - step - extra),
                    round(iter / max_steps * 100), my_seconds_to_period(time),
                    my_seconds_to_period(remainining), rejection_rate)
    cat(text)
    Sigma_e_inv <- solve(Sigma_e)
    Sigma_u_inv <- solve(Sigma_u)
    Sigma_v_inv <- solve(Sigma_v)
    Sigma_w_inv <- solve(Sigma_w)
    Xbeta <- x_covs %*% coeffs
    ## stochastic E step
    # sample Y_star
    Y_star <- sample_Y_star_3fac(Y_pos_ind, Y_zero_ind, Y_star, Sigma_e,
                                 Xbeta + V_all[j_ind,] + U_all[i_ind,] + W_all[jt_ind,])
    # cat(' sample Y done |')
    # sample U
    U_all <- sample_X_all1_3fac(U_all, Sigma_e_inv, Sigma_u_inv, t_len,
                                rowsum(Y_star - Xbeta - V_all[j_ind,] - W_all[jt_ind,], i_ind, reorder = T))
    # cat('sample U done |')
    # sample V
    V_all <- sample_X_all1_3fac(V_all, Sigma_e_inv, Sigma_v_inv, it_len,
                                rowsum(Y_star - Xbeta - U_all[i_ind,] - W_all[jt_ind,], j_ind, reorder = T))
    # cat('sample V done |')
    # sample W
    W_all <- sample_X_all1_3fac(W_all, Sigma_e_inv, Sigma_w_inv, i_len,
                                rowsum(Y_star - Xbeta - U_all[i_ind,] - V_all[j_ind,], jt_ind, reorder = T))
    # cat('sample W done')
    # update parameters
    params <- sample_params_3fac(x_covs, XtX, Y_star, U_all, V_all, W_all, coeffs,
                                 Sigma_e, Sigma_e_inv, Sigma_u, Sigma_v, Sigma_w, j_ind, i_ind, jt_ind, cor_step_size)
    coeffs <- params$coeffs
    Sigma_e <- params$Sigma_e
    Sigma_u <- params$Sigma_u
    Sigma_v <- params$Sigma_v
    Sigma_w <- params$Sigma_w
    # store results
    coeffs_all[iter,] <- coeffs
    Sigma_e_all[iter,] <- Sigma_e
    Sigma_u_all[iter,] <- Sigma_u
    Sigma_v_all[iter,] <- Sigma_v
    Sigma_w_all[iter,] <- Sigma_w

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
              'Sigma_v_all'=Sigma_v_all,
              'Sigma_w_all'=Sigma_w_all,
              'Y_star'=Y_star,
              'U_all'=U_all,
              'V_all'=V_all,
              'W_all'=W_all))
}

#' @export
mvreprobit_3factor_Gibbs_addon <- function(Y, x_covs, j_ind, i_ind, jt_ind, max_steps, res, cor_step_size = 0.02){
  K <- ncol(Y)
  rcd_num <- nrow(Y)
  ind_num <- max(i_ind)
  grp_num <- max(j_ind)
  jt_num <- max(jt_ind)
  cov_num <- ncol(x_covs)
  current_steps <- nrow(res$coeffs_all)

  t_len <- aggregate(i_ind, by=list(i_ind), length)$x
  it_len <- aggregate(j_ind, by=list(j_ind), length)$x
  i_len <- aggregate(jt_ind, by=list(jt_ind), length)$x
  ## initialize random variables
  Y_star <- res$Y_star
  U_all <- res$U_all
  V_all <- res$V_all
  W_all <- res$W_all
  ## initialize parameters
  coeffs <- matrix(res$coeffs_all[current_steps,], cov_num, K)
  Sigma_e <- matrix(res$Sigma_e_all[current_steps,], K, K)
  Sigma_v <- matrix(res$Sigma_v_all[current_steps,], K, K)
  Sigma_u <- matrix(res$Sigma_u_all[current_steps,], K, K)
  Sigma_w <- matrix(res$Sigma_w_all[current_steps,], K, K)

  coeffs_all <- matrix(0, max_steps, cov_num*K)
  Sigma_e_all <- Sigma_v_all <- Sigma_u_all <- Sigma_w_all <- matrix(0, max_steps, K*K)
  Y_pos_ind <- Y==1
  Y_zero_ind <- Y==0
  for(iter in 1:max_steps){
    cat('\r iter: ', iter, 'mu=',coeffs[1,])
    Sigma_e_inv <- solve(Sigma_e)
    Sigma_u_inv <- solve(Sigma_u)
    Sigma_v_inv <- solve(Sigma_v)
    Sigma_w_inv <- solve(Sigma_w)
    Xbeta <- x_covs %*% coeffs
    ## stochastic E step
    # sample Y_star
    Y_star <- sample_Y_star_3fac(Y_pos_ind, Y_zero_ind, Y_star, Sigma_e,
                                 Xbeta + V_all[j_ind,] + U_all[i_ind,] + W_all[jt_ind,])
    cat(' sample Y done |')
    # sample U
    U_all <- sample_X_all1_3fac(U_all, Sigma_e_inv, Sigma_u_inv, t_len,
                                rowsum(Y_star - Xbeta - V_all[j_ind,] - W_all[jt_ind,], i_ind, reorder = T))
    cat('sample U done |')
    # sample V
    V_all <- sample_X_all1_3fac(V_all, Sigma_e_inv, Sigma_v_inv, it_len,
                                rowsum(Y_star - Xbeta - U_all[i_ind,] - W_all[jt_ind,], j_ind, reorder = T))
    cat('sample V done |')
    # sample W
    W_all <- sample_X_all1_3fac(W_all, Sigma_e_inv, Sigma_w_inv, i_len,
                                rowsum(Y_star - Xbeta - U_all[i_ind,] - V_all[j_ind,], jt_ind, reorder = T))
    cat('sample W done')
    # update parameters
    params <- sample_params_3fac(x_covs, Y_star, U_all, V_all, W_all, coeffs,
                                 Sigma_e, Sigma_e_inv, Sigma_u, Sigma_v, Sigma_w, j_ind, i_ind, jt_ind, cor_step_size)
    coeffs <- params$coeffs
    Sigma_e <- params$Sigma_e
    Sigma_u <- params$Sigma_u
    Sigma_v <- params$Sigma_v
    Sigma_w <- params$Sigma_w
    # store results
    coeffs_all[iter,] <- coeffs
    Sigma_e_all[iter,] <- Sigma_e
    Sigma_u_all[iter,] <- Sigma_u
    Sigma_v_all[iter,] <- Sigma_v
    Sigma_w_all[iter,] <- Sigma_w
    if(iter %% 100 == 0 && K > 1) cat(' rejection rate:', mean(diff(Sigma_e_all[1:iter,2])==0))
  }
  return(list('coeffs_all'=coeffs_all,
              'Sigma_e_all'=Sigma_e_all,
              'Sigma_u_all'=Sigma_u_all,
              'Sigma_v_all'=Sigma_v_all,
              'Sigma_w_all'=Sigma_w_all,
              'Y_star'=Y_star,
              'U_all'=U_all,
              'V_all'=V_all,
              'W_all'=W_all))
}
#' @export
res_summary_3fac_Gibbs <- function(mvreprobit_res, xcov, burnin){
  K <- sqrt(ncol(mvreprobit_res$Sigma_e_all))
  p <- ncol(mvreprobit_res$coeffs_all) / K
  stem_len <- nrow(mvreprobit_res$coeffs_all)
  cov_names <- names(xcov)
  process_names <- NULL
  process_names <- c('tpprac', 'tpfin', 'fpprac', 'fpfin')[1:K]
  CORR_u_all <- matrix(0, stem_len - burnin, K*K)
  for(iter in (burnin+1):stem_len){
    CORR_u_all[iter-burnin,] <- cov2corr(mvreprobit_res$Sigma_u_all[iter,])
  }
  return(c(summary_res_basic(mvreprobit_res$coeffs_all[-(1:burnin),], 'coeffs', p, K, list(cov_names,process_names)),
           summary_res_basic(mvreprobit_res$Sigma_e_all[-(1:burnin),], 'Sigma_e', K, K, list(process_names,process_names)),
           summary_res_basic(mvreprobit_res$Sigma_u_all[-(1:burnin),], 'Sigma_u', K, K, list(process_names,process_names)),
           summary_res_basic(mvreprobit_res$Sigma_v_all[-(1:burnin),], 'Sigma_v', K, K, list(process_names,process_names)),
           summary_res_basic(mvreprobit_res$Sigma_w_all[-(1:burnin),], 'Sigma_w', K, K, list(process_names,process_names)),
           summary_res_basic(CORR_u_all, 'CORR_u', K, K, list(process_names,process_names))))
}
