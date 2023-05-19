
cavi_msfa <- function(X_s,
                    K,
                    J_s,
                    verbose = 1,
                    min_iter = 5,
                    tol = 5e-3,
                    max_iter = 100,
                    scale = TRUE) {
  
  S <- length(X_s)
  
  p_dims <- sapply(X_s, function(x) dim(x)[2])
  if(length(unique(p_dims)) > 1) {
    warning("Each study must have the same number of variables!")
  } else{
    P <- unique(p_dims)
  }
  
  N_s <- sapply(X_s, function(x) dim(x)[1])
  
  if(scale){
    X_scaled <- lapply(X_s, function(x) scale(x, scale = TRUE, center = TRUE))
  } else {
    X_scaled <-lapply(X_s, function(x) scale(x, scale = FALSE, center = TRUE))
  }
  
  # Update Priors ###########################################
  update_omega_lambda <- function(S, P, J, v_s, mean_lambda, var_lambda, dist_delta){
    lapply(1:S, function(s){
      tau <- cumprod(dist_delta[[s]][3,])
      lapply(1:P, function(p){sapply(1:J[s], function(j){
        alpha <- (3 + 1)/2
        beta <- (3 + tau[j]*(mean_lambda[[s]][[p]][j]^2 + var_lambda[[s]][[p]][j,j]))/2
        
        mean <- alpha/beta
        
        return(c(alpha, beta, mean))
      })
      })
    })
  }
  
  update_omega_phi <- function(P, K, v_s, mean_Phi, var_Phi, dist_delta){
    tau <- cumprod(dist_delta[3,])
    lapply(1:P, function(p){sapply(1:K, function(k) {
      alpha <- (3 + 1)/2
      beta <- (3 + tau[k]*(mean_Phi[[p]][k]^2 + var_Phi[[p]][k,k]))/2
      mean <- alpha/beta
      
      return(c(alpha, beta, mean))
    } )})
  }
  
  update_shrinkage <- function(P, dim, delta, prior_delta, omega, mean_loading, var_loading){
    delta_update <- lapply(1:dim, function(l){
      delta_l <- delta[3,]
      delta_l[l] <- 1
      
      tau <- cumprod(delta_l)
      
      alpha <- prior_delta[1,l] + 1/2 * P * (dim - l + 1)
      beta <- prior_delta[2,l] + 1/2 * sum(unlist(lapply(1:P, function(p){lapply(l:dim, function(k){
        omega[[p]][3,k] * tau[k] * (mean_loading[[p]][k]^2 + var_loading[[p]][k,k])
      })})))
      
      mean_delta <- alpha/beta
      
      
      delta[1,l] <- alpha
      delta[2,l] <- beta
      delta[3,l] <- alpha/beta
      
      
      return(c(alpha, beta, mean_delta))
    })
    
    result <- matrix(unlist(delta_update), nrow=3, ncol=dim, byrow=FALSE)
    
    return(result)
  }
  
  # Update latent factors ###################################
  update_factors <- function(S, sizes, X_samp,  mean_f_samp, mean_l_samp, P, K, J_s, mean_Psi, mean_Phi, var_Phi, mean_Lambda, var_Lambda, n_iter){
    phi_mat <- matrix(unlist(mean_Phi), nrow=P, ncol=K, byrow=TRUE)
    
    dists <- lapply(1:S, function(s){
      Lambda_s <-matrix(unlist(mean_Lambda[[s]]), nrow=P, ncol=J_s[s], byrow=TRUE)
      
      var_ls <- solve(
        diag(J_s[s]) + Reduce('+', lapply(1:P, function(p) mean_Psi[[s]][p] * (tcrossprod(mean_Lambda[[s]][[p]]) + var_Lambda[[s]][[p]])))
      )
      
      var_fs <- solve(
        diag(K) + Reduce('+', lapply(1:P, function(p) mean_Psi[[s]][p] * (tcrossprod(mean_Phi[[p]]) + var_Phi[[p]])))
      )
      
      mean_fs <- matrix(unlist(mean_f_samp[[s]]), nrow=K, ncol=sizes[s], byrow=FALSE)
      mean_ls <- matrix(unlist(mean_l_samp[[s]]), nrow=J_s[s], ncol=sizes[s], byrow=FALSE)
      
      
      mean_ls_update <- var_ls %*% t(Lambda_s) %*% diag(mean_Psi[[s]]) %*% (t(X_samp[[s]]) - phi_mat %*% mean_fs)
      mean_fs_update <- var_fs %*% t(phi_mat) %*% diag(mean_Psi[[s]]) %*% (t(X_samp[[s]]) - Lambda_s %*% mean_ls_update)
      
      mean_lis_update <- lapply(1:sizes[s], function(i){matrix(mean_ls_update[,i], ncol=1)})
      mean_fis_update <- lapply(1:sizes[s], function(i){matrix(mean_fs_update[,i], ncol=1)})
      
      return(list("mean_l"=mean_lis_update, "var_l"=var_ls, "mean_f"=mean_fis_update, "var_f"=var_fs))
    })
    return(dists)
  }
  
  # Update Factor Loadings  ################################
  update_Phi <- function(S, P, X, N_s, K, J_s, dist_omega_phi, dist_delta, mean_Psi, mean_f, var_f, mean_l, mean_lambda){
    
    tau <- cumprod(dist_delta[3,])
    
    # Caluclate variance for each lambda_ps
    var_phi <- lapply(1:P, function(p){
      Sum_matrix <- Reduce('+',lapply(1:S, function(s) {
        mean_Psi[[s]][p] * Reduce('+', lapply(1:N_s[s], function(n) var_f[[s]] + tcrossprod(mean_f[[s]][[n]])))
      }))
      
      
      Dp_matrix <- diag(dist_omega_phi[[p]][3,] * tau)
      
      
      var_phi_p <- solve(Dp_matrix + Sum_matrix)
    })
    
    mean_phi <- lapply(1:P, function(p){
      diff_term <- Reduce('+', 
                          lapply(1:S, function(s){
                            x_sp <- matrix(X[[s]][,p], ncol = 1, nrow=N_s[s])
                            f_s <- matrix(unlist(mean_f[[s]]), nrow = N_s[s], ncol=K, byrow=TRUE)
                            l_s <- matrix(unlist(mean_l[[s]]), nrow = N_s[s], ncol=J_s[s], byrow=TRUE)
                            mean_Psi[[s]][p] * t(f_s) %*% (x_sp - l_s %*% mean_lambda[[s]][[p]])
                          })
      )
      matrix(var_phi[[p]] %*% diff_term, nrow=K, ncol=1)
    })
    
    return(list(mean_phi, var_phi))
  }
  
  update_Lambda <- function(S, P, X, N_s, K, J_s, dist_omega_lambda, dist_delta, mean_Psi, mean_l, var_l, mean_f, mean_Phi){
    
    # Caluclate variance for each lambda_ps
    var_lambda <- lapply(1:S, function(s) lapply(1:P, function(p){
      solve(diag(dist_omega_lambda[[s]][[p]][3,]*cumprod(dist_delta[[s]][3,])) +
              mean_Psi[[s]][p] * Reduce('+', lapply(1:N_s[s], function(n) var_l[[s]] + tcrossprod(mean_l[[s]][[n]])))
      )
    }))
    
    mean_lambda <- lapply(1:S, function(s) lapply(1:P, function(p){
      x_sp <- matrix(X[[s]][,p], ncol = 1, nrow=N_s[s])
      f_s <- matrix(unlist(mean_f[[s]]), nrow = N_s[s], ncol=K, byrow=TRUE)
      l_s <- matrix(unlist(mean_l[[s]]), nrow = N_s[s], ncol=J_s[s], byrow=TRUE)
      
      diff_term <- mean_Psi[[s]][p] * t(l_s) %*% (x_sp - f_s %*% mean_Phi[[p]])
      
      matrix(var_lambda[[s]][[p]] %*% diff_term, ncol=1, nrow=J_s[s]) 
    }))
    
    return(list(mean_lambda, var_lambda))
  }
  
  # Update Errors  #########################################
  update_Psi <- function(x, n, p, Psi_prior, mean_Phi, var_Phi, mean_f, var_f, mean_lambda, var_lambda, mean_l, var_l){
    a_psi_update <- Psi_prior[1] + n/2
    
    b_psi_update <- Psi_prior[2] + 1/2 * sapply(1:p, function(p){
      sq_sum <- sum(
        sapply(1:n, function(i){
          (x[i,p] - t(mean_Phi[[p]]) %*% mean_f[[i]] - t(mean_lambda[[p]]) %*% mean_l[[i]])^2
        })
      )
      
      f_quad_term <- Reduce("+", lapply(1:n, function(i) tcrossprod(mean_f[[i]])))
      l_quad_term <- Reduce("+", lapply(1:n, function(i) tcrossprod(mean_l[[i]])))
      
      f_phi_quad_term <- t(mean_Phi[[p]]) %*% (n * var_f) %*% mean_Phi[[p]]
      l_lambda_quad_term <- t(mean_lambda[[p]]) %*% ( n * var_l)  %*% mean_lambda[[p]]
      
      tr_term_f <- sum(
        diag(
          (f_quad_term + n * var_f) %*% var_Phi[[p]]
        )
      )
      tr_term_l <- sum(
        diag(
          (l_quad_term + n * var_l) %*% var_lambda[[p]]
        )
      )
      
      sq_sum + f_phi_quad_term + l_lambda_quad_term + tr_term_f + tr_term_l
    })
    return(list(a_psi_update,b_psi_update,a_psi_update/b_psi_update))
  }
  
  # Helper Functions  ######################################
  convergence_criteria <- function(delta, n_iter, min_iter){
    num_params <- P*(2*K + 2*sum(J_s) + 1 + K^2 + sum(J_s^2)) 
    if(res <- n_iter < min_iter){1}else{res <- sqrt(delta[n_iter])/(num_params)}
    return(res)
  }
  
  start_vi_sparse <- function(X_s, P, K, J_s, alpha=1e-5, beta=1e-5){
    S = length(X_s)
    N_s = sapply(X_s, nrow)
    spca_phi <- sparsepca::rspca(Reduce(rbind, X_s), k = K, alpha=alpha, beta=beta)
    phi_est <- matrix(spca_phi$loadings, nrow=P, ncol=K)
    
    spca_scores <- spca_phi$scores
    scaling_factors <- apply(spca_scores, 2, sd)
    
    for (k in 1:K){
      spca_scores[,k] = spca_scores[,k] / scaling_factors[k]
      phi_est[,k] = phi_est[,k] * scaling_factors[k]
    }
    
    
    spca_residuals <- Reduce(rbind, X_s) - spca_scores %*% t(phi_est)
    
    f_initial <-  lapply(1:S, function(s, n_s){
      start <- 1 + ifelse(s>1, sum(n_s[1:(s-1)]) ,0)
      end <- sum(n_s[1:s])
      spca_scores[start:end,]
    },
    n_s=N_s
    )
    
    X_s_residuals <- lapply(1:S, function(s, n_s){
      start <- 1 + ifelse(s>1, sum(n_s[1:(s-1)]) ,0)
      end <- sum(n_s[1:s])
      spca_residuals[start:end,]
    },
    n_s=N_s
    )
    
    
    lambda_est <- list()
    l_initial <- list()
    
    for (s in 1:S){
      spca_lambda <- sparsepca::rspca(X_s_residuals[[s]], k=J_s[s], alpha=alpha, beta=beta)
      lambda_est[[s]] <- matrix(spca_lambda$loadings,  nrow=P, ncol=J_s[s])
      l_initial[[s]] <- spca_lambda$scores
      
      scaling_factors <- apply(l_initial[[s]], 2, sd)
      for(j in 1:J_s[s]){
        lambda_est[[s]][,j] <- lambda_est[[s]][,j] * scaling_factors[j]
        l_initial[[s]][,j] <- l_initial[[s]][,j] / scaling_factors[j]
      }
    }
    
    psi_est <-  lapply(1:S, function(s){
      abs(diag(cov(X_s[[s]]) - tcrossprod(phi_est) - tcrossprod(lambda_est[[s]]))) + 1e-5
    })
    
    return(list("Phi"=phi_est, "Lambda_s"=lambda_est, "Psi_s"=psi_est, "f" = f_initial, "l" = l_initial))
  }
  
  # Initialize Parameters ##################################
  prior_omega <- c(3/2, 3/2)
  prior_psi <- c(1, 0.3)
  
  dist_omega_lambda_prior <- lapply(1:S, function(s){lapply(1:P, function(p) matrix(c(rep(3/2, J_s[s]), rep(3/2, J_s[s]), rep(1, J_s[s])), nrow=3, ncol=J_s[s], byrow=TRUE))})
  dist_omega_lambda <- dist_omega_lambda_prior
  
  dist_omega_phi_prior <- lapply(1:P, function(p) matrix(c(rep(3/2, K), rep(3/2, K), rep(1, K)), nrow=3, ncol=K, byrow=TRUE))
  dist_omega_phi <- dist_omega_phi_prior
  
  a1 <- 2.1
  a2 <- 3.1
  
  dist_delta_phi_prior <- matrix(
    c(c(a1, rep(a2, K-1)), rep(1, K), c(a1, rep(a2, K-1))), nrow=3, ncol=K, byrow=TRUE
  )
  dist_delta_phi <- dist_delta_phi_prior
  
  dist_delta_lambda_prior <- lapply(1:S, function(s, j){
    matrix(
      c(c(a1, rep(a2, j[s]-1)), rep(1, j[s]), c(a1, rep(a2, j[s]-1))), nrow=3, ncol=j[s], byrow=TRUE
    )
  }, j = J_s)
  dist_delta_lambda <- dist_delta_lambda_prior
  
  initial_params <- start_vi_sparse(X_scaled, P, K, J_s)
  
  mean_Phi <- lapply(1:P, function(x) matrix(initial_params$Phi[x,], ncol=1, nrow=K))
  var_Phi_prior <- lapply(1:P, function(p) diag(dist_omega_phi[[p]][3,]^(-1)*cumprod(dist_delta_phi[3,])^(-1)))
  var_Phi <- var_Phi_prior
  
  mean_Lambda <- lapply(1:S, function(s) lapply(1:P, function(p) matrix(initial_params$Lambda_s[[s]][p,], nrow=J_s[s], ncol=1)))
  var_Lambda_prior <- lapply(1:S, function(s){
    lapply(1:P, function(p){
      diag(cumprod(dist_delta_lambda[[s]][3,])^(-1) * dist_omega_lambda[[s]][[p]][3,]^(-1))
    })
  })
  var_Lambda <- var_Lambda_prior
  
  mean_Psi <- lapply(1:S, function(s) initial_params$Psi_s[[s]]^(-1))
  dist_Psi <- lapply(1:S, function(s){
    alpha <- prior_psi[1] + N_s[s]/2
    beta <- alpha/mean_Psi[[s]]
    mean <- alpha/beta
    return(list(alpha,beta,mean))
  })
  
  mean_f <- lapply(1:S, function(s){
    lapply(1:N_s[[s]], function(i){
      matrix(initial_params$f[[s]][i,], ncol=1)
    })
  })

  mean_l <- lapply(1:S, function(s){
    lapply(1:N_s[[s]], function(i){
      matrix(initial_params$l[[s]][i,], ncol=1)
    })
  })

  
  var_l <- lapply(1:S, function(s){
    solve(
      diag(J_s[s]) + Reduce('+', lapply(1:P, function(p) mean_Psi[[s]][p] * (tcrossprod(mean_Lambda[[s]][[p]]) + var_Lambda[[s]][[p]]) ))
    )
  })
  
  var_f <- lapply(1:S, function(s){
    solve(
      diag(K) + Reduce('+', lapply(1:P, function(p) mean_Psi[[s]][p] * (tcrossprod(mean_Phi[[p]]) + var_Phi[[p]])))
    )
  })
  
  phi_update_initial <- update_Phi(
    S = S,
    P = P,
    X = X_scaled,
    N_s = N_s,
    K = K,
    J_s = J_s,
    dist_omega_phi = dist_omega_phi,
    dist_delta = dist_delta_phi,
    mean_Psi = mean_Psi,
    mean_f = mean_f,
    var_f = var_f,
    mean_l = mean_l,
    mean_lambda = mean_Lambda
  )
  
  mean_Phi = phi_update_initial[[1]]
  var_Phi = phi_update_initial[[2]]
  
  n_iter = 0
  delta_trace <- rep(NA, max_iter)
  
  while(TRUE){
    n_iter <- n_iter + 1
    if(n_iter > max_iter){break}
    
    Lambda_params_update <- update_Lambda(
      S = S,
      P = P,
      X = X_scaled,
      N_s = N_s,
      K = K,
      J_s = J_s,
      dist_omega_lambda = dist_omega_lambda,
      dist_delta = dist_delta_lambda,
      mean_Psi = mean_Psi,
      mean_l = mean_l,
      var_l = var_l,
      mean_f = mean_f,
      mean_Phi = mean_Phi
    )
    
    Phi_params_update <- update_Phi(
      S = S,
      P = P,
      X = X_scaled,
      N_s = N_s,
      K = K,
      J_s = J_s,
      dist_omega_phi = dist_omega_phi,
      dist_delta = dist_delta_phi,
      mean_Psi = mean_Psi,
      mean_f = mean_f,
      var_f = var_f,
      mean_l = mean_l,
      mean_lambda = Lambda_params_update[[1]]
    )
    
    Psi_params_update <- lapply(1:S, function(s) {
      update_Psi(
        x = X_scaled[[s]],
        n = N_s[s],
        p = P,
        Psi_prior = prior_psi,
        mean_Phi = Phi_params_update[[1]],
        var_Phi = Phi_params_update[[2]],
        mean_f = mean_f[[s]],
        var_f = var_f[[s]],
        mean_lambda = Lambda_params_update[[1]][[s]],
        var_lambda = Lambda_params_update[[2]][[s]],
        mean_l = mean_l[[s]],
        var_l = var_l[[s]]
      )
    })
    
    mean_Psi_update <- lapply(Psi_params_update, '[[', 3)
    
    factors_update <- update_factors(
      S = S,
      sizes = N_s,
      X_samp = X_scaled,
      mean_f_samp = mean_f,
      mean_l_samp = mean_l,
      P = P,
      K = K,
      J_s = J_s,
      mean_Psi = mean_Psi_update,
      mean_Phi = Phi_params_update[[1]],
      var_Phi = Phi_params_update[[2]],
      mean_Lambda = Lambda_params_update[[1]],
      var_Lambda = Lambda_params_update[[2]],
      n_iter = n_iter
    )
    
    mean_l_update <- lapply(factors_update, '[[', 1)
    var_l_update <-  lapply(factors_update, '[[', 2)
    mean_f_update <-  lapply(factors_update, '[[', 3)
    var_f_update <-  lapply(factors_update, '[[', 4)
    
    dist_omega_lambda_update <- update_omega_lambda(
      S = S,
      P = P,
      J = J_s,
      v_s = prior_omega,
      mean_lambda = Lambda_params_update[[1]],
      var_lambda = Lambda_params_update[[2]],
      dist_delta = dist_delta_lambda
    )
    
    dist_omega_phi_update <- update_omega_phi(
      P = P,
      K = K,
      v_s = prior_omega,
      mean_Phi = Phi_params_update[[1]],
      var_Phi = Phi_params_update[[2]],
      dist_delta = dist_delta_phi
    )
    
    dist_delta_phi_update <- update_shrinkage(
      P = P,
      dim = K,
      delta =  dist_delta_phi,
      prior_delta = dist_delta_phi_prior,
      omega = dist_omega_phi_update,
      mean_loading = Phi_params_update[[1]],
      var_loading = Phi_params_update[[2]]
    )
    
    dist_delta_lambda_update <- lapply(1:S, function(s) {
      update_shrinkage(
        P = P,
        dim = J_s[s],
        delta = dist_delta_lambda[[s]],
        prior_delta =  dist_delta_lambda_prior[[s]],
        omega = dist_omega_lambda_update[[s]],
        mean_loading = Lambda_params_update[[1]][[s]],
        var_loading = Lambda_params_update[[2]][[s]]
      )
    })
    
    sq_diff <- function(a,b){sum((a-b)^2)}
    
    change_Phi <- sq_diff(
      c(
        unlist(Phi_params_update[[1]]),
        unlist(Phi_params_update[[2]])
      ),
      c(
        unlist(mean_Phi),
        unlist(var_Phi)
      )
    )
    
    change_Lambda <- sq_diff(
      c(unlist(Lambda_params_update[[1]]), unlist(Lambda_params_update[[2]])),
      c(unlist(mean_Lambda), unlist(var_Lambda))
    )
    
    change_Psi <- sq_diff(
      c(unlist(lapply(Psi_params_update, "[[", 2))),
      c(unlist(lapply(dist_Psi, "[[", 2)))
    )
    
    change_omega <- sq_diff(
      c(
        unlist(lapply(1:P, function(p) dist_omega_phi_update[[p]][2,])),
        unlist(lapply(1:S, function(s) lapply(1:P, function(p) dist_omega_lambda_update[[s]][[p]][2,])))
      ),
      c(
        unlist(lapply(1:P, function(p) dist_omega_phi[[p]][2,])),
        unlist(lapply(1:S, function(s) lapply(1:P, function(p) dist_omega_lambda[[s]][[p]][2,])))
      )
    )
    
    change_delta <- sq_diff(
      c(
        c(dist_delta_phi_update[2,]),
        unlist(lapply(1:S, function(s) dist_delta_lambda_update[[s]][2,]))
      ),
      c(
        c(dist_delta_phi[2,]),
        unlist(lapply(1:S, function(s) dist_delta_lambda[[s]][2,]))
      )
    )
    
    delta <- change_Phi +
      change_Lambda +
      change_Psi + 
      change_omega + 
      change_delta
    
    delta_trace[n_iter] <- delta
    
    convergence <- convergence_criteria(delta_trace, n_iter, min_iter)
    
    if(verbose>1){
      print(paste0(
        "iter: ", as.integer(n_iter), " ",
        "convergence: ", convergence, " ",
        "change phi: ", round(change_Phi), " ",
        "change lambda: ", round(change_Lambda), " ",
        "change psi: ", round(change_Psi), " ",
        "change omega: ", round(change_omega), " ",
        "change delta: ", round(change_delta)
        
      ))
    } else if (verbose==1){
      paste0("iteration ", n_iter, " finished")
    }
    
    mean_f <- mean_f_update
    var_f <- var_f_update
    mean_l <- mean_l_update
    var_l <- var_l_update
    
    mean_Phi <- Phi_params_update[[1]]
    var_Phi <- Phi_params_update[[2]]
    mean_Lambda <- Lambda_params_update[[1]]
    var_Lambda <- Lambda_params_update[[2]]
    
    dist_Psi <- Psi_params_update
    mean_Psi <- mean_Psi_update
    dist_omega_phi <- dist_omega_phi_update
    dist_omega_lambda <- dist_omega_lambda_update
    
    dist_delta_phi <- dist_delta_phi_update
    dist_delta_lambda <- dist_delta_lambda_update
    
    if((convergence < tol | delta < 5e-3) & n_iter >= min_iter){
      if(verbose>0){
        print(paste("Algorithm converged in", n_iter, "iterations."))
      }
      break
    }
    
  }
  
  # Estimate for Phi
  Phi_estimated <- matrix(unlist(mean_Phi),
                          nrow = P,
                          ncol = K,
                          byrow = TRUE)
  
  # Estimate for Lambda
  Lambda_estimated <- lapply(1:S,
                             function(s) {
                               lambda_estimated <- matrix(
                                 unlist(mean_Lambda[[s]]),
                                 nrow = P,
                                 ncol = J_s[s],
                                 byrow = TRUE
                               )
                             })
  
  # Estimate for Psi
  Psi_estimated <- lapply(1:S, function(s, psi) mean_Psi[[s]]^-1, psi=mean_Psi)
  
  Psi_shape <- lapply(dist_Psi, "[[", 1)
  Psi_rate <- lapply(dist_Psi, "[[", 2)
  
  
  ests <- list("mean_phi"= Phi_estimated,
               "var_phi"= var_Phi,
               "mean_lambda_s"= Lambda_estimated,
               "var_lambda_s"= var_Lambda,
               "mean_psi_s" = Psi_estimated,
               "shape_psi_s" = Psi_shape,
               "rate_psi_s" = Psi_rate
  )
  return(ests)
}
