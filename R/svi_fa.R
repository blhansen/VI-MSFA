############################################################
# Variational Inference for standard Factor Analysis Model
# Author: Blake Hansen
# Spring 2022
############################################################

svi_fa <- function(X,
                   K,
                   batch_prop = 0.2,
                   fr = 0.75,
                   delay = 1,
                   verbose = 1,
                   tol = 1e-2,
                   min_iter = 5,
                   max_iter = 100,
                   scale = TRUE) {
  
  N <- dim(X)[1]
  P <- dim(X)[2]
  
  batch_size = floor(batch_prop*N)
  
  if(scale){
    X_scale <- scale(X, center=TRUE, scale=TRUE)
  } else {
    X_scale <- scale(X, center=TRUE, scale=FALSE)
  }
  
  # Update Priors ###########################################
  update_omega_phi <- function(P, K, v_s, mean_Phi, var_Phi, dist_delta){
    tau <- cumprod(dist_delta[3,])
    lapply(1:P, function(p){sapply(1:K, function(k) {
      alpha <- (3 + 1)/2
      beta <- (3 + tau[k]*(mean_Phi[[p]][k]^2 + var_Phi[[p]][k,k]))/2
      mean <- alpha/beta
      
      return(c(alpha, beta, mean))
    } )})
  }
  
  update_shrinkage <- function(p, dim, delta, prior_delta, omega, mean_loading, var_loading){
    delta <- lapply(1:dim, function(l){
      delta_l <- delta[3,]
      delta_l[l] <- 1
      
      tau <- cumprod(delta_l)
      
      alpha <- prior_delta[1,l] + 1/2 * p * (dim - l + 1)
      beta <- prior_delta[2,l] + 1/2 * sum(sapply(l:dim, function(k, tau){
        tau[k] * sum(sapply(1:P, function(p){
          omega[[p]][3,k] * (mean_loading[[p]][k]^2 + var_loading[[p]][k,k])
        }))
      },
      tau = tau
      )
      )
      
      delta[1,l] <- alpha
      delta[2,l] <- beta
      delta[3,l] <- alpha/beta
      
      
      mean_delta <- alpha/beta
      return(c(alpha, beta, mean_delta))
    })
    
    result <- matrix(unlist(delta), nrow=3, ncol=dim, byrow=FALSE)
    
    return(result)
  }
  
  # Update latent factors ###################################
  update_f <- function(N, X, P, K, mean_Psi, mean_Phi, var_Phi){
    phi_mat <-  matrix(unlist(mean_Phi), nrow=P, ncol=K, byrow=TRUE)
    
    var_fs <- solve(
      diag(K) + Reduce('+', lapply(1:P, function(p) mean_Psi[p] * (tcrossprod(mean_Phi[[p]]) + var_Phi[[p]])))
    )
    
    # Find conditional mean for all obs within study s
    mean_fis <- lapply(1:N, function(i){
      var_fs %*% t(phi_mat) %*% diag(mean_Psi) %*% matrix(X[i,], ncol=1)
    })
    
    return(list(mean_fis, var_fs))
  }
  
  # Update Factor Loadings  ################################
  
  update_Phi <- function(P, X, N, K, dist_omega_phi, dist_delta, mean_Psi, mean_f, var_f){
    
    tau <- cumprod(dist_delta[3,])
    
    # Caluclate variance for each lambda_ps
    var_phi <- lapply(1:P, function(p){
      Sum_matrix <- mean_Psi[p] * N/nrow(X)*Reduce('+', lapply(1:nrow(X), function(i) var_f + tcrossprod(mean_f[[i]])))
      
      diag_term <- dist_omega_phi[[p]][3,] * tau
      Dp_matrix <- diag(diag_term)
      var_phi_p <- mnormt::pd.solve(signif(Dp_matrix + Sum_matrix, digits=6))
    })
    
    mean_phi <- lapply(1:P, function(p){
      f <- matrix(unlist(mean_f), nrow=nrow(X), ncol=K, byrow=TRUE)
      x_p <- matrix(X[,p], ncol=1)
      diff_term <- N/nrow(X)*mean_Psi[p] * (t(f) %*% x_p)
      
      matrix(var_phi[[p]] %*% diff_term, ncol=1)
    })
    
    return(list(mean_phi, var_phi))
  }
  
  # Update Errors  #########################################
  
  # Calculate updates for alpha and beta parameters for Psi based on section 9 of Gibbs Sampling documentation
  update_Psi <- function(x, n, p, Psi_prior, mean_Phi, var_Phi, mean_f, var_f){
    a_psi_update <- Psi_prior[1] + n/2
    
    b_psi_update <- Psi_prior[2] + 1/2*sapply(1:p, function(p){
      sq_sum <- n/nrow(x)*sum(
        sapply(1:nrow(x), function(i){
          (x[i,p] - t(mean_Phi[[p]]) %*% mean_f[[i]])^2
        })
      )
      sum_latent_var <- Reduce("+", lapply(1:nrow(x), function(i){
        tcrossprod(mean_f[[i]])
      }))
      
      tr_sum <- sum(
        diag(
          n/nrow(x)*(sum_latent_var + var_f) %*% var_Phi[[p]] 
        )
      )
      
      quad_term <- t(mean_Phi[[p]]) %*% (n * var_f) %*% mean_Phi[[p]]
      
      result <- sq_sum + tr_sum + quad_term
      if(result<0){print(paste0("negative psi value for variable ", p))}
      return(result)
    })
    
    return(a_psi_update/b_psi_update)
  }
  
  # Helper Functions  ######################################
  convergence_criteria <- function(delta, n_iter, min_iter){
    if(n_iter<min_iter){res <- 1}else{
      res <- sqrt(delta[n_iter]) / (P*(K + 1) )
    }
    return(res)
  }
  
  start_vanilla_fa <- function(X, k, shrinkage_params = c(1e-5, 1e-5))
  {
    p <- dim(X)[2]
    
    spca_phi <- sparsepca::rspca(X, k=k, shrinkage_params[1], shrinkage_params[2])
    Phi <- as.matrix(spca_phi$loadings)
    psi <- diag(cov(X)-tcrossprod(Phi))
    scores <- spca_phi$scores
    
    for(col in 1:k){
      scaling_factor = sd(scores[,col])
      scores[,col] = scores[,col]/scaling_factor
      Phi[,col] = Phi[,col]*scaling_factor
    }
    
    
    out <- list(Phi=Phi, psi=psi, scores=scores)
    return(out)
  }
  
  # Initialize Parameters ##################################
  prior_omega <- c(3/2, 3/2)
  prior_psi <- c(1, 0.3)
  
  dist_omega_phi_prior <- lapply(1:P, function(p) matrix(c(rep(3/2, K), rep(3/2, K), rep(1, K)), nrow=3, ncol=K, byrow=TRUE))
  dist_omega_phi <- dist_omega_phi_prior
  
  a1 <- 2.1
  a2 <- 3.1
  
  dist_delta_common_prior <- matrix(
    c(c(a1, rep(a2, K-1)), rep(1, K), c(a1, rep(a2, K-1))), nrow=3, ncol=K, byrow=TRUE
  )
  dist_delta_common <- dist_delta_common_prior
  
  
  var_Phi_prior <- lapply(1:P, function(p) diag(dist_omega_phi[[p]][3,]^(-1)*cumprod(dist_delta_common[3,])^(-1)))
  var_Phi <- var_Phi_prior
  
  initial_params <- start_vanilla_fa(X_scale, K)
  mean_Phi <- lapply(1:P, function(x) matrix(initial_params$Phi[x,], ncol=1))
  mean_Psi <- initial_params$psi^(-1)
  Psi_params <- vector(mode="list", length=2)
  Psi_params[[1]] <- prior_psi[1] + N/2
  Psi_params[[2]] <- (prior_psi[1] + N/2) / mean_Psi 
  
  mean_f_store <- lapply(1:N, function(i) matrix(initial_params$scores[i,], ncol=1) )
  var_f <- solve(
    diag(K) + Reduce('+', lapply(1:P, function(p) mean_Psi[p] * (tcrossprod(mean_Phi[[p]]) + var_Phi[[p]])))
  )
  
  initial_phi <- update_Phi(P = P,
                            X = X_scale,
                            N = N,
                            K = K,
                            dist_omega_phi = dist_omega_phi,
                            dist_delta = dist_delta_common,
                            mean_Psi = mean_Psi,
                            mean_f = mean_f_store,
                            var_f = var_f
  )

  mean_Phi <- initial_phi[[1]]
  var_Phi <- initial_phi[[2]]
  
  # Perform CAVI  #########################################
  n_iter = 0
  
  delta_trace <- rep(NA, 1000)
  
  while(TRUE){
    n_iter <- n_iter + 1
    if(n_iter > max_iter){break}
    step_size <- (n_iter + delay)^-fr
    
    minibatch <- sample(1:N, size=batch_size)
    X_t <- X_scale[minibatch,]
    
    f_params_update <- update_f(N = batch_size,
                                K = K,
                                X = X_t,
                                P = P,
                                mean_Psi = mean_Psi,
                                mean_Phi = mean_Phi,
                                var_Phi = var_Phi
    )
    
    mean_f_batch <- f_params_update[[1]]
    var_f <- f_params_update[[2]]
    
    Phi_params_update <- update_Phi(P = P,
                                    X = X_t,
                                    N = N,
                                    K = K,
                                    dist_omega_phi = dist_omega_phi,
                                    dist_delta = dist_delta_common,
                                    mean_Psi = mean_Psi,
                                    mean_f = mean_f_batch,
                                    var_f = var_f
    )
    
    nat_phi_lam <- lapply(var_Phi, mnormt::pd.solve)
    nat_phi_nu <- lapply(1:P, function(p){
      nat_phi_lam[[p]] %*% mean_Phi[[p]]
    })
    
    nat_phi_lam_opt <- lapply(Phi_params_update[[2]],  mnormt::pd.solve)
    nat_phi_nu_opt <- lapply(1:P, function(p){
      nat_phi_lam_opt[[p]] %*% Phi_params_update[[1]][[p]]
    })
    
    nat_phi_lam_update <- lapply(1:P, function(p){
      step_size * nat_phi_lam_opt[[p]] + (1-step_size) * nat_phi_lam[[p]]
    })
    nat_phi_nu_update <- lapply(1:P, function(p){
      step_size * nat_phi_nu_opt[[p]] + (1-step_size) * nat_phi_nu[[p]]
    })
    
    var_Phi_update <- lapply(nat_phi_lam_update,  mnormt::pd.solve)
    mean_Phi_update <- lapply(1:P, function(p){
      matrix(var_Phi_update[[p]] %*% nat_phi_nu_update[[p]], ncol=1)
    })
    
    Psi_params_opt <- update_Psi(x = X_t,
                                 n = N,
                                 p = P,
                                 Psi_prior = prior_psi,
                                 mean_Phi = mean_Phi_update,
                                 var_Phi = var_Phi_update,
                                 mean_f = mean_f_batch,
                                 var_f = var_f
    )
    Psi_params_update <- vector(mode="list", length=2)
    Psi_params_update[[1]] <- Psi_params[[1]]
    Psi_params_update[[2]] <- step_size * Psi_params_opt + (1-step_size) * Psi_params[[2]]
    
    mean_Psi_update <- Psi_params_update[[1]]/Psi_params_update[[2]]
    
    dist_omega_phi_update <- update_omega_phi(P = P, 
                                              K = K,
                                              v_s = prior_omega,
                                              mean_Phi = mean_Phi_update,
                                              var_Phi = var_Phi_update,
                                              dist_delta = dist_delta_common
    )
    
    shrinkage_params_phi_update <- update_shrinkage(p = P,
                                                    dim = K,
                                                    delta =  dist_delta_common,
                                                    prior_delta = dist_delta_common_prior,
                                                    omega = dist_omega_phi_update,
                                                    mean_loading = mean_Phi_update,
                                                    var_loading = var_Phi_update
    )
    
    
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
    
    Psi_rate_new <- mean_Psi_update / (prior_psi[1] + N/2)
    Psi_rate_old <- mean_Psi / (prior_psi[1] + N/2)
    
    change_Psi <- sq_diff(
      c(Psi_rate_new),
      c(Psi_rate_old)
    )
    
    change_omega <- sq_diff(
      c(
        unlist(lapply(1:P, function(p) dist_omega_phi_update[[p]][3,]))
      ),
      c(
        unlist(lapply(1:P, function(p) dist_omega_phi[[p]][3,]))
      )
    )
    
    change_delta <- sq_diff(
      c(
        c(shrinkage_params_phi_update[3,])
      ),
      c(
        c(dist_delta_common[3,])
      )
    )
    
    delta <- change_Phi + change_Psi
    delta_trace[n_iter] <- delta
    convergence <- convergence_criteria(delta_trace, n_iter, min_iter)
    
    if(verbose>1){
      print(paste0(
        "iter: ", as.integer(n_iter), " ",
        "convergence: ", convergence, " ",
        "change phi: ", round(change_Phi), " ",
        "change psi: ", round(change_Psi), " ",
        "change omega: ", round(change_omega), " ",
        "change delta: ", round(change_delta)
        
      ))
    } else if (verbose==1){
      paste0("iteration ", n_iter, " finished")
    }
    
    
    
    if(convergence < tol){
      if(verbose>0){
        print(paste("Algorithm converged in", n_iter, "iterations."))
      }
      break
    }
    
    mean_Phi <- mean_Phi_update
    var_Phi <- var_Phi_update
    
    mean_Psi <- mean_Psi_update
    Psi_params <- Psi_params_update
    dist_omega_phi <- dist_omega_phi_update
    dist_delta_common <- shrinkage_params_phi_update
    
    print(paste0("iteration ", n_iter, " finished"))
  }
  
  # Estimate for Phi
  Phi_estimated <- matrix(unlist(mean_Phi),
                          nrow=P,
                          ncol=K,
                          byrow=TRUE
  )
  
  # Estimate for Psi
  Psi_estimated <- mean_Psi^(-1)
  
  estimates <- list("mean_lambda"=Phi_estimated, "mean_psi"=Psi_estimated, "n_iter"=n_iter, "var_lambda"=var_Phi)
  
  
  return(estimates)
}
