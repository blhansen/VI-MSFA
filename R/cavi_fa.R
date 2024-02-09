
#' Coordinate Ascent Variational Inference for Bayesian Factor Analysis
#'
#' The function implements the CAVI algorithm described in Hansen et al. (2023).
#' The code will return the relevant parameters of the mean-field
#' approximation described in the main paper.
#'
#' @param X A data matrix, of dimensions \eqn{P \times N}
#' @param J Number of factors.
#' @param verbose Desired verbosity of output printed out each iteration. A value greater than 1 will report the change in parameters for each iteration, a value equal to 1 will report when each iteration is finished, a value equal to 0 will return no output.
#' @param tol Convergence criteria. Default is 1e-3.
#' @param min_iter Minimum number of iterations before returning approximation. Default is 5.
#' @param max_iter Maximum number of iterations before returning approximation. Default is 100.
#' @param scale Controls whether the data should be scaled prior to performing CAVI.
#' @return A list  containing the the relevant parameters of the mean-field approximation.  The components of the list are:
#' \item{\code{mean_lambda}}{The mean of the factor loadings. A matrix of dimension \eqn{P \times J}.}
#' \item{\code{var_lambda}}{The variance of thefactor loadings. A list of length \eqn{P}, each entry is a matrix of dimension \eqn{K \times j} corresponding to the covariance of the common factor loadings of row \eqn{p}.}
#' \item{\code{shape_psi}}{The shape parameter of the gamma distribution about the inverse idiosyncratic errors. A vector of length \eqn{P}.}
#' \item{\code{rate_psi}}{The rate parameter of the gamma distribution about the inverse idiosyncratic errors. A vector of length \eqn{P}.}
#' \item{\code{shape_psi}}{The mean of the inverse of Psi, corresponding to idiosyncratic error. A vector of length \eqn{P}.}
#' \item{\code{mean_l}}{The mean of the study-specific latent factors associated with the common factor loadings. A list of length \eqn{S}, each entry is a list of length \eqn{N_s} where each entry is a vector of length \eqn{J}.}
#' @import sparsepca mnormt
#' @importFrom stats cov sd
#' @export
#' @references Hansen, B., Avalos-Pacheco, A., Massimiliano, R., De Vito, R. (2023).
#' Fast Variational Inference for Bayesian Factor Analysis in Single and Multi-Study Settings.
#' Submitted manuscript.
cavi_fa <- function(X,
                    J,
                    verbose = 1,
                    tol = 1e-2,
                    min_iter = 5,
                    max_iter = 100,
                    scale = TRUE) {

  N <- dim(X)[1]
  P <- dim(X)[2]

  if(scale){
    X_scale <- scale(X, center = TRUE, scale = TRUE)
  } else {
    X_scale <- scale(X, center = TRUE, scale = FALSE)
  }

  # Update Priors ###########################################
  update_omega_phi <- function(P, J, v_s, mean_Phi, var_Phi, dist_delta){
    tau <- cumprod(dist_delta[3,])
    lapply(1:P, function(p){sapply(1:J, function(j) {
      alpha <- (3 + 1)
      beta <- (3 + tau[j]*(mean_Phi[[p]][j]^2 + var_Phi[[p]][j,j]))
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
      beta <- prior_delta[2,l] + 1/2 * sum(sapply(l:dim, function(j, tau){
        tau[j] * sum(sapply(1:P, function(p){
          omega[[p]][3,j] * (mean_loading[[p]][j]^2 + var_loading[[p]][j,j])
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
  update_f <- function(N, X, P, J, mean_Psi, mean_Phi, var_Phi){
    phi_mat <-  matrix(unlist(mean_Phi), nrow=P, ncol=J, byrow=TRUE)

    var_f <- solve(
      diag(J) + Reduce('+', lapply(1:P, function(p) mean_Psi[p] * (tcrossprod(mean_Phi[[p]]) + var_Phi[[p]])))
    )

    # Find conditional mean for all obs within study s
    mean_fi <- lapply(1:N, function(i){
      var_f %*% t(phi_mat) %*% diag(mean_Psi) %*% matrix(X[i,], ncol=1)
    })

    return(list(mean_fi, var_f))
  }

  # Update Factor Loadings  ################################

  update_Phi <- function(P, X, N, J, dist_omega_phi, dist_delta, mean_Psi, mean_f, var_f){

    tau <- cumprod(dist_delta[3,])

    # Calclate variance for each lambda_p
    var_phi <- lapply(1:P, function(p){
      Sum_matrix <- mean_Psi[p] * Reduce('+', lapply(1:N, function(n) var_f + tcrossprod(mean_f[[n]])))
      Dp_matrix <- diag(dist_omega_phi[[p]][3,] * tau)
      var_phi_p <- solve(Dp_matrix + Sum_matrix)
    })

    mean_phi <- lapply(1:P, function(p){
      f <- matrix(unlist(mean_f), nrow=N, ncol=J, byrow=TRUE)
      x_p <- matrix(X[,p], ncol=1)
      diff_term <- mean_Psi[p] * (t(f) %*% x_p)

      matrix(var_phi[[p]] %*% diff_term, ncol=1)
    })

    return(list(mean_phi, var_phi))
  }

  # Update Errors  #########################################
  update_Psi <- function(x, n, p, Psi_prior, mean_Phi, var_Phi, mean_f, var_f){
    a_psi_update <- Psi_prior[1] + n/2

    b_psi_update <- Psi_prior[2] + 1/2*sapply(1:p, function(p){
      sq_sum <- sum(
        sapply(1:n, function(i){
          (x[i,p] - t(mean_Phi[[p]]) %*% mean_f[[i]])^2
        })
      )
      sum_latent_var <- Reduce("+", lapply(1:n, function(i){
        tcrossprod(mean_f[[i]])
      }))

      tr_sum <- sum(
        diag(
         (sum_latent_var + n * var_f) %*% var_Phi[[p]]
        )
      )

      quad_term <- t(mean_Phi[[p]]) %*% (n * var_f) %*% mean_Phi[[p]]

      sq_sum + tr_sum + quad_term
    })
    #return(a_psi_update/b_psi_update)
    return(list("shape"=a_psi_update, "rate"=b_psi_update))
  }

  # Helper Functions  ######################################
  convergence_criteria <- function(delta, n_iter, min_iter){
    if(n_iter<min_iter){res <- 1}else{
      res <- sqrt(delta[n_iter]) / (P*(J+1))
    }
    return(res)
  }

  start_vanilla_fa <- function(X, j, shrinkage_params = c(1e-5, 1e-5))
  {
    p <- dim(X)[2]

    spca_phi <- sparsepca::rspca(X, k=j, shrinkage_params[1], shrinkage_params[2])
    Phi <- as.matrix(spca_phi$loadings)
    psi <- diag(cov(X)-tcrossprod(Phi))
    scores <- spca_phi$scores

    for(col in 1:j){
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

  dist_omega_phi_prior <- lapply(1:P, function(p) matrix(c(rep(3/2, J), rep(3/2, J), rep(1, J)), nrow=3, ncol=J, byrow=TRUE))
  dist_omega_phi <- dist_omega_phi_prior

  a1 <- 2.1
  a2 <- 3.1

  dist_delta_common_prior <- matrix(
    c(c(a1, rep(a2, J-1)), rep(1, J), c(a1, rep(a2, J-1))), nrow=3, ncol=J, byrow=TRUE
  )
  dist_delta_common <- dist_delta_common_prior

  var_Phi_prior <- lapply(1:P, function(p) diag(dist_omega_phi[[p]][3,]^(-1)*cumprod(dist_delta_common[3,])^(-1)))
  var_Phi <- var_Phi_prior

  initial_params <- start_vanilla_fa(X_scale, J)
  mean_Phi <- lapply(1:P, function(x) matrix(initial_params$Phi[x,], ncol=1))
  mean_Psi <- initial_params$psi^-1

  mean_f <- lapply(1:N, function(i) matrix(initial_params$scores[i,], ncol=1) )
  var_f <- solve(
    diag(J) + Reduce('+', lapply(1:P, function(p) mean_Psi[p] * (tcrossprod(mean_Phi[[p]]) + var_Phi[[p]])))
  )

  initial_phi <- update_Phi(P = P,
                            X = X_scale,
                            N = N,
                            J = J,
                            dist_omega_phi = dist_omega_phi,
                            dist_delta = dist_delta_common,
                            mean_Psi = mean_Psi,
                            mean_f = mean_f,
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

    Phi_params_update <- update_Phi(P = P, X = X_scale, N = N, J = J,
                                    dist_omega_phi = dist_omega_phi, dist_delta = dist_delta_common, mean_Psi = mean_Psi,
                                    mean_f = mean_f, var_f = var_f
    )

    Psi_params_update <- update_Psi(x = X_scale, n = N, p = P, Psi_prior = prior_psi,
                                    mean_Phi = Phi_params_update[[1]], var_Phi = Phi_params_update[[2]],
                                    mean_f = mean_f, var_f = var_f
                                    )
    mean_Psi_update <- Psi_params_update$shape/Psi_params_update$rate

    f_params_update <- update_f(N = N, J = J, X = X_scale, P = P,
                                mean_Psi = mean_Psi_update, mean_Phi = Phi_params_update[[1]], var_Phi = Phi_params_update[[2]]
    )
    mean_f_update <- f_params_update[[1]]
    var_f_update <- f_params_update[[2]]

    dist_omega_phi_update <- update_omega_phi(P = P, J = J,
                                              v_s = prior_omega,
                                              mean_Phi = Phi_params_update[[1]], var_Phi = Phi_params_update[[2]],
                                              dist_delta = dist_delta_common
    )


    shrinkage_params_phi_update <- update_shrinkage(p = P, dim = J, delta =  dist_delta_common,
                                                    prior_delta = dist_delta_common_prior, omega = dist_omega_phi_update,
                                                    mean_loading = Phi_params_update[[1]], var_loading = Phi_params_update[[2]]
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

    change_Psi <- sq_diff(
      c((prior_psi[1] + N/2)/mean_Psi), # shape is known
      c(Psi_params_update$rate)
    )

    change_omega <- sq_diff(
      c(
        unlist(lapply(1:P, function(p) dist_omega_phi_update[[p]][2,]))
      ),
      c(
        unlist(lapply(1:P, function(p) dist_omega_phi[[p]][2,]))
      )
    )

    change_delta <- sq_diff(
      c(
        c(shrinkage_params_phi_update[2,])
      ),
      c(
        c(dist_delta_common[2,])
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


    if(convergence < tol & n_iter >= min_iter){
      if(verbose>0){
        print(paste("Algorithm converged in", n_iter, "iterations."))
      }
      break
    }


    mean_f <- mean_f_update
    var_f <- var_f_update

    mean_Phi <- Phi_params_update[[1]]
    var_Phi <- Phi_params_update[[2]]

    mean_Psi <- mean_Psi_update
    dist_omega_phi <- dist_omega_phi_update
    dist_delta_common <- shrinkage_params_phi_update

    print(paste0("iteration ", n_iter, " finished"))
  }

  # Estimate for Phi
  Phi_estimated <- matrix(unlist(mean_Phi),
                          nrow=P,
                          ncol=J,
                          byrow=TRUE
  )

  # Estimate for Psi
  Psi_estimated <- Psi_params_update$rate/(Psi_params_update$shape - 1)

  # Estimate for latent factors
  l_estimated <- do.call(cbind, mean_f)
  l_var <- var_f

  estimates <-
    list(
      "mean_lambda" = Phi_estimated,
      "var_lambda" = var_Phi,
      "shape_psi" = Psi_params_update$shape,
      "rate_psi" = Psi_params_update$rate,
      "mean_psi_inverse" = Psi_estimated,
      "mean_l" = l_estimated,
      "var_l" = l_var
    )

  return(estimates)
}
