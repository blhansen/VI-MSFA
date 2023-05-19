# VI-MSFA

author: Blake Hansen, Alejandra Avalos-Pacheco, Massimiliano Russo, Roberta De Vito

Approximate Bayesian Factor Analysis via [CAVI] and [SVI](#1-approximating-bayesian-FA).

## 1 Bayesian Factor Analysis via CAVI and SVI

The following example illustrates how to perform approximate inference for Bayesian Factor Analysis via CAVI and SVI,
using simulated data.

## Simulate some data

```{r sims, echo = TRUE, results = TRUE, tidy = TRUE}
P = 50
J = 4
N = 100

Lambda <- matrix(rnorm(P*J), nrow=P, ncol=J)
Psi <- runif(P, 0.1,1)
Sigma <- tcrossprod(Lambda) + diag(Psi)

X <- MASS::mvrnorm(N, mu=rep(0,P), Sigma=Sigma)
```

## Fitting Bayesian FA via CAVI
Now we can approximate the posterior distribution of Bayesian FA via CAVI and SVI

```{r get estimate, results = FALSE}
source("R/cavi_fa.R")
cavi_est <- cavi_fa(X, J, scale=FALSE)
```

The estimated parameter means can be accessed via
```{r analyze results, results = FALSE}
cavi_lambda <- cavi_est$mean_lambda
cavi_psi <- cavi_est$mean_psi
cavi_sigma <- tcrossprod(cavi_lambda) + diag(cavi_psi)

MatrixCorrelation::RV(cavi_sigma, Sigma)
```

We can do the same analysis using SVI instead:
```{r svi, results=FALSE}
source("R/svi_fa.R")

svi_est <- svi_fa(X, J, scale=FALSE)
svi_lambda <- svi_est$mean_lambda
svi_psi <- svi_est$mean_psi
svi_sigma <- tcrossprod(svi_lambda) + diag(svi_psi)

MatrixCorrelation::RV(svi_sigma, Sigma)
```

## 2 Bayesian Multistudy Factor Analysis via CAVI and SVI
We can perform the same analyses using multistudy methods via
```{r multistudy, results = FALSE}
S = 4
P = 50
K = 4
J_s = c(2,3,2,3)
N_s = c(50, 100, 60, 80)

Phi <- matrix(rnorm(P*K), nrow=P, ncol=K)
Lambda_s <- lapply(1:S, function(s){matrix(rnorm(P*J_s[s]), nrow=P, ncol=J_s[s])})
Psi_s <- lapply(1:S, function(s){runif(P, 0.1, 1)})
Sigma_s <- lapply(1:S, function(s){tcrossprod(Phi) + tcrossprod(Lambda_s[[s]]) + diag(Psi_s[[s]])})

X_s <- lapply(1:S, function(s) MASS::mvrnorm(N_s[s], mu=rep(0,P), Sigma=Sigma_s[[s]]))

source("R/cavi_msfa.R")

cavi_est <- cavi_msfa(X_s, K, J_s)
cavi_phi <- cavi_est$mean_phi
cavi_lambda_s <- cavi_est$mean_lambda_s
cavi_psi_s <- cavi_est$mean_psi_s
cavi_sigma_s <- lapply(1:S, function(s) tcrossprod(cavi_phi) + tcrossprod(cavi_lambda_s[[s]]) + diag(cavi_psi_s[[s]]))

MatrixCorrelation::RV(tcrossprod(cavi_phi), tcrossprod(Phi))

source("R/svi_msfa.R")

svi_est <- svi_msfa(X_s, K, J_s)
svi_phi <- svi_est$mean_phi
svi_lambda_s <- svi_est$mean_lambda_s
svi_psi_s <- svi_est$mean_psi_s
svi_sigma_s <- lapply(1:S, function(s) tcrossprod(svi_phi) + tcrossprod(svi_lambda_s[[s]]) + diag(svi_psi_s[[s]]))

MatrixCorrelation::RV(tcrossprod(svi_phi), tcrossprod(Phi))
```
