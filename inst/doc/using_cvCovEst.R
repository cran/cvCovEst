## ----brief-example------------------------------------------------------------
library(MASS)
library(cvCovEst)
set.seed(1584)

# generate a 50x50 covariance matrix with unit variances and off-diagonal
# elements equal to 0.5
sigma <- matrix(0.5, nrow = 50, ncol = 50) + diag(0.5, nrow = 50)

# sample 50 observations from multivariate normal with mean = 0, var = Sigma
dat <- mvrnorm(n = 50, mu = rep(0, 50), Sigma = sigma)

# run CV-selector
cv_cov_est_out <- cvCovEst(
    dat = dat,
    estimators = c(
      linearShrinkLWEst, denseLinearShrinkEst,
      thresholdingEst, poetEst, sampleCovEst
    ),
    estimator_params = list(
      thresholdingEst = list(gamma = c(0.2, 0.4)),
      poetEst = list(lambda = c(0.1, 0.2), k = c(1L, 2L))
    ),
    cv_loss = cvMatrixFrobeniusLoss,
    cv_scheme = "v_fold",
    v_folds = 5,
  )

# print the table of risk estimates
cv_cov_est_out$risk_df

# print a subset of the selected estimator's estimate
cv_cov_est_out$estimate[1:5, 1:5]

## ----toep-sim-----------------------------------------------------------------
library(MASS)
library(cvCovEst)
set.seed(1584)

toep_sim <- function(p, rho, alpha) {
    times <- seq_len(p)
    H <- abs(outer(times, times, "-")) + diag(p)
    H <- H^-(1 + alpha) * rho
    covmat <- H + diag(p) * (1 - rho)

    sign_mat <- sapply(
      times,
      function(i) {
        sapply(
          times,
          function(j) {
            (-1)^(abs(i - j))
          }
        )
      }
    )
    return(covmat * sign_mat)
}

# simulate a 100 x 100 covariance matrix
sim_covmat <- toep_sim(p = 100, rho = 0.6, alpha = 0.3)

# sample 75 observations from multivariate normal mean = 0, var = sim_covmat
sim_dat <-  MASS::mvrnorm(n = 75, mu = rep(0, 100), Sigma = sim_covmat)

# run CV-selector
cv_cov_est_sim <- cvCovEst(
  dat = sim_dat,
  estimators = c(
    linearShrinkEst, thresholdingEst,
    bandingEst, adaptiveLassoEst,
    sampleCovEst),
  estimator_params = list(
    linearShrinkEst = list(alpha = seq(0.25, 0.75, 0.05)),
    thresholdingEst = list(gamma = seq(0.25, 0.75, 0.05)),
    bandingEst = list(k = seq(2L, 10L, 2L)),
    adaptiveLassoEst = list(lambda = c(0.1, 0.25, 0.5, 0.75, 1),
                            n = seq(1, 5))),
  cv_loss = cvMatrixFrobeniusLoss,
  cv_scheme = "v_fold",
  v_folds = 5,
  center = TRUE,
  scale = TRUE
)

## ----summary-sim--------------------------------------------------------------
cv_sum <- summary(cv_cov_est_sim)
cv_sum$bestInClass

## ----hyperRisk-sim------------------------------------------------------------
cv_sum$hyperRisk$bandingEst

## ----risk-sim, fig.height = 3.5, fig.width=4----------------------------------
plot(cv_cov_est_sim, dat_orig = sim_dat, plot_type = "risk")

## ----multi-heat-sim, fig.height=3.5, fig.width=7------------------------------
plot(cv_cov_est_sim, dat_orig = sim_dat, plot_type = "heatmap",
     stat = c("min", "median", "max"))

## ----sign-multi-heat-sim, fig.height=3.5, fig.width=7-------------------------
plot(cv_cov_est_sim, dat_orig = sim_dat, plot_type = "heatmap",
     stat = c("min", "median", "max"), abs_v = FALSE)

## ----samp-heat-sim, fig.height=3.5, fig.width=5-------------------------------
plot(cv_cov_est_sim, dat_orig = sim_dat, plot_type = "heatmap",
     estimator = c("bandingEst", "sampleCovEst"),
     stat = c("min"), abs_v = FALSE)

## ----eigen-sim, fig.height = 3.5, fig.width=5---------------------------------
plot(cv_cov_est_sim, dat_orig = sim_dat, plot_type = "eigen",
     stat = c("min", "median", "max"))

## ----multi-eigen-sim, fig.height=4.5, fig.width=6-----------------------------
plot(cv_cov_est_sim, dat_orig = sim_dat, plot_type = "eigen",
     stat = c("min", "median", "max"),
     estimator = c("bandingEst", "thresholdingEst", "linearShrinkEst",
                    "adaptiveLassoEst"))

## ----plot-summary, fig.height=7, fig.width=7.5--------------------------------
plot(cv_cov_est_sim, dat_orig = sim_dat)

