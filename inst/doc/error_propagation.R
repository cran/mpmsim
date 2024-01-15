## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(mpmsim)
set.seed(42)

## -----------------------------------------------------------------------------
matU <- matrix(c(
  0.1, 0.0,
  0.2, 0.4
), byrow = TRUE, nrow = 2)

matF <- matrix(c(
  0.0, 3.0,
  0.0, 0.0
), byrow = TRUE, nrow = 2)

## -----------------------------------------------------------------------------
compute_ci(
  mat_U = matU, mat_F = matF, sample_size = 20,
  FUN = popdemo::eigs, what = "lambda"
)

## ---- fig.height = 4, fig.width = 6, fig.align = "center"---------------------
distLambda_20 <- compute_ci(
  mat_U = matU, mat_F = matF,
  sample_size = 20, FUN = popdemo::eigs, what = "lambda",
  dist.out = TRUE
)
hist(distLambda_20$estimates)

## ---- fig.height = 4, fig.width = 6, fig.align = "center"---------------------
distLambda_100 <- compute_ci(
  mat_U = matU, mat_F = matF,
  sample_size = 100, FUN = popdemo::eigs, what = "lambda",
  dist.out = TRUE
)

## ----compareDist2, fig.height = 6, fig.width = 6, fig.align = "center"--------
par(mfrow = c(2, 1))
hist(distLambda_20$estimates, xlim = c(0, 1.75))
hist(distLambda_100$estimates, xlim = c(0, 1.75))

## -----------------------------------------------------------------------------
observedLambda <- popdemo::eigs(matU + matF, what = "lambda")
reducedLambda <- observedLambda * 0.8


simDist <- compute_ci(
  mat_U = matU, mat_F = matF,
  sample_size = 50, FUN = popdemo::eigs, what = "lambda",
  dist.out = TRUE
)$estimates
hist(simDist)
abline(v = observedLambda, lty = 2, lwd = 2)
abline(v = reducedLambda, lty = 2, lwd = 2, col = "red")
sum(simDist < observedLambda) / length(simDist)

## -----------------------------------------------------------------------------
sample_size_mat <- matrix(c(
  20, 100,
  20, 20
), byrow = TRUE, nrow = 2)

distLambda_variable <- compute_ci(
  mat_U = matU, mat_F = matF,
  sample_size = sample_size_mat,
  FUN = popdemo::eigs, what = "lambda",
  dist.out = TRUE
)
hist(distLambda_variable$estimates)

## ----compareDist3, fig.height = 9, fig.width = 6, fig.align = "center"--------
par(mfrow = c(3, 1))
hist(distLambda_20$estimates, xlim = c(0, 1.75))
hist(distLambda_100$estimates, xlim = c(0, 1.75))
hist(distLambda_variable$estimates, xlim = c(0, 1.75))

