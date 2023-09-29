## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(42)

## ----message=FALSE------------------------------------------------------------
library(mpmsim)
library(Rage)
library(Rcompadre)
library(dplyr)
library(popbio)
library(ggfortify)

## -----------------------------------------------------------------------------
set.seed(42)
constrain_df <- data.frame(fun = "lambda", arg = NA, lower = 0.9, upper = 1.1)
sim_life_hist_1 <- generate_mpm_set(
  n = 50, n_stages = 3, fecundity = c(0, 6, 6), archetype = 1, split = TRUE,
  max_surv = 0.95, constraint = constrain_df
)

## -----------------------------------------------------------------------------
sim_life_hist_1 <- cdb_flag(sim_life_hist_1, checks = "check_irreducible") %>%
  filter(check_irreducible == TRUE)

## -----------------------------------------------------------------------------
# Put the matrices into the metadata
sim_life_hist_1$matA <- matA(sim_life_hist_1)
sim_life_hist_1$matU <- matU(sim_life_hist_1)
sim_life_hist_1$matF <- matF(sim_life_hist_1)

# Use cdb_metadata to turn this into a data frame
sim_life_hist_1 <- cdb_metadata(sim_life_hist_1)

## -----------------------------------------------------------------------------
# New functions to calculate generation time from life table.
# Function to calculate generation time from the life table
gt_lt <- function(matU, matF, start = 1, ...) {
  tempLT <- mpm_to_table(matU, matF, start = start, ...)
  return(sum(tempLT$x * tempLT$lxmx) / sum(tempLT$lxmx))
}

## -----------------------------------------------------------------------------
sim_life_hist_1$gt_lt <- mapply(
  gt_lt, sim_life_hist_1$matU,
  sim_life_hist_1$matF
)
sim_life_hist_1$longevity <- sapply(sim_life_hist_1$matU,
  Rage::longevity,
  x_max = 1000, lx_crit = 0.01
)
sim_life_hist_1$lifeExpect <- sapply(
  sim_life_hist_1$matU,
  Rage::life_expect_mean
)
sim_life_hist_1$entropy_d <- mapply(
  entropy_d,
  sim_life_hist_1$matU,
  sim_life_hist_1$matF
)

sim_life_hist_1$entropy_k <- mapply(entropy_k, sim_life_hist_1$matU)
sim_life_hist_1$nrr_R0 <- mapply(
  net_repro_rate, sim_life_hist_1$matU,
  sim_life_hist_1$matF
)

## -----------------------------------------------------------------------------
pcData <- sim_life_hist_1 %>%
  select(gt_lt, longevity, lifeExpect, entropy_d, entropy_k, nrr_R0) %>%
  na.omit()

## -----------------------------------------------------------------------------
PCA <- prcomp(pcData, scale = TRUE, center = TRUE)

# Add the PC data to the raw data.
pcData <- pcData %>%
  cbind(PCA$x[, 1:2])

## ----fig.height = 4, fig.width = 6, fig.align = "center", warning=FALSE-------
PCA_plot <- autoplot(
  object = PCA, alpha = 0, size = 4, fill = "#55616D60",
  loadings.colour = "#0072B2", shape = 16,
  loadings = TRUE, loadings.label = TRUE, loadings.label.colour = "red",
  loadings.label.size = 3, loadings.label.repel = TRUE,
  frame = FALSE, frame.type = "norm", scale = 0
)

PCA_plot$layers <- c(
  geom_point(
    aes_(
      x = pcData$PC1,
      y = pcData$PC2
    ),
    size = 2, alpha = 0.5
  ),
  PCA_plot$layers
)

PCA_plot

## ----fig.height = 4, fig.width = 6, fig.align = "center", warning=FALSE, echo = FALSE----
## -----------------------------------------------------------------------------
set.seed(42)
constrain_df <- data.frame(fun = "lambda", arg = NA, lower = 0.9, upper = 1.1)

sim_life_hist_1 <- generate_mpm_set(
  n = 50, n_stages = 3, fecundity = c(0, 6, 6), archetype = 4,
  split = TRUE, constraint = constrain_df
)


## -----------------------------------------------------------------------------
sim_life_hist_1 <- cdb_flag(sim_life_hist_1, checks = "check_irreducible") %>%
  filter(check_irreducible == TRUE)

## -----------------------------------------------------------------------------
# Put the matrices into the metadata
sim_life_hist_1$matA <- matA(sim_life_hist_1)
sim_life_hist_1$matU <- matU(sim_life_hist_1)
sim_life_hist_1$matF <- matF(sim_life_hist_1)

# Use cdb_metadata to turn this into a data frame
sim_life_hist_1 <- cdb_metadata(sim_life_hist_1)

## -----------------------------------------------------------------------------
# New functions to calculate generation time from life table.
# Function to calculate generation time from the life table
gt_lt <- function(matU, matF, start = 1, ...) {
  tempLT <- mpm_to_table(matU, matF, start = start, ...)
  return(sum(tempLT$x * tempLT$lxmx) / sum(tempLT$lxmx))
}

## -----------------------------------------------------------------------------
sim_life_hist_1$gt <- sapply(
  sim_life_hist_1$matA,
  popbio::generation.time
)
sim_life_hist_1$gt_lt <- mapply(
  gt_lt, sim_life_hist_1$matU,
  sim_life_hist_1$matF
)
sim_life_hist_1$longevity <- sapply(sim_life_hist_1$matU, Rage::longevity,
  x_max = 1000, lx_crit = 0.01
)
sim_life_hist_1$lifeExpect <- sapply(
  sim_life_hist_1$matU,
  Rage::life_expect_mean
)
sim_life_hist_1$entropy_d <- mapply(
  entropy_d, sim_life_hist_1$matU,
  sim_life_hist_1$matF
)

sim_life_hist_1$entropy_k <- mapply(entropy_k, sim_life_hist_1$matU)
sim_life_hist_1$nrr_R0 <- mapply(
  net_repro_rate, sim_life_hist_1$matU,
  sim_life_hist_1$matF
)
sim_life_hist_1$gt_Rage <- mapply(
  gen_time, sim_life_hist_1$matU,
  sim_life_hist_1$matF
)

## -----------------------------------------------------------------------------
pcData <- sim_life_hist_1 %>%
  select(gt_lt, longevity, lifeExpect, entropy_d, entropy_k, nrr_R0) %>%
  na.omit()

## -----------------------------------------------------------------------------
PCA <- prcomp(pcData, scale = TRUE, center = TRUE)

# Add the PC data to the raw data.
pcData <- pcData %>%
  cbind(PCA$x[, 1:2])

## -----------------------------------------------------------------------------
PCA_plot <- autoplot(
  object = PCA, alpha = 0, size = 4, fill = "#55616D60",
  loadings.colour = "#0072B2", shape = 16,
  loadings = TRUE, loadings.label = TRUE, loadings.label.colour = "black",
  loadings.label.size = 3, loadings.label.repel = TRUE,
  frame = FALSE, frame.type = "norm", scale = 0
)

PCA_plot$layers <- c(
  geom_point(
    aes_(
      x = pcData$PC1,
      y = pcData$PC2
    ),
    size = 2, alpha = 0.5
  ),
  PCA_plot$layers
)

PCA_plot

