## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE-----------------------------------------------------------
library(mpmsim)
library(dplyr)
library(Rage)

## -----------------------------------------------------------------------------
b_1_values <- seq(0.1, 0.9, 0.1)
lifeTables <- list()
for (i in 1:length(b_1_values)) {
  lifeTables[[i]] <- model_survival(
    params = c(b_0 = 0.1, b_1 = b_1_values[i]),
    model = "Gompertz"
  )
}

## -----------------------------------------------------------------------------
for (i in 1:length(lifeTables)) {
  lifeTables[[i]] <- lifeTables[[i]] |>
    mutate(stage = ifelse(x <= round(max(x) * 0.25), 1,
      ifelse(x <= round(max(x) * 0.75), 2, 3)
    ))
}

## -----------------------------------------------------------------------------
for (i in 1:length(lifeTables)) {
  lifeTables[[i]] <- lifeTables[[i]] |>
    mutate(fert = model_fertility(
      age = x, params = c(A = 3),
      maturity = min(x[stage == 2]),
      model = "step"
    ))
}

## -----------------------------------------------------------------------------
lifeTables[[5]]

## -----------------------------------------------------------------------------
leslie_matrices <- list()
for (i in 1:length(lifeTables)) {
  leslie_matrices[[i]] <- make_leslie_mpm(
    survival = lifeTables[[i]]$px,
    fertility = lifeTables[[i]]$fert,
    n_stages = nrow(lifeTables[[i]]), split = TRUE
  )
}

## -----------------------------------------------------------------------------
leslie_matrices[[5]]$matA

## -----------------------------------------------------------------------------
collapsed_matrices <- list()
for (i in 1:length(lifeTables)) {
  stages <- lifeTables[[i]]$stage
  matrices <- leslie_matrices[[i]]
  collapse_list <- split(stages, stages)
  # get the indices of each element in the original vector
  collapse_list <- lapply(collapse_list, function(x) which(stages %in% x))

  collapsed_matrices[[i]] <- Rage::mpm_collapse(
    matU = matrices$mat_U,
    matF = matrices$mat_F, collapse = collapse_list
  )
}

## -----------------------------------------------------------------------------
collapsed_matrices[[5]]$matA

## -----------------------------------------------------------------------------
recovered_life_tables <- list()
for (i in 1:length(lifeTables)) {
  m1 <- collapsed_matrices[[i]]
  recovered_life_tables[[i]] <- Rage::mpm_to_table(
    matU = m1$matU, matF = m1$matF,
    remove_final = TRUE
  )
}

## ---- fig.height = 4, fig.width = 6, fig.align = "center"---------------------
i <- 5
recovered_lt <- recovered_life_tables[[i]]

plot(0:(length(recovered_lt$lx) - 1),
  recovered_lt$lx,
  type = "l",
  xlab = "age", ylab = "survivorship", main = "Survivorship"
)
lines(lifeTables[[i]]$x, lifeTables[[i]]$lx, type = "l", col = "red")

## ---- fig.height = 4, fig.width = 6, fig.align = "center"---------------------
plot(0:(length(recovered_lt$lx) - 1),
  recovered_lt$hx,
  type = "l",
  xlab = "age", ylab = "mortality", main = "mortality", ylim = c(0, 2)
)
lines(lifeTables[[i]]$x, lifeTables[[i]]$hx, type = "l", col = "red")

## ---- fig.height = 4, fig.width = 6, fig.align = "center"---------------------
df1 <- data.frame(b_1_values, lifespan_lt = NA, lifespan_afs = NA)
for (i in 1:length(lifeTables)) {
  df1$lifespan_lt[i] <- max(lifeTables[[i]]$x)
  df1$lifespan_afs[i] <- max(recovered_life_tables[[i]]$x)
}

df1 <- df1 %>%
  mutate(lifespan_diff = lifespan_afs - lifespan_lt) %>%
  mutate(lifespan_diff_perc = 100 * (lifespan_diff / lifespan_lt))

plot(df1$b_1_values, df1$lifespan_diff_perc,
  type = "b",
  ylab = "Lifespan overestimation (%)", xlab = "Gompertz parameter"
)

