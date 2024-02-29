## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(
  fig.align = "center", 
  fig.width = 7, 
  fig.height = 7
)

## ----echo = FALSE, results = FALSE, warning = FALSE, message = FALSE----------
library(chouca)
library(plyr)
library(ggplot2)

## ----example_mod_def----------------------------------------------------------
kubo_model <- camodel(
  transition(from = "0",  to = "+", ~ alpha * p["+"]),         # reproduction of trees
  transition(from = "+", to = "0", ~ delta0 + delta * q["0"]), # death 
  parms = list(alpha = 0.6, delta0 = 0.05, delta = 0.4),
  wrap = TRUE,
  neighbors = 8
)

## -----------------------------------------------------------------------------
kubo_model

## -----------------------------------------------------------------------------
plot(kubo_model)

## ----initmat------------------------------------------------------------------
init_mat <- generate_initmat(kubo_model, 
                             pvec = c(0.5, 0.5),
                             nr = 128, nc = 256)

## ----run----------------------------------------------------------------------
out <- run_camodel(kubo_model, init_mat, times = seq(0, 100))

## ----display_console----------------------------------------------------------
summary(out)

## ----display_plot, fig.width = 12, fig.height = 6-----------------------------

oldpar <- par(mfrow = c(1, 2))
plot(out)
title("Global covers over time")
image(out)

# Restore graphical parameters
par(oldpar)

## -----------------------------------------------------------------------------
print(kubo_model)

## -----------------------------------------------------------------------------
mod <- camodel(
  transition(from = "dead", to = "live", ~ 0.1 + exp(0.2 * q["dead"] + p["dead"])),
  transition(from = "live", to = "dead", ~ 0.1),
  wrap = TRUE,
  neighbors = 8
)
print(mod)

## -----------------------------------------------------------------------------
mod <- camodel(transition(from = "dead", to = "live",
                          ~ 0.1 + sin(pi * q["dead"] * p["dead"] )),
               transition(from = "live", to = "dead",
                          ~ 0.1),
               wrap = TRUE,
               neighbors = 8)
print(mod)

## -----------------------------------------------------------------------------

sin_approx <- function(x) {
  x - x^3/factorial(3) + x^5/factorial(5)
}

mod <- camodel(transition(from = "dead", to = "live",
                          ~ 0.1 + sin_approx(pi * q["dead"] * p["dead"]) ),
               transition(from = "live", to = "dead",
                          ~ 0.1),
               wrap = TRUE,
               neighbors = 8)
print(mod) # negligible error this time, as we use an approximation of sin()


## -----------------------------------------------------------------------------
control_args <- list(engine = "compiled",
                     precompute_probas = FALSE) # see below for meaning of this parameter

out <- run_camodel(kubo_model,
                   initmat = init_mat,
                   times = seq(0, 512),
                   control = control_args)

## -----------------------------------------------------------------------------

control_args <- list(engine = "compiled",
                     precompute_probas = TRUE,
                     console_output_every = 128) # report less often on console

out <- run_camodel(kubo_model,
                   initmat = init_mat,
                   times = seq(0, 512),
                   control = control_args)


## -----------------------------------------------------------------------------

control_args <- list(engine = "compiled",
                     precompute_probas = TRUE,
                     cores = 2,
                     console_output_every = 128)

out <- run_camodel(kubo_model,
                   initmat = init_mat,
                   times = seq(0, 512),
                   control = control_args)


## -----------------------------------------------------------------------------

init <- generate_initmat(kubo_model, c(0.5, 0.5),
                         nr = 128)

run <- run_camodel(kubo_model, init, times = seq(0, 128, by = 1))

# Extract covers and display the last lines of the table
covers <- run[["output"]][["covers"]]
tail(covers)

## -----------------------------------------------------------------------------
landscapes <- run[["output"]][["snapshots"]]
spatialwarnings::display_matrix(landscapes)

## -----------------------------------------------------------------------------
# Save every eight iterations
run2 <- run_camodel(kubo_model, init, times = seq(0, 128, by = 8))
covers2 <- run2[["output"]][["covers"]]
tail(covers2)

plot(covers[ ,"t"], covers[ ,"+"], type = "l", col = "red")
lines(covers2[ ,"t"], covers2[ ,"+"], col = "blue")

## -----------------------------------------------------------------------------
ctrl <- list(save_covers_every = 8)
run3 <- run_camodel(kubo_model, init, times = seq(0, 128, by = 1),
                    control = ctrl)
covers3 <- run3[["output"]][["covers"]]
tail(covers3)

plot(covers[ ,"t"], covers[ ,"+"], type = "b", col = "red",
     xlab = "t", ylab = "forest c over")
lines(covers2[ ,"t"], covers2[ ,"+"], col = "blue")
points(covers2[ ,"t"], covers2[ ,"+"], col = "blue")
lines(covers3[ ,"t"], covers3[ ,"+"], col = "darkgreen")
points(covers3[ ,"t"], covers3[ ,"+"], col = "darkgreen")
legend(x = "right",
       legend = c("original", "reduced 'times' vector", "using 'save_covers_every'"),
       col = c("red", "blue", "darkgreen"),
       lty = c(1, 1, 1))


## -----------------------------------------------------------------------------

mod <- ca_library("forestgap")
init <- generate_initmat(mod, c(TREE = .5, EMPTY = .5), nr = 256)

# Define the function that computes AC at lag-1 using Moran's I
compute_aclag1 <- function(t, mat) {
  m <- matrix(mat == "TREE", nrow = nrow(mat), ncol = ncol(mat))
  data.frame(t = t, ac = spatialwarnings::raw_moran(mat))
}

control_args <- list(save_covers_every = 0,
                     custom_output_every = 1,
                     custom_output_fun = compute_aclag1)

out <- run_camodel(mod, init, times = seq(0, 64), control = control_args)


## -----------------------------------------------------------------------------
my_output <- out[["output"]][["custom"]]

## -----------------------------------------------------------------------------
my_tbl <- do.call(rbind, my_output)
tail(my_tbl, 10)

## ----trace_landscape, eval = FALSE--------------------------------------------
#  
#  mod <- ca_library("rock-paper-scissor")
#  
#  ctrl <- list(custom_output_every = 1,
#               precompute_probas = FALSE,
#               custom_output_fun = landscape_plotter(mod))
#  
#  init <- generate_initmat(mod, rep(1, 3)/3, 100, 178)
#  
#  out <- run_camodel(mod, init, seq(0, 128), control = ctrl)
#  

## ----trace_covers, eval = FALSE-----------------------------------------------
#  
#  ctrl <- list(custom_output_every = 1,
#               custom_output_fun = trace_plotter(mod, init, max_samples = 128),
#               console_output_every = 128)
#  
#  out <- run_camodel(mod, init, seq(0, 512), control = ctrl)
#  

