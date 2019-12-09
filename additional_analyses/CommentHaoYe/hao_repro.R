library(dplyr)
library(purrr)
library(ggplot2)
library(rEDM)

set.seed(42)

# functions
simulate_data <- function(burn_in = 500, n = 300, 
                          init_x = runif(1,0.1,0.7), 
                          init_y = runif(1,0.1,0.7))
{
    obs_idx <- burn_in + seq(n)
    
    x <- numeric(burn_in + n)
    y <- numeric(burn_in + n)
    x[1] <- init_x
    y[1] <- init_y
    for (t in 1:(burn_in + n - 1))
    {
        x[t + 1] <- x[t] * (3.8 - 3.8 * x[t])
        y[t + 1] <- y[t] * (3.5 - 3.5 * y[t])
    }
    
    x <- log(x[obs_idx])
    y <- log(y[obs_idx])
    x <- x - mean(x)
    y <- y - mean(y)
    data.frame(x, y)
}

do_ccm <- function(dat, lib_sizes = 100, ...)
{
    bind_rows(ccm(dat, E = 1, lib_column = 1, target_column = 2, 
                  replace = FALSE, silent = TRUE, 
                  lib_sizes = lib_sizes, ...) %>% 
                  ccm_means() %>%
                  mutate(dir = "x -> y"), 
              ccm(dat, E = 1, lib_column = 2, target_column = 1, 
                  replace = FALSE, silent = TRUE, 
                  lib_sizes = lib_sizes, ...) %>%
                  ccm_means() %>%
                  mutate(dir = "y -> x"))
}

init <- readRDS("init_values.RDS")
dat <- simulate_data(init_x = init$x, init_y = init$y)
out <- do_ccm(dat)

# random shuffle surrogates
x_surr <- make_surrogate_data(dat$x)
out_surr <- map_dfr(seq_len(NCOL(x_surr)), function(j) {
    temp_dat <- data.frame(x = x_surr[, j], y = dat$y)
    do_ccm(temp_dat)
})

# random shuffle, maintaining periodicity
make_periodic_surrogate <- function(ts, num_surr = 100, period = 2)
{
    n <- length(ts)
    matrix(unlist(
        lapply(seq(num_surr), function(i) {
            surr_ts <- numeric(n)
            for (phase in seq_len(period))
            {
                idx <- seq(from = phase, to = n, by = period)
                surr_ts[idx] <- ts[sample(idx)]
            }
            surr_ts
        })
    ), ncol = num_surr)
}
x_surr_periodic <- make_periodic_surrogate(dat$x, period = 8)
out_surr_periodic <- map_dfr(seq_len(NCOL(x_surr_periodic)), function(j) {
    temp_dat <- data.frame(x = x_surr_periodic[, j], y = dat$y)
    do_ccm(temp_dat)
})
