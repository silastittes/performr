#' Pipe
#'
#' Internal access of the pipe from magrittr
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @param lhs,rhs specify what lhs and rhs are
NULL



#' Performance data frame
#'
#' Constructs a tidy data frame from the output of stan_performance()
#' @importFrom purrr map map_lgl map_dbl map_df
#' @importFrom magrittr set_colnames
#' @importFrom tidyr gather
#' @importFrom dplyr group_by mutate select summarise summarise_all
#' @param stan_out The name of the object output of stan_performance().
#' @param species_order A character vector of IDs that match the output dimensions of each parameter matrix.
#' @return A tidy data frame, where each row is a draw from one of the input groups and each column is a model parameter. Additionally, several derived parameters, including the optimum, area, breadth, and area scaled by breadth (called special).
#' @export

perform_df <- function(stan_out, species_order){

  #stan_out output from gen_tolerance()
  #species_order vector of desired labels matching group_ids used in gen_tolerance()

  new_post <- stan_out %>% rstan::extract()

  params <- c("x_min", "x_max", "shape1", "shape2", "stretch", "nu")

  par_df <- params %>% map( ~{
    new_post[[.x]] %>%
      data.frame() %>%
      set_colnames(species_order) %>%
      gather("species", x) %>%
      set_colnames(c("species", .x)) %>%
      group_by(species) %>%
      mutate(draw = 1:n())
  }) %>%
    do.call(cbind, .) %>%
    select(-starts_with("species")) %>%
    select(-matches("draw[0-9]")) %>%
    mutate(
      maxima = (((shape1 - 1)/(shape1*shape2 - 1))^(1/shape1) * (x_max - x_min) + x_min),
      breadth = (x_max - x_min),
      area = stretch * breadth,
      special = area/breadth
    )

  par_df
}



#' Performance curve data frame
#'
#' Constructs a tidy data frame representing the full perforance curve (facilitates plotting).
#' @importFrom purrr map map_lgl map_dbl map_df
#' @importFrom magrittr set_colnames
#' @importFrom tidyr gather
#' @importFrom dplyr group_by mutate select summarise summarise_all
#' @param x A [0,1] vector of values to evaluate performance over
#' @param par_df A data frame produced by perform_df().
#' @return A tidy data frame, containing a points along the environmental axis, and corresponding points for the performance axis, for each group and posterior draw input.
#' @export

map_performance <- function(x = seq(0, 1, length.out = 100), par_df){
  1:nrow(par_df) %>% map_df(~{
    draw_x <- par_df$draw[.x]
    species <- par_df$species[.x]
    x_min <- par_df$x_min[.x]
    x_max <- par_df$x_max[.x]
    shape1 <- par_df$shape1[.x]
    shape2 <- par_df$shape2[.x]
    stretch <- par_df$stretch[.x]
    xs <- x*(x_max - x_min) + x_min
    mod_fit <- stretch*((shape1*shape2*x^(shape1-1)) * (1-x^shape1)^(shape2-1))
    data_frame(species = species, x = xs, y = mod_fit, draw = draw_x)
  })
}





#' Performance curve data frame
#'
#' Similar to map_performance(), but over a set of fixed points along the axis rather than being sampled relative to each species performance limits.
#' @importFrom purrr map map_lgl map_dbl map_df
#' @importFrom magrittr set_colnames
#' @importFrom tidyr gather
#' @importFrom dplyr group_by mutate select summarise summarise_all
#' @param x A vector of values to evaluate performance over.
#' @param par_df A data frame produced by perform_df().
#' @return A tidy data frame, containing a points along the environmental axis, and corresponding points for the performance axis, for each group and posterior draw input.
#' @export

map_performance_fixed <- function(x, par_df){
  1:nrow(par_df) %>% map_df(~{
    draw_x <- par_df$draw[.x]
    species <- par_df$species[.x]
    x_min <- par_df$x_min[.x]
    x_max <- par_df$x_max[.x]
    shape1 <- par_df$shape1[.x]
    shape2 <- par_df$shape2[.x]
    stretch <- par_df$stretch[.x]
    mod_fit <-performance_mu(x, shape1, shape2, stretch, x_min, x_max)
    data_frame(species = species, x = x, y = mod_fit, draw = draw_x)
  })
}



#' Generate psuedo-observed data from performance curve parameters
#'
#' Constructs a tidy data frame of generated data given a set of input parameters.
#' @importFrom purrr map map_lgl map_dbl map_df
#' @importFrom magrittr set_colnames
#' @importFrom tidyr gather
#' @importFrom dplyr group_by mutate select summarise summarise_all
#' @param x A vector of values to evaluate performance over
#' @param par_df A data frame like the one produced by perform_df(). MUST contain columns named: draw, species, x_min, x_max, shape1, shape2, stretch, and nu.
#' @return A tidy data frame, containing a points along the environmental axis, and corresponding points for the performance axis, for each group and posterior draw input.
#' @export
#'
posterior_predict <- function(x, par_df){
  if(missing(x)){
    x <- seq(min(par_df$x_min), max(par_df$x_max), length.out = 100)
  }

  1:nrow(par_df) %>% map_df(~{
    draw_x <- par_df$draw[.x]
    species <- par_df$species[.x]
    x_min <- par_df$x_min[.x]
    x_max <- par_df$x_max[.x]
    shape1 <- par_df$shape1[.x]
    shape2 <- par_df$shape2[.x]
    stretch <- par_df$stretch[.x]
    nu <- par_df$nu[.x]
    mu <- performance_mu(x, shape1, shape2, stretch, x_min, x_max)
    zero_idx <- x < x_min | x > x_max
    mu_species <- mu %>%
      map_dbl(function(x){
        rnorm(n = 1, mean = x, sd = (1+x)^2*1/nu) %>%
          (function(z) ifelse(z < 0, 0, z))
      }) %>%
      replace(zero_idx, 0)

    tibble(x = x,
           trait = mu_species,
           mu = mu,
           species = rep(species, length(mu)),
           draw = rep(draw_x, length(mu))
    )
  })
}


#' Generate psuedo-observed quantiles from performance curve parameters
#'
#' Constructs a tidy data frame of generated data given a set of input parameters.
#' @importFrom purrr map map_lgl map_dbl map_df
#' @importFrom magrittr set_colnames
#' @importFrom tidyr gather
#' @importFrom dplyr group_by mutate select summarise summarise_all
#' @importFrom stringr str_glue
#' @param x A vector of values to evaluate performance over
#' @param par_df A data frame like the one produced by perform_df(). MUST contain columns named: draw, species, x_min, x_max, shape1, shape2, stretch, and nu.
#' @param p Vector of probability values passed to qnorm -- the amount of probability density right of the returned quantiles
#' @return A tidy data frame, containing upper and lower quantiles along the environmental axis, and corresponding points for the performance axis, for each group and posterior draw input.
#' @export
#'
posterior_quantile <- function(x, p, par_df){
  if(missing(x)){
    x <- seq(min(par_df$x_min), max(par_df$x_max), length.out = 100)
  }

  if(min(p) < 0.5){
    warning(str_glue("p < 0.5 are not sensible here, your lowest is {min(p)}."))
  }

  low_q <- (1 - p)/2
  hi_q <- 1 - low_q

  1:nrow(par_df) %>% map_df(~{
    draw_x <- par_df$draw[.x]
    species <- par_df$species[.x]
    x_min <- par_df$x_min[.x]
    x_max <- par_df$x_max[.x]
    shape1 <- par_df$shape1[.x]
    shape2 <- par_df$shape2[.x]
    stretch <- par_df$stretch[.x]
    nu <- par_df$nu[.x]
    mu <- performance_mu(x, shape1, shape2, stretch, x_min, x_max)
    lower <- mu %>%
      map_dbl(function(x){
        qnorm(p = low_q, mean = x, sd = (1+x)^2*1/nu) %>%
          (function(z) ifelse(z < 0, 0, z))
      })
    upper <- mu %>%
      map_dbl(function(x){
        qnorm(p = hi_q, mean = x, sd = (1+x)^2*1/nu) %>%
          (function(z) ifelse(z < 0, 0, z))
      })
    tibble(x = x,
           lower = lower,
           upper = upper,
           mu = mu,
           species = rep(species, length(mu)),
           draw = rep(draw_x, length(mu))
    )
  })
}



#' Prediction quantiles
#'
#' Generate performance prediction quantiles over multiple posterior draws. Good for visualzation.
#' @importFrom purrr map map_lgl map_dbl map_df
#' @importFrom magrittr set_colnames
#' @importFrom tidyr gather
#' @importFrom dplyr group_by mutate select summarise summarise_all
#' @param spp The species to produce predictions over, one at a time is recommended.
#' @param par_df A data frame like the one produced by perform_df(). MUST contain columns named: draw, species, x_min, x_max, shape1, shape2, stretch, and nu.
#' c p Vector of probability values passed to qnorm -- the amount of probability density right of the returned quantiles
#' @param x A vector of values to evaluate performance over
#' @return A tidy data frame, containing averaged upper and lower quantiles along the environmental axis, and corresponding points for the performance axis.
#' @export
#'

predictions <- function(x, spp, par_df, x_draws, p){

  if(missing(x)){
    x <- seq(min(par_df$x_min), max(par_df$x_max), length.out = 100)
  }

  sub_df <- par_df %>%
    filter(species == spp, draw %in% x_draws)

  preds <- p %>%
    map(~{
      posterior_quantile(x = x, par_df = sub_df, p = .x) %>%
        group_by(species, x) %>%
        summarise_all(.funs = mean) %>%
        mutate(level = .x) %>%
        arrange(x) %>%
        select(-draw)
    })
  names(preds) <- paste0("l", round(p*100, 0))
  preds
}


#' Improved method to ccalculate prediction quantiles
#'
#' Generate performance prediction quantiles over multiple posterior draws. Good for visualzation.
#' @importFrom purrr map map_lgl map_dbl map_df
#' @importFrom magrittr set_colnames
#' @importFrom tidyr gather
#' @importFrom dplyr group_by mutate select summarise summarise_all
#' @param spp The species to produce predictions over, one at a time is recommended.
#' @param par_df A data frame like the one produced by perform_df(). MUST contain columns named: draw, species, x_min, x_max, shape1, shape2, stretch, and nu.
#' c p Vector of probability values passed to qnorm -- the amount of probability density right of the returned quantiles
#' @param x A vector of values to evaluate performance over
#' @return A tidy data frame, containing averaged upper and lower quantiles along the environmental axis, and corresponding points for the performance axis.
#' @export
#'

predict_interval <- function (x, spp, par_df, x_draws, p){
  if (missing(x)) {
    x <- seq(min(par_df$x_min), max(par_df$x_max), length.out = 100)
  }
  sub_df <- par_df %>% filter(draw %in% x_draws)

  p %>% map_df(~{
    posterior_quantile(x = x, par_df = sub_df, p = .x) %>%
      group_by(species, x) %>%
      summarise_all(.funs = mean) %>%
      mutate(level = .x) %>%
      arrange(x) %>%
      select(-draw) %>%
      mutate(level = round(.x * 100, 0))
  }) %>%
    arrange(species, x)
}




#' Draw parameter values from priors for simulation
#'
#' Intended for internal use only.
#' @importFrom purrr map map_lgl map_dbl map_df
#' @importFrom magrittr set_colnames
#' @importFrom tidyr gather
#' @importFrom dplyr group_by mutate select summarise summarise_all
#' @importFrom truncnorm rtruncnorm
#' @export
#'

generate_parameters <- function(
  n_spp,
  mu_sig,
  shape1_pr_mu,
  shape2_pr_mu,
  stretch_pr_mu,
  min_pr_mu,
  max_pr_mu,
  nu_pr_shape,
  nu_pr_scale,
  pr_sig = 1
){


  mu_shape1 <- rnorm(1, shape1_pr_mu, mu_sig)
  shape1 <- truncnorm::rtruncnorm(n_spp, a = 2, mean = mu_shape1, sd = pr_sig)

  mu_shape2 <- rnorm(1, shape2_pr_mu, mu_sig)
  shape2 <- truncnorm::rtruncnorm(n_spp, a = 2, mean = mu_shape2, sd = pr_sig)

  mu_stretch <- rnorm(1, stretch_pr_mu, mu_sig)
  stretch <- truncnorm::rtruncnorm(n_spp, a = 0, mean = mu_stretch, sd = pr_sig)

  mu_min <- rnorm(1, min_pr_mu, mu_sig)
  x_min <- rnorm(n_spp, mu_min, pr_sig)

  mu_max <- rnorm(1, max_pr_mu, mu_sig)
  x_max <- rnorm(n_spp, mu_max, pr_sig)

  ed_test <- 1:n_spp %>% map_lgl(~x_min[.x] > x_max[.x]) %>% sum()
  if(ed_test > 0){
    stop("parameter d is greater than parameter e")
  }

  mu_nu <- truncnorm::rtruncnorm(1, a = 0, mean = nu_pr_scale, sd = 1)
  nu <- rgamma(n_spp, shape = nu_pr_shape, scale = mu_nu)
  #nu <- truncnorm::rtruncnorm(n_spp, a = 0, mean = mu_nu, sd = 1)


  #assume same nu for each group
  # nu <- rgamma(1, shape = 10, scale = 1) %>%
  #   rep(n_spp)
  #
  list(shape1 = shape1, shape2 = shape2, stretch = stretch, x_min = x_min, x_max = x_max, nu = nu,
       mu_shape1 = mu_shape1, mu_shape2 = mu_shape2, mu_stretch = mu_stretch,
       mu_min = mu_min, mu_max = mu_max, mu_nu = mu_nu
  )
}


#' Generates simulated data that can be used for practice and model validation
#'
#' Generates correctly formatted input data for the performance_stan() function.
#' @importFrom purrr map map_lgl map_dbl map_df
#' @importFrom magrittr set_colnames
#' @importFrom tidyr gather
#' @importFrom dplyr group_by mutate select summarise summarise_all
#' @param n_spp The number of groups (i.e. species) to be simulated
#' @param n_axis The number of sampling locations along environmental axis
#' @param n_reps The number of times each location along the environmental axis should be sampled
#' @param q_end Quantile [0,1] to determine how close to the ends of the axis sampling should start.
#' @details All other arguments are hyper prior values used to generate species-level parameters. The arguments for the two shape parameters should remain above two, stretch above zero.
#' @return A list containing the true parameters for each species (named true_params), and tidy data frame, containing columns for environmental axis, response trait, mean response (handy for plotting the curve), and groups (named sim_data).
#' @export
#' @examples
#' simulate_data()

simulate_data <- function(
  n_spp = 2,
  n_axis = 20,
  n_reps = 2,
  q_end = 0.9,
  mu_sig = 1,
  shape1_pr_mu = 2,
  shape2_pr_mu = 2,
  stretch_pr_mu = 2,
  min_pr_mu = -5,
  max_pr_mu = 5,
  nu_pr_shape = 8,
  nu_pr_scale = 3,
  pr_sig = 1
){
  true_params <- generate_parameters(
    n_spp = n_spp,
    mu_sig = mu_sig,
    shape1_pr_mu = shape1_pr_mu,
    shape2_pr_mu = shape2_pr_mu,
    stretch_pr_mu = stretch_pr_mu,
    min_pr_mu = min_pr_mu,
    max_pr_mu = max_pr_mu,
    nu_pr_shape = nu_pr_shape,
    nu_pr_scale = nu_pr_scale,
    pr_sig = pr_sig
  )

  #where along axis to sample
  q_end
  x_lo <- quantile(true_params$x_min, q_end)
  x_hi <- quantile(true_params$x_max, 1-q_end)

  #in case max is less than min
  if(x_lo > x_hi){
    x_lo <- x_lo + x_hi
    x_hi <- x_lo - x_hi
    x_lo <- x_lo - x_hi
  }

  xseq <- seq(x_lo, x_hi, length.out = n_axis) %>%
    rep(each = n_reps)

  #generated data
  sim_data <- 1:n_spp %>%
    map_df(~{
      true_df <- tibble(
        shape1 = true_params$shape1[.x],
        shape2 = true_params$shape2[.x],
        stretch = true_params$stretch[.x],
        x_min = true_params$x_min[.x],
        x_max = true_params$x_max[.x],
        nu = true_params$nu[.x],
        species = .x,
        draw = 1
      )
      posterior_predict(xseq, true_df) %>%
        select(-draw)
    })

  list(true_params = true_params, sim_data = sim_data)
}

