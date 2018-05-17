#' Performance data frame
#'
#' Constructs a tidy data frame from the output of stan_performance()
#' @import tidyverse
#' @import magrittr
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
      gather("Species", x) %>%
      set_colnames(c("Species", .x)) %>%
      group_by(Species) %>%
      mutate(draw = 1:n())
  }) %>%
    do.call(cbind, .) %>%
    select(-starts_with("Species")) %>%
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
#' @import tidyverse
#' @param x A [0,1] vector of values to evaluate performance over
#' @param par_df A data frame produced by perform_df().
#' @return A tidy data frame, containing a points along the environmental axis, and corresponding points for the performance axis, for each group and posterior draw input.
#' @export

map_performance <- function(x = seq(0, 1, length.out = 100), par_df){
  #generate model fits where zeros affect parameters directly
  1:nrow(par_df) %>% map_df(~{
    draw_x <- par_df$draw[.x]
    spp <- par_df$Species[.x]
    x_min <- par_df$x_min[.x]
    x_max <- par_df$x_max[.x]
    shape1 <- par_df$shape1[.x]
    shape2 <- par_df$shape2[.x]
    stretch <- par_df$stretch[.x]
    xs <- x*(x_max - x_min) + x_min
    mod_fit <- stretch*((shape1*shape2*x^(shape1-1)) * (1-x^shape1)^(shape2-1))
    data_frame(Species = spp, x = xs, y = mod_fit, draw = draw_x)
  })
}



#' Generate perforamcne curves data from parameters
#'
#' Constructs a tidy data frame of generated data given a set of input parameters.
#' @import tidyverse
#' @param x A vector of values to evaluate performance over
#' @param par_df A data frame like the one produced by perform_df(). MUST contain columns named: draw, Species, x_min, x_max, shape1, shape2, stretch, and nu.
#' @return A tidy data frame, containing a points along the environmental axis, and corresponding points for the performance axis, for each group and posterior draw input.
#' @export
posterior_predict <- function(x, par_df){
  if(missing(x)){
    x <- seq(x_min, x_max, length.out = 100)
  }

  1:nrow(par_df) %>% map_df(~{
    draw_x <- par_df$draw[.x]
    species_id <- par_df$Species[.x]
    x_min <- par_df$x_min[.x]
    x_max <- par_df$x_max[.x]
    shape1 <- par_df$shape1[.x]
    shape2 <- par_df$shape2[.x]
    stretch <- par_df$stretch[.x]
    nu <- par_df$nu[.x]
    mu <- performance_mu(x, shape1, shape2, stretch, x_min, x_max)
    zero_idx <- x < x_min | x > x_max
    mu_spp <- mu %>%
      map_dbl(function(x){
        rnorm(n = 1, mean = x, sd = (1+x^2)*1/nu) %>%
        #rnorm(n = 1, mean = x, sd = (1+x)*1/nu^2) %>%
          (function(z) ifelse(z < 0, 0, z))
      }) %>%
      replace(zero_idx, 0)

    tibble(x = x,
           trait = mu_spp,
           mu = mu,
           spp = rep(species_id, length(mu)),
           draw = rep(draw_x, length(mu))
    )
  })
}

