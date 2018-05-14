#' Performance curve inference with Stan
#'
#' Runs performance/tolerance curve model
#' @import tidyverse
#' @param df Data frame containing response, treatment, and groups of interest.
#' @param response the unquoted column name containing the response trait (zero and positive real numbers expected).
#' @param treatment the unquoted column name containing the treatment values (real values expected).
#' @param group_ids the unquoted column name containing the fixed groups (i.e. species, populations, etc.).
#' @param file_id Optional file prefix to write Stan samples to (stan() sample_file argument).
#' @param rseed random seed (default = 123).
#' @export

stan_performance <- function(df, response,
                             treatment,
                             group_ids,
                             file_id = NULL,
                             rseed = 123,
                             max_treedepth = 10,
                             iter = 2000,
                             thin = 1,
                             adapt_delta = 0.95,
                             chains = 4,
                             shape1_pr_mu = 4,
                             shape1_pr_sig = 1,
                             shape2_pr_mu = 4,
                             shape2_pr_sig = 1,
                             stretch_pr_mu = 0,
                             stretch_pr_sig = 1,
                             min_pr_mu,
                             min_pr_sig = 1,
                             max_pr_mu,
                             max_pr_sig = 1){


  response <- enquo(response)
  treatment <- enquo(treatment)
  group_ids <- enquo(group_ids)


  sel_pull <- function(column){
    df %>% select(!! column) %>% pull()
  }

  if(missing(min_pr_mu)){
    min_pr_mu <- sel_pull(treatment) %>% min()
  }

  if(missing(max_pr_mu)){
    max_pr_mu <- sel_pull(treatment) %>% max()
  }

  #stan_out <- c(0,1) %>% map(~{ #don't need both versions at this point
  stan_in <- list(
    N = sel_pull(response) %>% length(),
    y = sel_pull(response),
    x = sel_pull(treatment),
    numSpp = sel_pull(group_ids) %>% unique() %>% length(),
    sppint = sel_pull(group_ids),
    shape1_pr_mu = shape1_pr_mu,
    shape1_pr_sig = shape1_pr_sig,
    shape2_pr_mu = shape2_pr_mu,
    shape2_pr_sig = shape2_pr_sig,
    stretch_pr_mu = stretch_pr_mu,
    stretch_pr_sig = stretch_pr_sig,
    min_pr_mu = min_pr_mu,
    min_pr_sig = min_pr_sig,
    max_pr_mu = max_pr_mu ,
    max_pr_sig = max_pr_sig
  )

  if(is.null(file_id)){
    stan_out <- sampling(
      stanmodels$perform,
      data = stan_in,
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
      iter = iter, chains = chains, thin = thin,
      seed = rseed
    )
  } else {
    stan_out <- sampling(
      stanmodels$perform,
      data = stan_in,
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
      iter = iter, chains = chains, thin = thin,
      sample_file = paste0(.x, eval(substitute(file_id)), ".samples"),
      seed = rseed
    )
  }
  stan_out
}
