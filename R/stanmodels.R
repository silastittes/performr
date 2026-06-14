# Internal environment for cached CmdStan models.
.pkg.env <- new.env(parent = emptyenv())

get_perform_model <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "performr requires the cmdstanr package. ",
      "Install it with install.packages('cmdstanr') and then run cmdstanr::install_cmdstan().",
      call. = FALSE
    )
  }

  if (!is.null(.pkg.env$perform_model)) {
    return(.pkg.env$perform_model)
  }

  stan_file <- system.file("stan", "perform.stan", package = "performr")
  if (stan_file == "") {
    stop("Could not locate perform.stan in the installed performr package.", call. = FALSE)
  }

  .pkg.env$perform_model <- cmdstanr::cmdstan_model(stan_file = stan_file)
  .pkg.env$perform_model
}
