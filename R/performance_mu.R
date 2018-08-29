#' Performance curve
#'
#' Deterministic function representing an individual performance/tolerance curve.
#' @import purrr
#' @param xs Vector of environmental axis values
#' @param shape1 First parameter controlling asymmetry along environmental axis.
#' @param shape2 Second parameter controlling asymmetry along environmental axis.
#' @param stretch Parameter controlling stretch along the performance axis.
#' @param x_min Parameter controlling where perforamance falls to zero below optimum.
#' @param x_max Parameter to controlling where perforamance falls to zero above optimum.
#' @return Vector of performance values.
#' @export
#' @examples
#' curve(performance_mu(xs = x, 2, 3, 1, 1, 5), from = 1, to = 5)

performance_mu <- function(xs, shape1, shape2, stretch, x_min, x_max){
  if(shape1 < 2 | shape2 < 2){
    stop("shape1 and shape2 must greater than 2")
  }
  if(x_max <= x_min){
    stop("x_max must be greater than x_min")
  }
  x <- (xs - x_min)/(x_max - x_min)
  x %>% map_dbl(~ {
    if(.x > 0 & .x < 1){
      stretch * (shape1 * shape2 * .x ^ (shape1 - 1) ) * (1 - .x ^ shape1) ^ (shape2 - 1)
    } else{
      0
    }
  })
}
