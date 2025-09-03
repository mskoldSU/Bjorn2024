LpolycubNB <- function(mu, x, y, k, sigma, region, lambda, theta){
  #' Log-Likelihood for the Negative Binomial Spatial Model
  #'
  #' Computes the log-likelihood contribution of an individual under a spatial Negative Binomial model,
  #' where the integration over space is approximated using the `polyCub` package.
  #'
  #' @param mu A numeric vector of length 2 giving the center of the home range (east, north).
  #' @param x Observed east coordinate of the individual.
  #' @param y Observed north coordinate of the individual.
  #' @param k Count for the individual (number of samples or detections).
  #' @param sigma Standard deviation of the home range kernel (in meters).
  #' @param region An `owin` object specifying the spatial integration region.
  #' @param lambda Mean of the Negative Binomial distribution.
  #' @param theta Dispersion parameter of the Negative Binomial distribution.
  #'
  #' @return A numeric value representing the log-likelihood contribution.
  #'
  #' @importFrom polyCub polyCub.SV
  #' @export
  f <- function(s, sigma = 5, mu){
    x <- s[,1]
    y <- s[,2]
    exp(-(x - mu[1])^2 / (2 * sigma^2)) *
      exp(-(y - mu[2])^2 / (2 * sigma^2)) /
      (2 * pi * sigma^2)
  }
  P_mu <- polyCub::polyCub.SV(region, f, sigma = sigma, mu = mu)
  - (theta + k) * log(theta / lambda + P_mu) -
    k * ((x - mu[1])^2 + (y - mu[2])^2) / (2 * sigma^2)
}



fit_center_NB <- function(east, north, k, sigma, region_name, lambda, theta){
  #' Estimate Home Range Center from Spatial Count Using NB Likelihood
  #'
  #' Finds the most likely home range center for an individual given its spatial location and count, under a spatial Negative Binomial model.
  #'
  #' @param east Observed east coordinate.
  #' @param north Observed north coordinate.
  #' @param k Count value (e.g., number of detections).
  #' @param sigma Standard deviation of the home range kernel (meters).
  #' @param region_name Name of the region of interest (not currently used; placeholder for later extension).
  #' @param lambda Mean of the Negative Binomial distribution.
  #' @param theta Dispersion parameter of the Negative Binomial distribution.
  #'
  #' @return A list with components:
  #' \describe{
  #'   \item{mu_east}{Estimated east coordinate of the home range center.}
  #'   \item{mu_north}{Estimated north coordinate of the home range center.}
  #'   \item{p}{Proportion of the region occupied by the estimated home range.}
  #' }
  #'
  #' @importFrom sf st_bbox st_as_sfc st_crs st_intersection
  #' @importFrom spatstat.geom as.owin
  #' @importFrom polyCub polyCub.SV
  #' @export
  clipbox <- sf::st_bbox(c(
    xmin = east - 10 * sigma,
    xmax = east + 10 * sigma,
    ymax = north + 10 * sigma,
    ymin = north - 10 * sigma
  ), crs = sf::st_crs(survey_region)) |>
    sf::st_as_sfc()

  region <- sf::st_intersection(clipbox, survey_region) |>
    spatstat.geom::as.owin()

  mu_hat <- optim(
    c(east, north),
    function(mu) -LpolycubNB(mu, east, north, k, sigma, region, lambda, theta),
    lower = c(east, north) - 4 * sigma,
    upper = c(east, north) + 4 * sigma,
    method = "L-BFGS-B",
    control = list(parscale = c(1000, 1000))
  )

  f <- function(s, sigma = 5, mu){
    x <- s[,1]
    y <- s[,2]
    exp(-(x - mu[1])^2 / (2 * sigma^2)) *
      exp(-(y - mu[2])^2 / (2 * sigma^2)) /
      (2 * pi * sigma^2)
  }

  P_mu <- polyCub::polyCub.SV(region, f, sigma = sigma, mu = mu_hat$par)

  list(mu_east = mu_hat$par[1], mu_north = mu_hat$par[2], p = P_mu)
}


fit_ztNB <- function(n){
  #' Fit Zero-Truncated Negative Binomial Distribution
  #'
  #' Estimates parameters of a zero-truncated Negative Binomial distribution via numerical likelihood maximization.
  #'
  #' @param n A numeric vector of observed counts (assumed > 0).
  #'
  #' @return A list with:
  #' \describe{
  #'   \item{theta}{Estimated dispersion parameter.}
  #'   \item{lambda}{Estimated mean count.}
  #' }
  #'
  #' @importFrom stats dnbinom pnbinom nlm
  #' @export
  loglik <- function(pars, n){
    -sum(
      dnbinom(n, pars[1], mu = pars[2], log = TRUE) -
        pnbinom(0, pars[1], mu = pars[2], log.p = TRUE, lower.tail = FALSE)
    )
  }
  fit <- nlm(loglik, c(5, mean(n)), n = n)
  theta <- fit$estimate[1]
  lambda <- fit$estimate[2]
  list(theta = theta, lambda = lambda)
}

sigma_table <- function(data, trim = 20000){
  #'  Calculate three versions of sigma by sex
  #'
  #' @param trim threshold for the trimmed mean
  #'
  #' @return data frame with three versions of sigma (ac_sd) and the method used to calculate it (sd_method)
  #'
  #' @export
  bind_rows(
    summarise(data, ac_sd = median((sd_east + sd_north)/ 2, na.rm = TRUE),
              sd_method = "median", .by = "sex"),
    summarise(data |> filter(sqrt((sd_east^2 + sd_north^2) / 2) < trim),
              ac_sd = sqrt(sum((n_samples - 1) * (sd_east^2 + sd_north^2) / 2, na.rm = TRUE) / sum(n_samples - 1)),
              sd_method = "trimmed", .by = "sex"),
    summarise(data |> filter(n_samples > 1),
              ac_sd = sqrt(sum((n_samples -1) * (sd_east^2 + sd_north^2) / 2, na.rm = TRUE) / sum(n_samples - 1)),
              sd_method = "pooled", .by = "sex")
  )
}

