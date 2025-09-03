#' Fit Multiple Closed Population Models with Heterogeneity
#'
#' Fits several closed population models to capture-recapture data using \pkg{RMark}, with different structures for detection (`p`), initial class probabilities (`pi`), and abundance (`f0`).
#'
#' @param data A data frame containing individual capture records with columns `id`, `sex`, and `week`. Will be passed to `process_data()`.
#'
#' @return A tibble summarizing model results, including estimated abundance for males, females, and total population with approximate 95% CIs, model AICc, and ΔAICc.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Defines four model structures using combinations of covariates for `p`, `pi`, and `f0`.
#'   \item Calls `process_data()` to prepare capture histories and fits each model using `RMark::mark()`.
#'   \item Extracts and summarizes estimates of population size by sex and in total, computing approximate confidence intervals using the delta method.
#'   \item Returns a tidy tibble of results with model names and ΔAICc.
#' }
#'
#' @importFrom dplyr rowwise mutate ungroup
#' @importFrom tibble tibble
#' @importFrom stringr str_remove_all
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   id = c(1, 1, 2, 2, 3, 4, 4, 5),
#'   sex = c("F", "F", "M", "M", "F", "M", "M", "F"),
#'   week = c(1, 2, 1, 3, 2, 1, 2, 3)
#' )
#' fit_models(df)
#' }
#'
#' @export
fit_models <- function(data){
  # Define model components
  p.time.mixture <- list(formula = ~time + mixture, share = TRUE)
  p.time <- list(formula = ~time, share = TRUE)
  p.time.sex <- list(formula = ~time + sex, share = TRUE)
  pi.1 <- list(formula = ~1)
  pi.sex <- list(formula = ~sex)
  f0.sex <- list(formula = ~sex)

  # Fit models with different combinations of structure
  fit1 <- process_data(data, "FullHet") %>%
    RMark::mark(model = "HetClosed",
                model.parameters = list(pi = pi.1, p = p.time.mixture, f0 = f0.sex),
                output = FALSE, silent = TRUE, delete = TRUE)

  fit2 <- process_data(data, "FullHet") %>%
    RMark::mark(model = "HetClosed",
                model.parameters = list(pi = pi.sex, p = p.time.mixture, f0 = f0.sex),
                output = FALSE, silent = TRUE, delete = TRUE)

  fit3 <- process_data(data, "Closed") %>%
    RMark::mark(model = "Closed",
                model.parameters = list(p = p.time.sex, f0 = f0.sex),
                output = FALSE, silent = TRUE, delete = TRUE)

  fit4 <- process_data(data, "Closed") %>%
    RMark::mark(model = "Closed",
                model.parameters = list(p = p.time, f0 = f0.sex),
                output = FALSE, silent = TRUE, delete = TRUE)

  # Combine model fits into a tibble
  all_fit <- tibble(fit = list(fit1, fit2, fit3, fit4))

  # Extract and summarize key results
  table <- all_fit %>%
    rowwise() %>%
    mutate(
      # Extract population size estimates and bounds
      nm = round(fit[["results"]][["derived"]][["N Population Size"]][["estimate"]][1]),
      nf = round(fit[["results"]][["derived"]][["N Population Size"]][["estimate"]][2]),
      nfm = nm + nf,
      nm_l = fit[["results"]][["derived"]][["N Population Size"]][["lcl"]][1],
      nm_u = fit[["results"]][["derived"]][["N Population Size"]][["ucl"]][1],
      nf_l = fit[["results"]][["derived"]][["N Population Size"]][["lcl"]][2],
      nf_u = fit[["results"]][["derived"]][["N Population Size"]][["ucl"]][2],

      # Approximate CI for total using delta method
      tot_var = sum(fit[["results"]][["derived.vcv"]][["N Population Size"]]),
      nfm_l = nfm / exp(1.96 * sqrt(log(1 + tot_var / nfm^2))),
      nfm_u = nfm * exp(1.96 * sqrt(log(1 + tot_var / nfm^2))),

      # AICc and model name cleanup
      AICc = fit[["results"]][["AICc"]],
      model = str_remove_all(fit[["model.name"]], "~|c\\(\\)"),

      # Formatted output columns
      females = paste0(nf, " (", round(nf_l), ", ", round(nf_u), ")"),
      males   = paste0(nm, " (", round(nm_l), ", ", round(nm_u), ")"),
      all     = paste0(nfm, " (", round(nfm_l), ", ", round(nfm_u), ")")
    ) %>%
    ungroup() %>%
    mutate(dAICc = round(AICc - min(AICc), 1))  # ΔAICc relative to best model

  return(table)
}

#' @export
pretty_ci <- function(estimate, lower, upper, decimals = 0){
  paste0(round(estimate, decimals), " (", round(lower, decimals), ", ", round(upper, decimals), ")")
}
