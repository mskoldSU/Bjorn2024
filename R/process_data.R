#' Prepare Capture History Data for RMark
#'
#' Processes individual-level capture data into capture histories suitable for analysis with \pkg{RMark}.
#'
#' @param data A data frame containing columns `id`, `sex`, and `week`, where each row corresponds to one observation (i.e., a capture).
#' @param model A character string specifying the model type to be passed to `RMark::process.data()`.
#'
#' @return An object of class `"RMark"` as returned by `RMark::process.data()`, containing the processed capture histories and grouping structure.
#'
#' @details
#' This function performs the following:
#' \enumerate{
#'   \item Groups data by `id`, `sex`, and `week`, creating one row per capture.
#'   \item Spreads the data into wide format with one column per week.
#'   \item Collapses weekly columns into a single capture history string (`ch`).
#'   \item Converts `sex` to a factor and uses it as a grouping variable for the RMark model.
#' }
#'
#' @importFrom dplyr group_by summarise ungroup arrange mutate
#' @importFrom tidyr pivot_wider unite
#' @importFrom RMark process.data
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   id = c(1, 1, 2, 2, 3),
#'   sex = c("F", "F", "M", "M", "F"),
#'   week = c(1, 2, 1, 3, 2)
#' )
#' process_data(df, model = "CJS")
#' }
#'
#' @export
process_data <- function(data, model){
  data %>%
    group_by(id, sex, week) %>%
    summarise(n = 1, .groups = "drop") %>%  # Ensure one record per capture
    ungroup() %>%
    arrange(week) %>%
    pivot_wider(
      id_cols = c("id", "sex"),
      names_from = week, values_from = n,
      values_fill = c(n = 0)                 # Fill non-captures with 0
    ) %>%
    unite("ch", -c("id", "sex"), sep = "") %>%  # Collapse into capture history
    mutate(sex = as.factor(sex)) %>%
    RMark::process.data(model = model, groups = "sex")  # Send to RMark
}


#' @import tidyverse
#' @import sf
#'
#' @export
prepare_data <- function(sample_file, dead_file, survey_year, region){
  hunt <- readxl::read_excel(dead_file, guess_max = 10000) |>
    janitor::clean_names() |>
    filter(year(dodsdato) == survey_year,
           str_detect(bakgrunn_arsak, "Licens")) |>
    mutate(id = str_sub(individ, 1, 8),
           date = as.Date(dodsdato),
           .keep = "none") |>
    na.omit()
  survey <- readxl::read_excel(sample_file, guess_max = 10000) |>
    janitor::clean_names() |>
    filter((funnetdato >= as.Date(paste0(survey_year, "-08-21")) & funnetdato <= as.Date(paste0(survey_year, "-10-31"))),
           analysert_av %in% c("Naturhistoriska riksmuseet", "Sveriges Lantbruksuniversitet (Umeå)"),
           kjonn_individ %in% c("Hane", "Hona")) |>
    mutate(year = survey_year,
           id = str_sub(individ, 1, 8),
           sex = kjonn_individ,
           week = isoweek(funnetdato),
           date = as.Date(funnetdato),
           north = nord_utm33_sweref99_tm,
           east = ost_utm33_sweref99_tm,
           dead = id %in% hunt$id,
           .keep = "none")
  counties <- mfomaps::SWEcounties
  regions <- beartools::survey_regions
  coords <- survey |>
    summarise(mean_east = mean(east),
              mean_north = mean(north),
              sd_east = sd(east),
              sd_north = sd(north),
              n_samples = n(),
              .by = c("id", "sex")) |>
    st_as_sf(coords = c("mean_east", "mean_north"), crs = st_crs(regions), remove = FALSE) |>
    rowwise() |>
    mutate(dist_to_border = st_distance(geometry, regions |> filter(survey_region != region) |> st_union()) |> as.numeric()) |>
    as_tibble() |>
    select(-geometry)

  NB_pars <- coords |>
    filter(dist_to_border > 20000) |>
    pull(n_samples) |>
    fit_ztNB()

  regions <- coords |>
    mutate(ac_sd =  median((sd_east[dist_to_border > 20000] + sd_north[dist_to_border > 20000])/2, na.rm = TRUE), .by = "sex") |>
    rowwise() |>
    mutate(mu = ifelse(dist_to_border > 20000,
                       list(list(mu_east = mean_east, mu_north = mean_north, p = 1)),
                       list(fit_center_NB(east = mean_east, north = mean_north,
                                          k = n_samples,
                                          sigma = ac_sd,
                                          region_name = survey_region,
                                          lambda = NB_pars$lambda, theta = NB_pars$theta)
                       )
    ),
    mu_east = mu$mu_east,
    mu_north = mu$mu_north,
    p = mu$p
    ) |> st_as_sf(coords = c("mu_east", "mu_north"), crs = st_crs(counties), remove = FALSE) |>
    st_join(counties, join = st_nearest_feature) |>
    st_join(regions, join = st_nearest_feature)
  survey |> left_join(regions |> as_tibble() |> select(id, mu_north, mu_east, survey_region, county = LnNamn)) |>
    mutate(county = paste(county, "län"), in_survey_region = survey_region == region) |>
    select(id, sex, week, date, north, east, mu_north, mu_east, survey_year = year, in_survey_region, shot_during_survey = dead, county_residence = county)
}
