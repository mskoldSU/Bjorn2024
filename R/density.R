#' Kernel Density Estimation Within County Boundaries
#'
#' Computes a spatial kernel density estimate of observed individual locations
#' (east/north coordinates) restricted to the boundary of a given county.
#' The density is scaled to represent population size estimates.
#'
#' @param east Numeric vector of easting coordinates (Sweref99 TM or similar).
#' @param north Numeric vector of northing coordinates.
#' @param county_no An integer or character vector of county codes (LnKod)
#'   used to extract the spatial boundary from \code{SWEcounties}.
#' @param scale Numeric scaling factor applied to the density surface,
#'   typically derived from population estimates (e.g., `nm` or `nf`).
#' @param adjust Numeric bandwidth adjustment passed to
#'   \code{\link[spatstat.geom]{density.ppp}}. Default is 1.
#'
#' @return A data frame with columns:
#'   - `x`: easting coordinate (numeric)
#'   - `y`: northing coordinate (numeric)
#'   - `value`: scaled density estimate at the grid cell
#'
#' @details
#' - The function converts county polygons (\code{SWEcounties}) to a
#'   \code{spatstat.geom::owin} observation window, removing internal holes
#'   (\code{nngeo::st_remove_holes}).
#' - Kernel density is estimated using Scott’s rule (\code{bw.scott}) for
#'   bandwidth selection.
#' - The density is multiplied by \code{100000000 * scale} to rescale from
#'   probability density to estimated number of individuals per square unit.
#'
#' @note Requires access to external object \code{SWEcounties}.
#'
#' @importFrom sf st_union
#' @importFrom nngeo st_remove_holes
#' @importFrom spatstat.geom as.owin ppp
#' @importFrom spatstat.explore density.ppp bw.scott
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#'
#' @examples
#' \dontrun{
#' get_density(east = runif(50, 500000, 600000),
#'             north = runif(50, 6500000, 6600000),
#'             county_no = "01",
#'             scale = 1)
#' }
get_density <- function(east, north, county_no, scale, adjust = 1){
  map <- SWEcounties[SWEcounties$LnKod %in% county_no, ]
  map_owin <- map |>
    st_union() |>
    nngeo::st_remove_holes() |>
    as_Spatial() |>
    spatstat.geom::as.owin()


  d <- density.ppp(ppp(east, north, window = map_owin), bw.scott, adjust = adjust)
  mat <- d$v |> as.data.frame()
  colnames(mat) <- d$xcol
  mutate(mat, y = d$yrow) |>
    pivot_longer(-y, names_to = "x") |>
    mutate(x = as.numeric(x),
           value = value * 100000000 * scale)
}

#' Plot Estimated Population Density by Sex and Year
#'
#' This function creates contour maps of estimated population density,
#' stratified by sex and survey year, for one or more Swedish counties.
#' Density is estimated from individual location data and scaled using
#' model-derived estimates of total population size by sex.
#'
#' @param samples A data frame of individual-level samples. Must contain:
#'   - `id`: unique identifier for individuals
#'   - `sex`: sex of the individual (e.g., `"Hane"`, `"Hona"`)
#'   - `mu_east`, `mu_north`: coordinates (east/north)
#'   - `survey_year`: year of sampling
#'   - `shot_during_survey`: logical indicating if individual was shot during survey
#'
#' @param fit A data frame with population size estimates by survey year.
#'   Must include:
#'   - `survey_year`: year corresponding to estimates
#'   - `nm`: estimated number of males
#'   - `nf`: estimated number of females
#'
#' @param county_no An integer or character vector of county codes (LnKod)
#'   defining the spatial extent of the plot.
#'
#' @return A `ggplot` object showing filled density contours for estimated
#'   population density, stratified by sex and year, with administrative
#'   boundaries and observed sample locations overlaid.
#'
#' @details
#' - Density estimation is performed by \code{get_density()}, which should
#'   return a data frame with columns `x`, `y`, and `value`.
#' - The density surfaces are scaled so that the integral over sampled
#'   individuals matches the model-estimated population size (`nm` or `nf`).
#' - White polygons mask areas outside the county boundary.
#'
#' @note
#' Requires external spatial data objects \code{SWEcounties} and
#' \code{SWEmunicipalities}, as well as the helper function
#' \code{get_density()}.
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import sf
#' @importFrom stringr str_sub
#'
#' @examples
#' \dontrun{
#' plot_density(samples_df, fit_df, county_no = "01")
#' }
#' @export
plot_density <- function(samples, fit, county_no){
  data <- samples |>
    filter(shot_during_survey == FALSE) |>
    select(id, sex, east = mu_east, north = mu_north, survey_year) |>
    distinct()

  scales_male <- fit |> rowwise() |>
    mutate(n_found = n_distinct(data |> filter(sex == "Hane") |> pull(id)),
           scale = nm / n_found,
           sex = "Hane") |>
    select(survey_year, sex, scale)

  scales_female <- fit |> rowwise() |>
    mutate(n_found = n_distinct(data |> filter(sex == "Hona") |> pull(id)),
           scale = nf / n_found,
           sex = "Hona") |>
    select(survey_year, sex, scale)

  scales <- bind_rows(scales_male, scales_female)

  map <- SWEcounties[SWEcounties$LnKod %in% county_no, ]
  region_bbox <- sf::st_make_grid(sf::st_buffer(sf::st_union(map),
                                                dist = 10000), n = c(1, 1))
  region_cover <- sf::st_difference(region_bbox, sf::st_union(map))


  density_table <- data |> nest_by(sex, survey_year) |>
    left_join(scales) |>
    mutate(d = list(get_density(east = data$east, north = data$north, county_no, scale))) |>
      unnest(d)

  density_table |> ggplot() + geom_contour_filled(aes(x = x, y = y, z = value)) +
    geom_sf(data = map, fill = NA, color = "white") +
    geom_sf(data = region_cover, fill = "white") +
    theme_bw() + labs(fill = "", title = "Uppskattad populationstäthet (individer/kvadratmil) uppdelad på kön") +
    geom_sf(data = SWEmunicipalities |> filter(str_sub(KnKod, 1, 2) %in% county_no), fill = NA, color = "white") +
    geom_point(data = data, aes(x = east, y = north), color = "red", size = .01) +
    facet_grid(survey_year~sex) +
    scale_x_continuous(expand = expansion(0, 0)) +
    scale_y_continuous(expand = expansion(0, 0)) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
}
