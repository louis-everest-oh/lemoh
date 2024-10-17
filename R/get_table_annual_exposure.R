#' Stratify Person Table With Time Varying Annual Exposure
#'
#' `get_table_annual_exposure` reads in a data.frame/tibble containing basic demographic information
#' for each person of the cohort and stratifies the person-time and deaths into age,
#' calendar period, race, and sex strata. See `Details` for information on how the
#' person file must be formatted.
#'
#'
#' The person_data tibble must contain the (case specific) variables:
#' * id,
#' * dob (date),
#' * pybegin (date),
#' * dlo	(date),
#' * hire_date	(date)
#' * term_date (date)
#' The exposure_data tibble must contain the (case specific) variables:
#' * id,
#' * year (date),
#' * exposure_variable_name (numeric),
#'
#' @param person_data data.frame like object containing one row per person with the required demographic information
#' @param exposure_data data.frame like object containing one row per person per year. See
#' Details for required variables.
#' @param outcomes character vector identifying the outcome variables (found in persondf)
#' where 1 represents a case and 0 represents no case.
#' @param cat_vars character vector identifying the categorical variables in person_data on which to stratify
#' @param cont_vars a list containing continuous variables in person_data on which to stratify,
#' modeled as exp_strata objects created by exp_strata()
#' @param exp_vars a list containing continuous variables in exposure_data on which to stratify,
#' modeled as exp_strata objects created by exp_strata()
#' @param bin a number specifying how many persons to stratify at a time. Default is 500.
#' @param mid_year a number specifying the proportion of the year to define as the middle when applying a lag. Default is 0.4516
#'
#' @return
#' A data.frame with a row for each strata containing the number of observed
#' deaths within each of the defined codes and the number of person days (pdays).
#'
#' @examples
#' library(LTASR)
#' library(tidyr)
#' library(powerjoin)
#' library(lubridate)
#'
#' outcomes <- c("LUNG", "CIRCULATORY", "CVD","RESPIRATORY", "CANCER")
#' cat_vars <- c("SEX")
#' exp1 <- exp_strata(var = 'ds',
#'   cutpt =  c(0,1e-10, 5,10,20,30,40,50,60,70,80,90,100,125,150,175,200,250,
#'     300,500,750,1000,1250,1500,2000,3000,4000,5000,6000,Inf),
#'   lag = 5) # The second cutpt value should be minimally above zero if you wish for no exposure to be in a separate group
#' exp_vars <- list(exp1)
#'
#' cont1 <- exp_strata(var = 'age',
#'   cutpt = c(-Inf, seq(20,90, 5), Inf),
#'    lag = 0)
#'cont2 <- exp_strata(var = 'year',
#'   cutpt = c(seq(1950,1995,5), Inf),
#'   lag = 0)
#'cont3 <- exp_strata(var = 'ndur',
#'   cutpt = c(0,0.5,1,5,10,20,30, Inf),
#'   lag = 0)
#'cont_vars <- list(cont1, cont2,cont3)
#'
#'get_table_annual_exposure(person_data = person_data_example,
#'                        exposure_data = exposure_data_example,
#'                        outcomes = outcomes,
#'                        cat_vars = cat_vars,
#'                        cont_vars = cont_vars,
#'                        exp_vars = exp_vars)
#' @export


get_table_annual_exposure <- function(person_data, exposure_data,
                          outcomes, cat_vars, cont_vars, exp_vars,
                          bin = 500, mid_year = 0.4516){
  ## Input data checks
  if(!all(c("id","pybegin","dob","dlo","hire_date","term_date")
          %in% colnames(person_data))){
    stop("ERROR: missing one or more mandatory variables: id, pybegin, do, dlo, hire-date, term_date")}

  exp_var <- unlist(lapply(exp_vars, function(x) x$var))
  exp_lag <- unlist(lapply(exp_vars, function(x) x$lag))
  exp_cut <- (lapply(exp_vars, function(x) x$cutpt))

  cont_var <- unlist(lapply(cont_vars, function(x) x$var))
  cont_lag <- unlist(lapply(cont_vars, function(x) x$lag))
  cont_cut <- (lapply(cont_vars, function(x) x$cutpt))

  cat_vars <- c(cat_vars, paste0(c(exp_var, cont_var),"Cat"))
  table_out <- as.data.frame(matrix(ncol = (length(cat_vars) +
                                              length(cont_var) +
                                              length(exp_var) +
                                                length(outcomes) + 1)))
  colnames(table_out) <- c(cat_vars, cont_var,exp_var, outcomes, "pdays")
  table_out <- table_out[-1,]

  breaks = c(0,rev(seq(nrow(person_data),0,-bin)))
  pb = txtProgressBar(min = 0, max = length(breaks)-1, initial = 0, style = 3)

  for (i in 1:(length(breaks)-1)) {
    init[i] <- Sys.time()
    cohort_break <- person_data[breaks[i]:(breaks[i+1]-1),]

    table_long <- cohort_break %>%
      rowwise() %>%
      mutate(age_sfu = floor(as.numeric(difftime(pybegin, dob,
                                                 unit = "days"))/365.25),
             day = list(seq(pybegin, dlo, by = "days"))) %>%
      unnest(cols = day) %>%
      mutate(year = year(day),
             age = as.numeric(difftime(day,dob, unit = "days")/365.25),
             ndur = pmin(as.numeric(difftime(day, hire_date, unit = "days")),
                         as.numeric(difftime(term_date, hire_date, unit = "days")))/365.25)


    table_long <- table_long %>%
      power_left_join(exposure_data, by = c("id", "year"), fill =0) %>%
      group_by(id) %>%
      mutate(
        across(all_of(cont_var), ~ cut(.x, breaks =
          cont_cut[cur_column()== cont_var][[1]], right = F),
          .names = "{col}Cat"),
        across(all_of(exp_var), ~lag(.x, floor((mid_year+
                                               as.numeric(exp_lag[
                                                 cur_column()==exp_var]))*365.25),
                                  default = 0)),
        across(all_of(exp_var), ~cut(.x, breaks =
          exp_cut[cur_column()== exp_var][[1]], right = F),
          .names = "{col}Cat"),
        across(all_of(outcomes), ~case_when(dlo==day ~ .x, T ~ 0)))

    output <- suppressMessages(table_long %>%
      group_by(across(all_of(cat_vars))) %>%
      summarise(pdays = n(),
                across(all_of(c(outcomes, cont_var, exp_var)), ~sum(.x))))
    end[i] <- Sys.time()
    setTxtProgressBar(pb,i)
    time <- round(seconds_to_period(sum(end - init)), 0)
    est <- (length(breaks)-1) * (mean(end[end != 0] -
                                        init[init != 0])) - time
    remainining <- round(seconds_to_period(est), 0)

    cat(paste(" // Execution time:", time,
              " // Estimated time remaining:", remainining), "")
    rm(table_long)
    table_out <- rbind(table_out, output)
  }
  close(pb)

  combined_output <- suppressMessages(table_out %>%
    group_by(across(all_of(cat_vars))) %>%
    summarise(
              pyears = sum(pdays)/365.25,
              across(all_of(outcomes), ~sum(.x)),
              across(all_of(c(cont_var, exp_var)), ~sum(.x)/sum(pdays))))
  return(combined_output)
}
