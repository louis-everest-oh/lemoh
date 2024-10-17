#' Stratify Person Table Rapidly
#'
#' `get_table_rapid` reads in a data.frame/tibble containing basic demographic information
#' for each person of the cohort and stratifies the person-time and deaths into age,
#' calendar period, race, and sex strata. See `Details` for information on how the
#' person file must be formatted.
#'
#' This method uses the Maxaluso's method B which makes use of the regular intervals of age and
#' calendar-period to speed up computation time.
#'
#' The persondf tibble must contain the variables:
#' * id,
#' * gender (character: 'M'/'F'),
#' * race (character: 'W'/'N'),
#' * dob (date),
#' * pybegin (date),
#' * dlo	(date),
#' * code (character: identifying outcomes of interest)
#'
#' @param persondf data.frame like object containing one row per person with the required demographic information
#' @param strata character vector identify additional variables (found in persondf) to stratify
#' person-time. By default, gender and race are included.
#' @param break_yr atomic numeric identifying the number of years to stratify age and calendar periods.
#'
#' @return
#' A data.frame with a row for each strata containing the number of observed
#' deaths within each of the defined codes and the number of person days (pdays).
#' @export
#'
#' @examples
#' library(LTASR)
#' library(dplyr)
#'
#'Import example person file
#'person <- person_example %>%
#'     mutate(dob = as.Date(dob, format='%m/%d/%Y'),
#'     pybegin = as.Date(pybegin, format='%m/%d/%Y'),
#'     dlo = as.Date(dlo, format='%m/%d/%Y'))
#'
#'Stratify person table
#'py_table <- get_table_rapid(person, strata = stratify = c("gender", "race"))


get_table_rapid <- function(persondf,
                            strata = c(),
                            break_yr=5){

  # Define helper values
  data = persondf
  stratify = append(c("gender", "race"), strata)
  time_break = break_yr
  ymin = min(year(data$pybegin))
  ymax = max(year(data$dlo))
  amin = as.numeric(min(difftime(data$pybegin,
                                 data$dob, unit = "days")/365.25))
  amax = as.numeric(max(difftime(data$dlo,
                                 data$dob, unit = "days")/365.25))
  # Create empty matrix
  col_names <- seq(floor(ymin/time_break)*time_break, floor(ymax/time_break)*time_break, time_break)
  row_names <- seq(floor(amin/time_break)*time_break, floor(amax/time_break)*time_break, time_break)
  matrix_obj <- matrix(0, ncol = length(col_names),
                       nrow = length(row_names))
  colnames(matrix_obj) <- col_names
  rownames(matrix_obj) <- row_names
  data$start2 <-  as.numeric(difftime(data$pybegin,
                                      data$dob, unit = "days")/365.25) + (1/365.25)
  data$stop2 <-  as.numeric(difftime(data$dlo,
                                     data$dob, unit = "days")/365.25) + (1/365.25)
  data$start1 <- decimal_date(data$pybegin)
  data$stop1 <- decimal_date(data$dlo)

  data$start1 <- ifelse(data$start1%%1!=0,data$start1,data$start1+0.000001)
  data$start2 <- ifelse(data$start2%%1!=0,data$start2,data$start2+0.000001)
  data$stop1 <- ifelse(data$stop1%%1!=0,data$stop1,data$stop1+0.000001)
  data$stop2 <- ifelse(data$stop2%%1!=0,data$stop2,data$stop2+0.000001)

  # calculate the intermediate person-times
  data$dura1 <- ceiling(data$start2/time_break)*time_break - data$start2
  data$dura2 <- ceiling(data$start1/time_break)*time_break - data$start1

  data$first <- pmin(data$dura1, data$dura2)
  data$last <- pmin(data$stop1 - (floor(data$stop1/time_break)*time_break),
                    data$stop2 - (floor(data$stop2/time_break)*time_break))  + (1/365.25)
  data$even <- abs(data$dura1 - data$dura2)
  data$odd <- time_break - abs(data$dura1 - data$dura2)

  data <- data %>%
    mutate(code = paste0(code,"_", sprintf("%02d", rev)))

  # Stratify the data by the selected variables
  data_strat <- data %>%
    group_by(across(all_of(stratify))) %>%
    group_split()

  # Stratify the data by the selected variables
  data_output <- data %>%
    filter(!is.na(rev)) %>%
    group_by(code, across(all_of(stratify))) %>%
    group_split()


  # Assign the person-years to the associated variable combination
  print("Step 1/2 count person-time")
  py_mat <- pblapply(data_strat, function(x){matrix_cells_vec(x, time_break, matrix_obj)})

  print("Step 2/2 count events")
  # Assign the case counts to the associated variable combination
  output_mat <- pblapply(data_output, function(x){matrix_counts_vec(x, time_break, matrix_obj)})

  # Get the stratification variable names
  data_names <- suppressMessages(data %>%
                                   group_by(across(all_of(stratify))) %>%
                                   summarise())

  strat_names <- suppressMessages(data %>%
                                    filter(!is.na(rev)) %>%
                                    group_by(code, across(all_of(stratify))) %>%
                                    summarise())

  ## Combine person-year tables
  py_mat_comb <-  lapply(py_mat, combine_matrix)
  py_mat_comb <- Map(cbind, py_mat_comb,
                     split(data_names, 1:nrow(data_names)))
  py_mat_comb <- rbindlist(py_mat_comb)
  py_mat_comb <- py_mat_comb[py_mat_comb$value !=0,]

  output_mat_comb <-  lapply(output_mat, combine_matrix)
  output_mat_comb <- Map(cbind, output_mat_comb,
                         split(strat_names, 1:nrow(strat_names)))
  output_mat_comb <- rbindlist(output_mat_comb)
  output_mat_comb <- output_mat_comb[!is.na(code),]
  output_mat_comb <- dcast(output_mat_comb, eval(paste0("Var1+Var2+",
                                                        paste0(stratify, collapse = "+"),
                                                        "~ code", collapse = " ") ),
                           value.var = "value", fun.aggregate = sum)
  # Clean output

  output_combined <- output_mat_comb %>%
    power_full_join(py_mat_comb, by = c("Var1","Var2",stratify))
  output_combined <- output_combined[!is.na(output_combined$value),]
  output_combined[is.na(output_combined)] <- 0
  colnames(output_combined)[1:2] <- c("ageCat","CPCat")
  output_combined$ageCat <- as.numeric(as.character(output_combined$ageCat))
  output_combined$CPCat <- as.numeric(as.character(output_combined$CPCat))
  output_combined$ageCat <- paste0("[", output_combined$ageCat,
                                   ",", output_combined$ageCat + time_break, ")")
  output_combined$CPCat <- paste0("[", output_combined$CPCat,
                                  ",", output_combined$CPCat + time_break, ")")

  colnames(output_combined)[(3+length(stratify)):(ncol(output_combined)-1)] <-
    paste0("_o", rateobj$mapping$minor[match(
      colnames(output_combined)[(3+length(stratify)):(ncol(output_combined)-1)],
      paste0(rateobj$mapping$code, "_", rateobj$mapping$rev))])
  output_combined$pdays <- round(output_combined$value *365.25, 2)
  output_combined <- output_combined[,!c("value")]

  output_combined <- output_combined %>%
    relocate(pdays, .after = stratify[length(stratify)])

  return(output_combined)
}
