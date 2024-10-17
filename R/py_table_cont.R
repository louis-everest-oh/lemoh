## Generate a person-years table with continuous variables
library(LTASR)
library(tidyverse)
library(powerjoin)
library(lubridate)


py_table_cont <- function(data, exposure_data,
                          outcomes, cat_vars, cont_vars, exp_vars,
                          bin = 500, mid_year = 0.4516){
  options(dplyr.summarise.inform = FALSE)
  ## Input data checks
  if(!all(c("id","start_fu","dob","dlo","hire_date","term_date")
          %in% colnames(data))){
    stop("ERROR: missing one or more mandatory variables: id, start_fu, do, dlo, hire-date, term_date")}

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

  breaks = c(0,rev(seq(nrow(data),0,-bin)))
  pb = txtProgressBar(min = 0, max = length(breaks)-1, initial = 0, style = 3)

  for (i in 1:(length(breaks)-1)) {
    init[i] <- Sys.time()
    cohort_break <- data[breaks[i]:(breaks[i+1]-1),]

    table_long <- cohort_break %>%
      rowwise() %>%
      mutate(age_sfu = floor(as.numeric(difftime(start_fu, dob,
                                                 unit = "days"))/365.25),
             day = list(seq(start_fu, dlo, by = "days"))) %>%
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

    output <- table_long %>%
      group_by(across(all_of(cat_vars))) %>%
      summarise(pdays = n(),
                across(all_of(c(outcomes, cont_var, exp_var)), ~sum(.x)))
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

  combined_output <- table_out %>%
    group_by(across(all_of(cat_vars))) %>%
    summarise(
              pyears = sum(pdays)/365.25,
              across(all_of(outcomes), ~sum(.x)),
              across(all_of(c(cont_var, exp_var)), ~sum(.x)/sum(pdays)))
  return(combined_output)
}
