#' Helpers for the `get_table_rapid()` function
#'
#' @param data data from get_table_rapid()
#' @param time_break number of years defined by user
#' @param matrix_obj data.frame of total person years. Arranged by years (columns) and age (rows)
#'
matrix_cells_vec <- function(data, time_break, matrix_obj){

  # Count the number of person-years
  cmin <- min(as.numeric(colnames(matrix_obj)))
  rmin <- min(as.numeric(rownames(matrix_obj)))

  data$s1 = floor((data$start1 - cmin)/time_break)
  data$s2 = floor((data$start2 - rmin)/time_break)
  data$l1 = ceiling((data$stop1 - cmin)/time_break)
  data$l2 = ceiling((data$stop2 - rmin)/time_break)
  data$count0start = pmin(data$s1,data$s2)
  data$count0last = pmin(ncol(matrix_obj) -data$l1,
                         nrow(matrix_obj) -data$l2)

  data$ov = ifelse(data$start1 - floor(data$start1/time_break)*time_break>
                     data$start2 - floor(data$start2/time_break)*time_break, -1, 1)
  data$mv = ifelse(data$s1 < data$s2, -1, 1)
  data$ev = ifelse(data$s2==0, -1, 0)
  data$even_last = !(data$l1 -data$s1 == data$l2 -data$s2)
  data$count0lastodd = ifelse(data$even_last,
                              pmin(ncol(matrix_obj) -data$l1 - (data$ov-1)/2,
                                   nrow(matrix_obj) -data$l2 + (data$ov+1)/2),
                              pmin(ncol(matrix_obj) -data$l1,
                                   nrow(matrix_obj) -data$l2))
  data$count0lasteven = ifelse(data$even_last,
                               pmin(ncol(matrix_obj) -data$l1,
                                    nrow(matrix_obj) -data$l2),
                               pmin(ncol(matrix_obj) -data$l1 + (data$ov+1)/2,
                                    nrow(matrix_obj) -data$l2 - (data$ov-1)/2))
  data$count0starteven <- ifelse(data$ov ==-1,
                                 pmin(data$s1+1,data$s2),
                                 pmin(data$s1,data$s2+1))

  data$oos1 = ifelse(data$mv==-1,
                     data$mv*data$s1,
                     data$mv*data$s2)
  data$oos2 = ifelse(data$mv==-1,
                     data$mv*data$s2,
                     data$mv*data$s1)
  data$eos1 = ifelse(data$mv==-1,
                     data$mv*data$s1+ data$ov,
                     data$mv*data$s2)
  data$eos2 = ifelse(data$mv==-1, data$mv*data$s2,
                     data$mv*data$s1+ (data$ov)*-1)

  obj_list <- with(data, mapply(list_helper, s1,s2,l1,l2,
                                first, odd, even, last, even_last,
                                count0start,count0starteven,count0lasteven,count0lastodd,
                                ov, ev,start1, stop1))
  odd_list <- unlist(obj_list[1,], recursive = F)
  even_list <- unlist(obj_list[2,], recursive = F)

  matrix_out <- mapply(function(a,b,c,d,e,f){matrix_helper(a,b,c,d,e,f, matrix_obj)},
                       data$eos1,data$eos2, data$oos1, data$oos2,
                       even_list, odd_list)

  matrix_out <- Reduce("+", matrix_out)
  return(matrix_out)
}

#' list_helper
#'
#' @param s1 start col
#' @param s2 start row
#' @param l1 end col
#' @param l2 end row
#' @param first first cell person years
#' @param odd odd count cell person years
#' @param even even count cell person years
#' @param last last cell person years
#' @param even_last bool, the last observation on an even count cell
#' @param count0start starting col or row
#' @param count0starteven offset if starting count cell is even
#' @param count0lasteven offset if last count cell is even
#' @param count0lastodd offset if last count cell is odd
#' @param ov offset to determine the cell count order
#' @param ev offset if starting from the first cell
#' @param start1 first year of follow-up
#' @param stop1 last year of follow-up
#'
#'
list_helper <- function(s1,s2,l1,l2, first, odd, even, last, even_last,
                        count0start,count0starteven,
                        count0lasteven,count0lastodd,
                        ov, ev, start1, stop1){
  # Assign the observed person years to the strata
  if(l1-s1==1 & l2 -s2 ==1){
    odd_list <- c(rep(0,count0start),
                  stop1 - start1 + (1/365.25),
                  rep(0,count0lastodd))
    even_list <- c(rep(0,count0starteven),
                   rep(0,count0lasteven))
  }else{
    odd_list <- c(rep(0,count0start),
                  first,
                  rep(odd,max(min(l1 - s1, l2 - s2)-1- !even_last,0)),
                  last[!even_last],
                  rep(0,max(count0lastodd,0)))
    even_list <- c(rep(0,max(count0starteven,0)) ,
                   rep(even, max(min(l1 - s1, l2 - s2)-1,0)),
                   last[even_last],
                   rep(0,max(count0lasteven,0)))}
  return(list(list(odd_list),list(even_list)))
}



#' matrix_counts_vec
#'
#' @param data data from get_table_rapid()
#' @param time_break number of years defined by user
#' @param matrix_obj data.frame of total observed events. Arranged by years (columns) and age (rows)
#'
matrix_counts_vec <- function(data, time_break, matrix_obj){

  # Count the number of observed cases
  cmin <- min(as.numeric(colnames(matrix_obj)))
  rmin <- min(as.numeric(rownames(matrix_obj)))

  data$s1 = floor((data$start1 - cmin)/time_break)
  data$s2 = floor((data$start2 - rmin)/time_break)
  data$l1 = ceiling((data$stop1 - cmin)/time_break)
  data$l2 = ceiling((data$stop2 - rmin)/time_break)
  data$count0start = pmin(data$s1,data$s2)

  data$ov = ifelse(data$start1 - floor(data$start1/time_break)*time_break>
                     data$start2 - floor(data$start2/time_break)*time_break, -1, 1)
  data$mv = ifelse(data$s1 < data$s2, -1, 1)
  data$ev = ifelse(data$s2==0, -1, 0)
  data$even_last = !(data$l1 -data$s1 == data$l2 -data$s2)
  data$count0lastodd = ifelse(data$even_last,
                              pmin(ncol(matrix_obj) -data$l1 - (data$ov-1)/2,
                                   nrow(matrix_obj) -data$l2 + (data$ov+1)/2),
                              pmin(ncol(matrix_obj) -data$l1,
                                   nrow(matrix_obj) -data$l2))
  data$count0lasteven = ifelse(data$even_last,
                               pmin(ncol(matrix_obj) -data$l1,
                                    nrow(matrix_obj) -data$l2),
                               pmin(ncol(matrix_obj) -data$l1 + (data$ov+1)/2,
                                    nrow(matrix_obj) -data$l2 - (data$ov-1)/2))
  data$count0starteven <- ifelse(data$ov ==-1,
                                 pmin(data$s1+1,data$s2),
                                 pmin(data$s1,data$s2+1))

  data$oos1 = ifelse(data$mv==-1,
                     data$mv*data$s1,
                     data$mv*data$s2)
  data$oos2 = ifelse(data$mv==-1,
                     data$mv*data$s2,
                     data$mv*data$s1)
  data$eos1 = ifelse(data$mv==-1,
                     data$mv*data$s1+ data$ov,
                     data$mv*data$s2)
  data$eos2 = ifelse(data$mv==-1, data$mv*data$s2,
                     data$mv*data$s1+ (data$ov)*-1)

  data$first=0
  data$even=0
  data$odd=0
  data$last=1

  obj_list <- with(data, mapply(list_count_helper, s1,s2,l1,l2,
                                first, odd, even, last, even_last,
                                count0start,count0starteven,count0lasteven,count0lastodd,
                                ov, ev,start1, stop1))
  odd_list <- unlist(obj_list[1,], recursive = F)
  even_list <- unlist(obj_list[2,], recursive = F)

  matrix_out <- mapply(matrix_helper,
                       data$eos1,data$eos2, data$oos1, data$oos2,
                       even_list, odd_list, list(matrix_obj))
  matrix_out <- Reduce("+", matrix_out)
  return(matrix_out)
}


#' list_count_helper
#'
#' @param s1 start col
#' @param s2 start row
#' @param l1 end col
#' @param l2 end row
#' @param first first cell events
#' @param odd odds cell counts events
#' @param even even cell counts events
#' @param last last cell events
#' @param even_last bool, is the last cell even
#' @param count0start starting col or row
#' @param count0starteven offset if starting count cell is even
#' @param count0lasteven offset if last count cell is even
#' @param count0lastodd offset if last count cell is odd
#' @param ov offset to determine the cell count order
#' @param ev offset if starting from the first cell
#' @param start1 first year of follow-up
#' @param stop1 last year of follow-up
#'
list_count_helper <- function(s1,s2,l1,l2, first, odd, even, last, even_last,
                              count0start,count0starteven,
                              count0lasteven,count0lastodd,
                              ov, ev, start1, stop1){
  ## Assign the observed cases per strata

  if(l1-s1==1 & l2 -s2 ==1){
    odd_list <- c(rep(0,count0start),
                  last,
                  rep(0,count0lastodd))
    even_list <- c(rep(0,count0starteven),
                   rep(0,count0lasteven))
  }else{
    odd_list <- c(rep(0,count0start),
                  first,
                  rep(odd,max(min(l1 - s1, l2 - s2)-1- !even_last,0)),
                  last[!even_last],
                  rep(0,max(count0lastodd,0)))
    even_list <- c(rep(0,max(count0starteven,0)) ,
                   rep(even, max(min(l1 - s1, l2 - s2)-1,0)),
                   last[even_last],
                   rep(0,count0lasteven))}
  return(list(list(odd_list),list(even_list)))
}


#' matrix_helper
#'
#' @param eos1 offset cols, odd data
#' @param eos2 offset rows, even data
#' @param oos1 offset cols, odd data
#' @param oos2 offset rows, odd data
#' @param even_list even cells of person-years or events count
#' @param odd_list odd cells of person-years or events count
#' @param matrix_obj person-years or events count data.frame. Arranged by years (columns) and age (rows)
#'
matrix_helper <- function(eos1, eos2,oos1, oos2,even_list, odd_list, matrix_obj){
  matrix_obj[row(matrix_obj) +eos2 == col(matrix_obj) + eos1] <-
    matrix_obj[row(matrix_obj) +eos2 == col(matrix_obj) + eos1] + even_list
  matrix_obj[row(matrix_obj) +oos2 == col(matrix_obj) + oos1] <-
    matrix_obj[row(matrix_obj) +oos2 == col(matrix_obj) + oos1] + odd_list
  return(list(matrix_obj))

}


#' combine_matrix
#'
#' @param data data from matrix_counts_vec() and matrix_cells_vec()
#'
combine_matrix <- function(data){
  # Convert matrix to a long-form data table
  data_dt <- cbind(
    as.data.table(expand.grid(rownames(data), colnames(data))
    ), value = as.vector(data))

  return(data_dt)
}
