#' Title
#' # WE MUST NOT FORGET TO PUT source("99_functions.R") IN ALL QMDs WHERE
#' # WHERE WE ARE USING FUCTIONS!!!!!
#' 
#'
#' @param mu1 mean of the first cell type
#' @param mu2 mean of the second cell type
#' @param n1 number of replicas in the first cell type
#' @param n2 number of replicas in the second cell type
#' @param s1 standard deviation of the first cell type
#' @param s2 standard deviation of the second cell type
#'
#' @return result of a T test
#' @export
#'
#' @examples
pval <- function(mu1,mu2,n1,n2,s1,s2){
  #Calculating test parameters
  tobs <- ((mu1-mu2)-0)/sqrt((s1^2)/n1+(s2^2)/n2)
  #Calculating degrees of freedom using the Welch-Satterthwaite equation (not sure???)
  df <- (((s1^2)/n1+(s2^2)/n2)^2)/((((s1^2)/n1)^2)/(n1-1)+(((s2^2)/n2)^2)/(n2-1))
  #Computing p-value
  pvalue <- 2*(1-pt(abs(tobs), df=df))
  #Function returns p-value
  return(pvalue)
}

#' This function creates 
#'
#' @param df is the dataframe we input
#' @param string1 3 character string of cell type 1
#' @param string2 3 character sting of cell type 2
#'
#' @return a new dataset, that is ready for a visualisation 
#' @export
#'
#' @examples
volcano_augment <- function(df, later_cell, earlier_cell){
  data_set_for_visualisation <- df |> 
    group_by(protein_groups,cell_type) |> 
    filter(cell_type == later_cell | cell_type == earlier_cell) |> 
    summarise(
      mean = mean(intensity),
      sd = sd(intensity)) |> 
    pivot_wider(names_from = cell_type, values_from = c(mean, sd)) |> 
    #!!sym() is used to evaluate the result as a column name.. 
    mutate(fold_log2 = log2(!!sym(paste0("mean_", later_cell)) /!!sym(paste0("mean_", earlier_cell))), 
           p_val = pval(!!sym(paste0("mean_", later_cell)), !!sym(paste0("mean_", earlier_cell)), 4, 4, !!sym(paste0("sd_", later_cell)), !!sym(paste0("sd_", earlier_cell))),
           q_val = (p.adjust(p_val))) |> 
    mutate(expression = case_when(fold_log2 > 0 & q_val <= 0.05~ "overexpressed",
                                  fold_log2 < 0 & q_val <= 0.05 ~ "underexpressed",
                                  q_val >0.05  ~ "not significant")) |> 
    select(protein_groups, fold_log2, q_val, expression)
  return(data_set_for_visualisation)
}






