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
#' @param data_set universal dataset
#' @param string1 3 character string of cell type 1
#' @param string2 3 character sting of cell type 2
#'
#' @return a new dataset, that is ready for a visualisation 
#' @export
#'
#' @examples
augmenting_data <- function(data_set, string1, string2){
  data_set_for_visualisation <- data_set |> 
    group_by(protein_groups,cell_type) |> 
    filter(cell_type == string1 | cell_type == string2) |> 
    summarise(
      mean = mean(intensity),
      sd = sd(intensity)) |> 
    pivot_wider(names_from = cell_type, values_from = c(mean, sd)) |> 
    mutate(fold_log2 = log2(paste0("mean_", string2) /paste0("mean_", string1)), 
           p_val = pval(paste0("mean_", string1), paste0("mean_", string2), 4, 4, paste0("sd_", string1), paste0("sd_", string2)),
           q_val = (p.adjust(p_val))) |> 
    mutate(expression = case_when(fold_log2 > 0 & q_val <= 0.05~ "overexpressed",
                                  fold_log2 < 0 & q_val <= 0.05 ~ "underexpressed",
                                  q_val >0.05  ~ "not_significant")) |> 
    select(protein_groups, fold_log2, q_val, expression)
  return(data_set_for_visualisation)
  
}




