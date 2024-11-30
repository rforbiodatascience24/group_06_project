#' Title
#' # WHERE WE ARE USING FUCTIONS!!!!! Creating a function that finds the p-value based on the students t.test. This 'manual' pvalue computation is used because a dataframe that consists of standard deviations and mean values.
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
#' @param df tibble containing values used to create volcano plots
#' @param later_cell name of vector containing later stage cell type specific values, compared to the earlier stage
#' @param earlier_cell name of vector containing earlier tage cell type specific values
#' @param n_later the number of replicates for the later cell type
#' @param n_earlier the number of replicates for the earlier cell type
#'
#' @return a new dataset, that is ready for a visualisation 
#' @export
#'
#' @examples
volcano_augment <- function(df, later_cell, earlier_cell, n_later, n_earlier){
  data_set_for_visualisation <- df |> 
    ungroup() |> # Ungrup the data_frame to avoid miscalculations.
    select(c(protein_groups,
             !!sym(paste0("mean_", earlier_cell)), 
             !!sym(paste0("mean_", later_cell)), 
             !!sym(paste0("sd_", later_cell)), 
             !!sym(paste0("sd_", earlier_cell)))) |> 
    #!!sym() is used to evaluate the result as a column name.. 
    mutate(fold_log2 = log2(!!sym(paste0("mean_", later_cell)) /!!sym(paste0("mean_", earlier_cell))), 
           p_val = pval(!!sym(paste0("mean_", later_cell)),
                        !!sym(paste0("mean_", earlier_cell)),
                        n_later, n_earlier,
                        !!sym(paste0("sd_", later_cell)),
                        !!sym(paste0("sd_", earlier_cell))),
           q_val = (p.adjust(p_val))) |> 
    mutate(expression = case_when(fold_log2 > 0 & q_val <= 0.05 ~ "overexpressed",
                                  fold_log2 < 0 & q_val <= 0.05 ~ "underexpressed",
                                  q_val > 0.05  ~ "not significant")) |> 
    select(protein_groups, fold_log2, q_val, expression)
  return(data_set_for_visualisation)
}

