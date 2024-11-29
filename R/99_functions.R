#' Title
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
  #Calculating test parametrs
  tobs <- ((mu1-mu2)-0)/sqrt((s1^2)/n1+(s2^2)/n2)
  #Calculating degrees of freedom using the Welch-Satterthwaite equation (not sure???)
  df <- (((s1^2)/n1+(s2^2)/n2)^2)/((((s1^2)/n1)^2)/(n1-1)+(((s2^2)/n2)^2)/(n2-1))
  #Computing p-value
  pvalue <- 2*(1-pt(abs(tobs), df=df))
  #Function returns p-value
  return(pvalue)
}

