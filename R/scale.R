#' Scaling Function
#'
#' The main goal of this function is to scale the data to avoid any unneccesary narrowing. 
#' @param Raw data matrix UnImputed gene expression data, which represent gene expressions for each cell/sample. rows = genes and cols= cells. 
#' @param Imputed data matrix WNMF/sNMF gene expression data, which represent gene expressions for each cell/sample. rows = genes and cols= cells. 
#' @keywords imputation
#' @examples
#' scale( Raw , WNMF)
#' @export
#'


scale <- function(Raw,Imputed){
    
    Scaled = Imputed
    for( i in c(1:dim(Imputed)[1])){
        if(length(which(Imputed[i,]!=0))>0){
            n = length(which(Raw[i,]!=0))
            IMS = sort(Imputed[i,],decreasing = TRUE)[1:n]
        Scaled[i,] = (Imputed[i,]*mean(quantile(Raw[i,which(Raw[i,]!=0)], c(.9,.95))))/mean(as.numeric(quantile(IMS, c(.9,.95))))
        }
    }
    return(Scaled)


}