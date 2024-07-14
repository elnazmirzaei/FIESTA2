#' 19v Normalization Function
#'
#' The main goal of this function is normalize the scRNA-seq data  
#' @param Data matrix UnImputed gene expression data, which represent gene expressions for each cell/sample. rows = genes and cols= cells. 
#' @keywords Norm19V
#' @examples
#' Norm19V ( Data )
#' @export


q19 <- function(X){
    Xpos = X[which(X>0)]
    q_90 <- quantile(Xpos, 0.9)
	q_95 <- quantile(Xpos, 0.95)
	return(sum(Xpos[Xpos >= q_90 & Xpos <= q_95]))
}

Norm19V <- function(Data){
    q = apply(Data, 2, q19)
    A_norm <- sweep(Data, 2, q, '/');
	if(length(which(q<=10))>0){
		A_norm[,which(q<=10)] = 0
	}
    zscore <- (q - mean(q)) / sd(q)
    q1 = q[which(abs(zscore)<2)]
	if(length(which(q<=10))>0){
		q1 = q1[which(q1>10)]
	}
    c = quantile(q1,.25)
    A_norm <- A_norm * c
    return(A_norm)
}
