#' findK Function
#'
#' The main goal of this function is to find a reasonable choice for k. 
#' @param A data matrix A, which represent gene expressions for each cell/sample. rows = genes and cols= cells, first column should be gene names, first row should be cell names. 
#' @param DimSen  is the Dimension Sensitivity, which can be "low", "medium" or "high", default is "medium"
#' @keywords findK 
#' @examples
#' findK(A)
#' findK(A, DimSen="high")
#' @export
#'

findK <- function (A , DimSen = 0.97 ) {

    rsvd_out <- svd(A)
    diffs <- rsvd_out$d[1:(length(rsvd_out$d)-1)] - rsvd_out$d[2:length(rsvd_out$d)]
	
	if(DimSen == "high"){
		diffsMinusOutliers <- diffs[-which(diffs>(mean(diffs)+0.5*sd(diffs)))]
		k <- length(which(diffs>mean(diffsMinusOutliers)+0.25*sd(diffsMinusOutliers)))
	}
	
	else if(DimSen == "medium"){
		diffsMinusOutliers <- diffs[-which(diffs>(mean(diffs)+0.25*sd(diffs)))]
		k <- length(which(diffs>mean(diffsMinusOutliers)+0.125*sd(diffsMinusOutliers)))
	}
	
	else if(DimSen == "low"){
		diffsMinusOutliers <- diffs[-which(diffs>(mean(diffs)+0.125*sd(diffs)))]
		k <- length(which(diffs>mean(diffsMinusOutliers)+0.062*sd(diffsMinusOutliers)))
	}
	else if(DimSen<1 & DimSen>0)
	{
		k = 1
		while(sum((diffs[1:k]))<DimSen*sum((diffs))){
			k = k + 1
		}
	
	}
	else
	{
		stop('accepted strings for DimSen is low, medium or high')
	}
	
    return (k)
}
