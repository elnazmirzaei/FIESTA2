#' Mix Normal thresholding Function
#'
#' The main goal of this function is to apply a thresholding step to avoid over imputation  
#' @param Raw data matrix UnImputed gene expression data, which represent gene expressions for each cell/sample. rows = genes and cols= cells. 
#' @param Imputed data matrix WNMF/sNMF gene expression data, which represent gene expressions for each cell/sample. rows = genes and cols= cells. 
#' @keywords imputation
#' @examples
#' thresholding( Raw , scaledWNMF)
#' @import mixtools fitdistrplus
#' @export

thresholding <- function(Raw,Imputed){
    
    print("Percentage of zero values in the Raw data")
    print(length(which(Raw==0))/(dim(Raw)[1]*dim(Raw)[2]))
    print("Percentage of zero values in the Imputed data")
    print(length(which(Imputed==0))/(dim(Raw)[1]*dim(Raw)[2]))
    dim(Raw)
    CellMin = rep(0 , dim(Raw)[2])


    n = length(CellMin)
    for(i in c(1:n)){

        cellFit=fitdist(ceiling(Raw[,i]),"nbinom",method="mse")
        pZ = pnbinom(0,size=list(cellFit[[1]])[[1]][1],mu=list(cellFit[[1]])[[1]][1])
        CellMin[i] = round(pZ*dim(Raw)[1])

    }


   #Vectorized for
    fitRec <- function(l = 0.1){

        if(l>=1){
            return(0)
        }
        else{
            tryCatch(fit1 <<- normalmixEM2comp(x,lambda=l,sigsqrd=1,mu=c(0,1)), error=function(e) fit1 <<- fitRec(l+0.1))
            return(fit1)
        }
    }

    MixN <- function(x) {
        print(rownames(Raw)[ii])
        ii <<- ii+1
        if(!(length(unique(x))==1)){
        tryCatch(fit1 <- normalmixEM2comp(x,lambda=.5,sigsqrd=1,mu=c(0,1)), error=function(e) fit1 <- fitRec(0.1))

            if(is.list(fit1)){ 
                Mini = 1
                if (fit1$mu[1] > fit1$mu[2]){
                    Mini = 2
                }

                round((fit1$lambda[Mini]*pnorm(0, fit1$mu[Mini], fit1$sigma))*length(x))
            } 
            else{
                0
            }
        }
        else{
            0
        }
    }

    ii <- 1
    GeneMin = apply(Raw,1,MixN)


    CthreshFind <- function(X) {
        if(X[1] != 0){
        sort(X[2:length(X)])[X[1]]
    }
    else{0}
    }


    CM = apply(rbind(CellMin,Imputed),2,CthreshFind)

    GM = apply(cbind(GeneMin,Imputed),1,CthreshFind)

    A = matrix(rep(GM,length(CM)),ncol = length(CM))
    A = as.matrix(A)
    thresholds <- function(X) {
    X[which(X[2:length(X)]>X[1])] = X[1]
    X[2:length(X)]
    }

    B = apply(rbind(CM,A),2,thresholds)
    B = as.matrix(B)
    Imputed = as.matrix(Imputed)

    Imputed[which(Imputed < B)] = 0


    print("Percentage of zero values in the Imputed data after Thresholding")
    print(length(which(Imputed==0))/(dim(Raw)[1]*dim(Raw)[2]))

    return(Imputed)

}
