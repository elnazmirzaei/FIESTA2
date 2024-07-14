#' Zero Model thresholding Function
#'
#' The main goal of this function is to apply a thresholding step to avoid over imputation  
#' @param Raw data matrix UnImputed gene expression data, which represent gene expressions for each cell/sample. rows = genes and cols= cells. 
#' @param Imputed data matrix WNMF/sNMF gene expression data, which represent gene expressions for each cell/sample. rows = genes and cols= cells. 
#' @keywords imputation
#' @examples
#' ZeroModelThresholding ( Raw , scaledWNMF, clusters, , zeroModel = "Known", size0 = 1, mu0 = 50, buffer = 0)
#' @import mixtools fitdistrplus
#' @export


ZeroModelThresholding <- function(Raw,imputed_scaled, clusters, zeroModel = "Unknown", size0 = -1, mu0 = -1, buffer = 0.05, skew = 1){
    
    print("Percentage of zero values in the Raw data")
    print(length(which(Raw==0))/(dim(Raw)[1]*dim(Raw)[2]))
    print("Percentage of zero values in the imputed_scaled data")
    print(length(which(imputed_scaled==0))/(dim(Raw)[1]*dim(Raw)[2]))
    dim(Raw)

    if( zeroModel == "Unknown" ){
        GeneList = rownames(imputed_scaled)
        df = data.frame(size=0,mu=0)
        count = 1
        for( j in c(1:length(GeneList))){
        if (max(imputed_scaled[j,])>0){
            G = as.numeric((imputed_scaled[j,]/max(imputed_scaled[j,]))*1000)
                
                for(i in names(table(clusters))){
                    
                    z = Raw[j,which(clusters==i)]
                    x = round(G[which(clusters==i)])
                    if(max(z)==0){
                        if(max(x)>0){
                            fit1<- fitdist(x,distr = "nbinom",method = "mle",discrete=FALSE)
                            if(gofstat(fit1)$kstest == "not rejected"){
                                df = rbind(df,fit1$estimate)
                                count = count + 1
                                rownames(df)[count] = paste(j,i,sep='-')
                            }
                        }
                        
                        
                    }
                }
            }
        }

        df = df[-1,]
        #zeroModel
        size0 = exp(mean(log(df[,1])))
        mu0 = exp(mean(log(df[,2])))
    }
    else if( zeroModel == "Known" ){
        if(size0 == -1 | mu0 == -1){
            stop('size0 and/or mu0 is not provided')
        }
    }
    print("size0")
    print(size0)
    print("mu0")
    print(mu0)
    
    SWGM05 = imputed_scaled
    for( j in c(1:dim(SWGM05)[1])){
        if (max(SWGM05[j,])>0){
            G = as.numeric((SWGM05[j,]/max(SWGM05[j,]))*1000)
            G0 = SWGM05[j,]
            for(i in names(table(clusters))){
                D = G[which(clusters==i)]
                D0 = G0[which(clusters==i)]
                Data = data.frame(Sample = D,
                    zeroSample = rnbinom(length(D),size=size0, mu=mu0))

                y = base::rank(Data$Sample,ties.method = "first")
                x = base::rank(Data$zeroSample,ties.method = "first")

                for(t in c(1:dim(Data)[1])){
                    
                    #if((Data$Sample[t] - Data$zeroSample[which(x==y[t])]) < 0){
                    if((Data$Sample[t] - skew*Data$zeroSample[which(x==y[t])]) < (max(Data$Sample)*buffer)){
                    D0[t] = 0
                    }
                }
                SWGM05[j,which(clusters==i)] = D0
            }
        }
    }
    print("Percentage of zero values in the imputed_scaled_Thresholded data")
    print(length(which(SWGM05==0))/(dim(Raw)[1]*dim(Raw)[2]))
    
    return(SWGM05)

}
