#' Zero Model thresholding Function
#'
#' The main goal of this function is to apply a thresholding step to avoid over imputation  
#' @param Raw data matrix UnImputed gene expression data, which represent gene expressions for each cell/sample. rows = genes and cols= cells. 
#' @param Imputed data matrix WNMF/sNMF gene expression data, which represent gene expressions for each cell/sample. rows = genes and cols= cells. 
#' @keywords imputation
#' @examples
#' ThresholdingClustered ( Raw,imputed_scaled, clusters)
#' @import mixtools fitdistrplus
#' @export


ThresholdingClustered <- function(Raw,imputed_scaled, clusters){

    t = AutoThresh(Raw = Raw, Imputed_Scaled = imputed_scaled, clusters = clusters)
    print("Percentage of zero values in the Raw data")
    print(length(which(Raw==0))/(dim(Raw)[1]*dim(Raw)[2]))
    print("Percentage of zero values in the imputed_scaled data")
    print(length(which(imputed_scaled==0))/(dim(Raw)[1]*dim(Raw)[2]))
    dim(Raw)

    SWGM05 = imputed_scaled
    for( j in c(1:dim(SWGM05)[1])){
        print(paste("j is",j,sep = " "))
        if (max(SWGM05[j,])>0){
            G = as.numeric((SWGM05[j,]))#/max(SWGM05[j,]))*1000)
            G0 = SWGM05[j,]
            GR = Raw[j,]
            for(i in names(table(clusters))){
                print(paste("i is",i,sep = " "))
                D = G[which(clusters==i)]
                D0 = G0[which(clusters==i)]
                DR = GR[which(clusters==i)]
                x = D
                l = 0.5
                threshIndex = 0

                #if(length(which(DR>0))>0){ 
                    list0 = which(DR == 0)
                    z = D0[list0]
                    meanZ = mean(z)
                    if(is.na(meanZ)){
                        meanZ = 0
                    }
                    sdZ = sd(z)
                    if(is.na(sdZ)){
                        sdZ = 1
                    }

                    listpos = which(DR>0)
                    Pos = D0[listpos]
                    meanPos = mean(Pos)
                    if(is.na(meanPos)){
                        meanPos = 1
                    }
                    sdPos = sd(Pos)
                    if(is.na(sdPos)){
                        sdPos = 1
                    }
                    
                    #Vectorized for
                    fitRec <- function(l = 0.1,meanZ,meanPos,sdZ,sdPos,flag = 0){

                        if(l>=1){
                            if(flag == 0){
                                flag = 1
                                l = 0.1
                                meanZ = 0
                                sdZ = 1
                                tryCatch(fit1 <<- normalmixEM2comp(x,lambda=l,sigsqrd=c(sdZ,sdPos),mu=c(meanZ,meanPos)), error=function(e) fit1 <<- fitRec(l+0.1,meanZ,meanPos,sdZ,sdPos,flag))
                            }
                            else{
                                return(0)
                            }
                        }
                        else{
                            tryCatch(fit1 <<- normalmixEM2comp(x,lambda=l,sigsqrd=c(sdZ,sdPos),mu=c(meanZ,meanPos)), error=function(e) fit1 <<- fitRec(l+0.1,meanZ,meanPos,sdZ,sdPos,flag))
                            return(fit1)
                        }
                    }

                    #fit1 <<- normalmixEM2comp(x,lambda=0.5,sigsqrd=c(sdZ,sdPos),mu=c(meanZ,meanPos))
                    tryCatch(fit1 <- normalmixEM2comp(x,lambda=.5,sigsqrd=c(sdZ,sdPos),mu=c(meanZ,meanPos)), error=function(e) fit1 <- fitRec(l = 0.1,meanZ,meanPos,sdZ,sdPos,flag = 0))
               #}

                if(is.list(fit1)){ 
                    Mini = 1
                    if (fit1$mu[1] > fit1$mu[2]){
                        Mini = 2
                    }
                    if(max(fit1$mu)<t){
                        D0 = rep(0,length(D0))
                    } else if(fit1$mu[Mini]<t){# & fit1$sigma[1]<0.5){
                        threshIndex = round((fit1$lambda[Mini])*length(x))
                    }
                    #else{
                    #    threshIndex = round((fit1$lambda[Mini]*pnorm(0.5, fit1$mu[Mini], fit1$sigma[Mini]))*length(x))
                    #}
                } 
                if(threshIndex>0){
                    thresh = sort(D0)[threshIndex]
                    D0[D0<thresh] = 0
                }


                SWGM05[j,which(clusters==i)] = D0
            }
        }
    }

    #In order to be faithfull to the original Raw single cell Rna seq data we return NonZera Values in Raw
    SWGM05[which(Raw>0)] = imputed_scaled[which(Raw>0)]

    print("Percentage of zero values in the imputed_scaled_Thresholded data")
    print(length(which(SWGM05==0))/(dim(Raw)[1]*dim(Raw)[2]))
    return(SWGM05)

}


normalmixEM2comp <- function(x, lambda, mu, sigsqrd, eps=1e-8, maxit=1000, verb=FALSE) {
  arbvar <- (length(sigsqrd)==2)
  mu1 <- mu[1]; mu2 <- mu[2]
  sigsqrd1 <- sigsqrd[1]; sigsqrd2 <- sigsqrd[arbvar+1]
  mx <- mean(x)
  const <- length(x) * 0.918938533204673 # i.e., times log(2*pi)/2
  dl <- 1 + eps
  iter<-0
  ll <- rep(0, maxit+1)
  a1<-(x-mu1)^2; b1<-(lambda/sqrt(sigsqrd1))*exp(-a1/2/sigsqrd1)
  a2<-(x-mu2)^2; b2<-((1-lambda)/sqrt(sigsqrd2))*exp(-a2/2/sigsqrd2)
  l <- sum(log(b1+b2+1))

  while (dl>eps && iter<maxit) {
    iter<-iter+1
    ll[iter] <- l
    postprobs <- b1/(b1+b2)
    postprobs[which(is.na(postprobs))] = 0
    lambda<-mean(postprobs)
    mu1<-mean(postprobs*x)/lambda
    mu2<-(mx-lambda*mu1)/(1-lambda)
    if (arbvar)  {   
      sigsqrd1<-mean(postprobs*a1)/lambda
      sigsqrd2<-mean((1-postprobs)*a2)/(1-lambda)
    } else {
      sigsqrd1 <- sigsqrd2 <- mean(postprobs*a1 + (1-postprobs)*a2) 
    }
    a1<-(x-mu1)^2; b1<-(lambda/sqrt(sigsqrd1))*exp(-a1/2/sigsqrd1)
    a2<-(x-mu2)^2; b2<-((1-lambda)/sqrt(sigsqrd2))*exp(-a2/2/sigsqrd2)

    oldl<-l    
    l <- sum(log(b1+b2+1))
    dl<-l-oldl
    if (verb) {
      cat("iteration =", iter, " log-lik diff =", dl, " log-lik =", 
          l-const, "\n")
    }
  }
  cat("number of iterations=", iter, "\n")
  iter <- iter+1
  ll[iter] <- l
  postprobs <- cbind(postprobs, 1-postprobs)
  colnames(postprobs) <- c(paste("comp", ".", 1:2, sep = ""))
  out <- list(x=x, lambda = c(lambda,1-lambda), mu = c(mu1, mu2), 
       sigma = sqrt(c(sigsqrd1, sigsqrd2)[1:(1+arbvar)]), 
       loglik = l - const, posterior = postprobs, 
       all.loglik=ll[1:iter] - const, 
       restarts=0, ft="normalmixEM")
  class(out) <- "mixEM"
  out
}


AutoThresh <- function(Raw = A_norm, Imputed_Scaled, clusters){
    MixNorm = c(0,0,0,0,0)
    for( j in c(1:dim(Imputed_Scaled)[1])){
            #print(paste("j is",j,sep = " "))
            if (max(Imputed_Scaled[j,])>0){
                G = as.numeric((Imputed_Scaled[j,]))#/max(Imputed_Scaled[j,]))*1000)
                G0 = Imputed_Scaled[j,]
                GR = Raw[j,]

                for(i in names(table(clusters))){
                    #print(paste("i is",i,sep = " "))
                    D = G[which(clusters==i)]
                    D0 = G0[which(clusters==i)]
                    DR = GR[which(clusters==i)]
                    x = D
                    l = 0.5
                    threshIndex = 0

                    if(max(DR) == 0){ 

                        fitRec <- function(l = 0.1,meanZ,meanPos,sdZ,sdPos,flag = 0){

                            if(l>=1){
                                if(flag == 0){
                                    flag = 1
                                    l = 0.1
                                    meanZ = 0
                                    sdZ = 1
                                    tryCatch(fit1 <<- normalmixEM2comp(x,lambda=l,sigsqrd=c(sdZ,sdPos),mu=c(meanZ,meanPos)), error=function(e) fit1 <<- fitRec(l+0.1,meanZ,meanPos,sdZ,sdPos,flag))
                                }
                                else{
                                    return(0)
                                }
                            }
                            else{
                                tryCatch(fit1 <<- normalmixEM2comp(x,lambda=l,sigsqrd=c(sdZ,sdPos),mu=c(meanZ,meanPos)), error=function(e) fit1 <<- fitRec(l+0.1,meanZ,meanPos,sdZ,sdPos,flag))
                                return(fit1)
                            }
                        }

                        
                        tryCatch(fit1 <- normalmixEM2comp(x,lambda=.5,sigsqrd=c(1,1),mu=c(0,1)), error=function(e) fit1 <- fitRec(l = 0.1,0,1,1,1,flag = 0))
                

                        if(is.list(fit1)){ 
                            Mini = 1
                            Maxi = 2
                            if (fit1$mu[1] > fit1$mu[2]){
                                Mini = 2
                                Maxi = 1
                            }


                        MixNorm = rbind(MixNorm,c(fit1$lambda[1],fit1$mu[Mini],fit1$sigma[Mini],fit1$mu[Maxi],fit1$sigma[Maxi]))

                        }

                    }
                }
            }
    }
    MixNorm.zero = MixNorm[-1,]
    colnames(MixNorm.zero) = c("Lambda","MuZero","sdZero","MuPos","sdPos")
    summary(MixNorm.zero[,"MuZero"])
    summary(MixNorm.zero[,"MuPos"])

    #Positive
    MixNorm = c(0,0,0,0,0)
    for( j in c(1:dim(Imputed_Scaled)[1])){
            #print(paste("j is",j,sep = " "))
            if (max(Imputed_Scaled[j,])>0){
                G = as.numeric((Imputed_Scaled[j,]))#/max(Imputed_Scaled[j,]))*1000)
                G0 = Imputed_Scaled[j,]
                GR = Raw[j,]

                for(i in names(table(clusters))){
                    #print(paste("i is",i,sep = " "))
                    D = G[which(clusters==i)]
                    D0 = G0[which(clusters==i)]
                    DR = GR[which(clusters==i)]
                    x = D
                    l = 0.5
                    threshIndex = 0

                    if(length(which(DR>0)) >= 10){ 
                        list0 = which(DR == 0)
                        z = D0[list0]
                        meanZ = mean(z)
                        if(is.na(meanZ)){
                            meanZ = 0
                        }
                        sdZ = sd(z)
                        if(is.na(sdZ)){
                            sdZ = 1
                        }

                        listpos = which(DR>0)
                        Pos = D0[listpos]
                        meanPos = mean(Pos)
                        if(is.na(meanZ)){
                            meanZ = 1
                        }
                        sdPos = sd(Pos)
                        if(is.na(sdPos)){
                            sdPos = 1
                        }
                        
                        fitRec <- function(l = 0.1,meanZ,meanPos,sdZ,sdPos,flag = 0){

                            if(l>=1){
                                if(flag == 0){
                                    flag = 1
                                    l = 0.1
                                    meanZ = 0
                                    sdZ = 1
                                    tryCatch(fit1 <<- normalmixEM2comp(x,lambda=l,sigsqrd=c(sdZ,sdPos),mu=c(meanZ,meanPos)), error=function(e) fit1 <<- fitRec(l+0.1,meanZ,meanPos,sdZ,sdPos,flag))
                                }
                                else{
                                    return(0)
                                }
                            }
                            else{
                                tryCatch(fit1 <<- normalmixEM2comp(x,lambda=l,sigsqrd=c(sdZ,sdPos),mu=c(meanZ,meanPos)), error=function(e) fit1 <<- fitRec(l+0.1,meanZ,meanPos,sdZ,sdPos,flag))
                                return(fit1)
                            }
                        }
                        
                        tryCatch(fit1 <- normalmixEM2comp(x,lambda=.5,sigsqrd=c(sdZ,sdPos),mu=c(meanZ,meanPos)), error=function(e) fit1 <- fitRec(l = 0.1,meanZ,meanPos,sdZ,sdPos,flag = 0))
                
                        if(is.list(fit1)){ 
                            Mini = 1
                            Maxi = 2
                            if (fit1$mu[1] > fit1$mu[2]){
                                Mini = 2
                                Maxi = 1
                            }


                            MixNorm = rbind(MixNorm,c(fit1$lambda[1],fit1$mu[Mini],fit1$sigma[Mini],fit1$mu[Maxi],fit1$sigma[Maxi]))
                        }

                    }
                }
            }
    }
    MixNorm.Pos = MixNorm[-1,]
    colnames(MixNorm.Pos) = c("Lambda","MuZero","sdZero","MuPos","sdPos")
    summary(MixNorm.Pos[,"MuZero"])
    summary(MixNorm.Pos[,"MuPos"])

    # Example data
    groupA <- MixNorm.zero[,"MuZero"]
    groupB <- MixNorm.Pos[,"MuZero"]

    threshold2 = (mean(groupA) + mean(groupB))/2

    print("threshold FOR this data is")
    print(threshold2)

    df = data.frame(Method = c(rep("zero",dim(MixNorm.zero)[1]), rep("Pos",dim(MixNorm.Pos)[1])), muFirstCurve = c(MixNorm.zero[,"MuZero"],MixNorm.Pos[,"MuZero"]))
    jpeg("test07132024.jpg",height=1200,width=800,res = 300)
    ggplot(df,aes(x = Method, y = muFirstCurve,color = Method, fill=Method)) + geom_violin() + ggtitle("MuFirstCurve") + 
    geom_hline(yintercept = threshold2, col = "darkblue") +
    geom_text(aes(1.5,threshold2,label = round(threshold2,digits = 3)),col = "darkblue", size=3)
    dev.off()

    return(threshold2)

}



