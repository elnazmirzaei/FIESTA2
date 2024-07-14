#' Imputation Function
#'
#' The main goal of this function is to Impute zero values of your data. 
#' @param A data matrix A, which represent gene expressions for each cell/sample. rows = genes and cols= cells, first column should be gene names, first row should be cell names. 
#' @param K is the number of features or estimation of number of sub-cell types, if you have no estimation for k you can initialize it to "low" "medium" or "high" and the findK algorithm will find a good k for you, default is "medium"
#' @param W is the weight matrix, which is only needed when using WNMF
#' @param loged should be TRUE if your data is loged, default is FALSE
#' @param method could be set to "WNMF" or "sNMF". default is "sNMF"
#' @keywords imputation
#' @examples
#' impute(A)
#' impute(A, method="sNMF" )
#' @import NMF
#' @export
#'
 
impute <- function(A = 0, k = 0.97, loged = FALSE, W = 0, method = "sNMF" ){

    print("preparing data..")
	if(!loged){
		A=log(A+1,2)
	}
	print("Data is prepared")
	#if(k == "medium" | k == "high" | k == "low"){
		print("preparing k..")
		k = findK(A, DimSen = k)
	#}
	print("k is ")
	print(k)
	
	print("Imputation Algorithm you chose")
	
	if(method=="brunet"){
		print("brunet")
		print("is running..")
		
		
		nonzeroA=A
		if(length(which(rowSums(A)==0))>0)
		{
			nonzeroA=A[-which(rowSums(A)==0),]
		}
		nmA=NMF::nmf(nonzeroA, k,method="brunet",seed="nndsvd",nrun=5)
		Data=basis(nmA)%*%coef(nmA)

		C1=c(1:length(A[,1]))
		if(length(which(rowSums(A)==0))>0)
		{
			C1=C1[-which(rowSums(A)==0)]
		}

		B=A
		B[C1,]=Data
		print("Done")
		return(B)
	}
	else if(method=="lee"){
		print("lee")
		print("is running..")
		
		
		nonzeroA=A
		if(length(which(rowSums(A)==0))>0)
		{
			nonzeroA=A[-which(rowSums(A)==0),]
		}
		nmA=NMF::nmf(nonzeroA, k,method="lee",seed="nndsvd",nrun=5)
		Data=basis(nmA)%*%coef(nmA)

		C1=c(1:length(A[,1]))
		if(length(which(rowSums(A)==0))>0)
		{
			C1=C1[-which(rowSums(A)==0)]
		}

		B=A
		B[C1,]=Data
		print("Done")
		return(B)
	}
	else if(method=="sNMF"){
		print("sNMF")
		print("is running..")
		
		
		nonzeroA=A
		if(length(which(rowSums(A)==0))>0)
		{
			nonzeroA=A[-which(rowSums(A)==0),]
		}
		nmA=NMF::nmf(nonzeroA, k,method="snmf/l",seed="nndsvd",nrun=5)
		Data=basis(nmA)%*%coef(nmA)

		C1=c(1:length(A[,1]))
		if(length(which(rowSums(A)==0))>0)
		{
			C1=C1[-which(rowSums(A)==0)]
		}

		B=A
		B[C1,]=Data
		print("Done")
		return(B)
	}
	else if(method=="WNMF"){
		print("WNMF")
		print("is running..")
		
		
		nonzeroA=A
		if(length(which(rowSums(A)==0))>0)
		{
			nonzeroA=A[-which(rowSums(A)==0),]
		}
	
		#Weight matrix W
		cellnonzero = colSums(nonzeroA != 0)
		lowQCellsLowerBound = mean(cellnonzero)-sd(cellnonzero)
		lowQCellsUpperBound = mean(cellnonzero)+sd(cellnonzero)
		
		genenonzero = rowSums(nonzeroA != 0)
		lowQGenes = mean(genenonzero)-sd(genenonzero)
		midQGenes = mean(genenonzero)-(0.5*sd(genenonzero))
		
		cellSum=colSums(nonzeroA)
		lowQCellsSum = mean(cellSum)-sd(cellSum)

		geneSum=rowSums(nonzeroA)
		lowQGenesSum = mean(geneSum)-sd(geneSum)
		midQGenesSum = mean(geneSum)-(0.5*sd(geneSum))
		
		W=array(3,dim=dim(nonzeroA))
		W[which(genenonzero < midQGenes),]=2
		#W[which(geneSum < midQGenesSum),]=2
		
		W[which(genenonzero < lowQGenes),]=1
		#W[which(geneSum < midQGenesSum),]=1

		W[,which(cellnonzero<lowQCellsLowerBound)]=1
		W[,which(cellnonzero>lowQCellsUpperBound)]=1
		
		#####

		
		#nmA=nmf(A, 50,method="snmf/l",seed="nndsvd",nrun=5)
		nmA=NMF::nmf(nonzeroA, rank = k, method="ls-nmf", weight=W)

		Data=basis(nmA)%*%coef(nmA)

		C1=c(1:length(A[,1]))
		if(length(which(rowSums(A)==0))>0)
		{
			C1=C1[-which(rowSums(A)==0)]
		}

		B=A
		B[C1,]=Data
		print("Done")
		return(B)
	}

}
