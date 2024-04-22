### ptm
### TWO- AND THREE-WAY PANEL THRESHOLD MODEL FOR CHANGE POINT DETECTION
##
##  The authors of this software are
##  Francesca Marta Lilja Di Lascio, and
##  Selene Perazzini, Copyright (c) 2024

##  Permission to use, copy, modify, and distribute this software for any
##  purpose without fee is hereby granted, provided that this entire notice
##  is included in all copies of any software which is or includes a copy
##  or modification of this software and in all copies of the supporting
##  documentation for such software.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  A copy of the GNU General Public License is available at
##  http://www.r-project.org/Licenses/

## ***************************************************************************************************
setClass("ptm3",
         representation(threshold="matrix",
			param="matrix",
			cov.="array",
			residuals. = "list",
			test.lin. = "list"
                        ),
         prototype = list(threshold=matrix(0,0,0),
				param = matrix(0,0,0),
				cov. = array(0,c(0,0,0)),
				residuals. = list(),
				test.lin. = list()
        )
)

## ***************************************************************************************************
ptm3 <- function(data., nameI, nameT, nameJ, nameY, nameTV=NULL, nameXendo=NULL, nameXexo=NULL, 
                            nameIV=NULL, trimrate=0.4, ngrid=100, h0=1.5, Iweight=FALSE, 
                            test.lin=TRUE, B=1000){
  #data. = data.frame with variables: I, T., J, Y, and if not null, TV, X, IV
  #nameI = the name of the (numerical) variable that identifies the individuals i.
  #nameT = the name of the (numerical) time variable t.
  #nameJ = the name of the (numerical) variable that indicates the third dimension j.
  #nameY = the name of the (numerical) dependent variable.
  #nameTV = the name of the (numerical) transition variable. If not specified, the first lag of Y is taken.
  #nameXendo = the names of the (numerical) endogenous independent variables (if any).
  #nameXexo = the names of the (numerical) exogenous independent variables (if any).
  #nameIV = the names of the (numerical) instrumental variables (if any).
  #trimrate = the trim rate when constructing the grid for estimating the threshold. The default value is set to 0.4.
  #ngrid = the number of grid points to estimate the threshold. The default is set to 100.
  #h0 = Silverman s rule of thumb parameter for kernel estimation. The default is set to 1.5.
  #Iweight = the 1st step weight matrix. If TRUE, the identity matrix is used. If FALSE, the 1st step weight matrix is constructed from the instrumental variables.
  #test.lin = run linearity test? TRUE or FALSE
  #B = bootstrap replication for linearity test
    #
	#check
	if(!is.data.frame(data.)){stop("data. should be a data.frame")}
	data. <- as.data.frame(data.)
	if(ncol(data.)<(4+length(nameTV)+length(nameXexo)+length(nameXendo)+length(nameIV))){stop("data.frame structure not compatible")}	
	if(any(is.na(data.))){stop("data. should not contain missing values")}
	if(length(nameI)!=1){ stop("unique identifier nameI should be specified")}
	if(length(nameT)!=1){ stop("unique identifier nameT should be specified")}
	if(length(nameJ)!=1){ stop("unique identifier nameJ should be specified")}
	if(length(nameY)!=1){ stop("unique dependent variable nameY should be specified")}
	if(!(nameI %in% colnames(data.))){stop("unable to find the selected nameI variables in the data.")}
	if(!(nameT %in% colnames(data.))){stop("unable to find the selected nameT variables in the data.")}
	if(!(nameJ %in% colnames(data.))){stop("unable to find the selected nameJ variables in the data.")}
	if(!(nameY %in% colnames(data.))){stop("unable to find the selected nameY variable in the data.")}
	# Prepare dataset
	colnames(data.)[which(colnames(data.)==nameI)] <- "i"
	colnames(data.)[which(colnames(data.)==nameJ)] <- "j"
	colnames(data.)[which(colnames(data.)==nameT)] <- "t"
	colnames(data.)[which(colnames(data.)==nameY)] <- "y"    
	# More checks	
	if(length(data.$t)%%length(unique(data.$t))!=0){stop("A balanced panel is required")}
	if(length(data.$i)%%length(unique(data.$i))!=0 | length(data.$j)%%length(unique(data.$j))!=0){stop("Statistical units should be the same for each time and level j")}
	if(is.numeric(data.$t)==FALSE | is.numeric(data.$y)==FALSE){stop("Variables T. and Y should be numeric")}
	if(is.numeric(data.$i)==FALSE | is.numeric(data.$j)==FALSE){stop("Variables nameI and nameJ should be numeric; if they are categorical, please convert them to numeric.")}
	if(!is.null(nameXendo) & is.numeric(data.[,which(colnames(data.)==nameXendo)])==FALSE){stop("Variables nameXendo should be numeric")}
	if(!is.null(nameXexo) & is.numeric(data.[,which(colnames(data.)==nameXexo)])==FALSE){stop("Variables nameXexo should be numeric")}
	if(!is.null(nameIV) & is.numeric(data.[,which(colnames(data.)==nameIV)])==FALSE){stop("Variable(s) nameIV should be numeric")}
	if(!is.null(nameTV) & is.numeric(data.[,which(colnames(data.)==nameTV)])==FALSE){stop("Variable(s) nameTV should be numeric")}
	data. <- data.[order(data.$j,data.$i,data.$t),]	
	# Set transition variable
	if(!is.null(nameTV)){
		if(length(nameTV)!=1){ stop("Unique dependent variable nameTV should be specified")}
		if(!(nameTV %in% colnames(data.))){stop("Unable to find the selected nameTV variable in the data.")}
		colnames(data.)[which(colnames(data.)==nameY)] <- "q"
	}
	if(is.null(nameTV)){
		tv <- c(NA,data.$y[1:(nrow(data.)-1)])
		data.$q <- tv
		min(data.$t)
		data. <- data.[which(data.$t!=min(data.$t)),]
	}
	# Prepare for loop: repeat estimation over j
	js <-  sort(unique(data.$j))
	nj <- length(js)
	if(!is.null(nameXexo) & !is.null(nameXendo)){nameX <- c(nameXexo,nameXendo)}
	if(!is.null(nameXexo) & is.null(nameXendo)){nameX <- c(nameXexo)}
	if(is.null(nameXexo) & !is.null(nameXendo)){nameX <- c(nameXendo)}
	if(is.null(nameXexo) & is.null(nameXendo)){nameX <- NULL}
	# A. Estimation: model with constants only
	if(is.null(nameX)){
		# Prepare the data
		data_res <- data.[,c("i", "t", "j", "y", "q")]
		# Create vectors of results
		th <- matrix(,nj,3)
		colnames(th) <-  c("j","coef", "p.val")
		delta_c <- matrix(,nj,3)
		colnames(delta_c) <-  c("j","delta.c", "p.val")
		test.lin.vec <- matrix(,nrow=nj, ncol=2)
		colnames(test.lin.vec) <- c("j","p.val")
		covmat <- array(NA,dim = c(2,2,nj))	
		residuals. <- vector(mode='list', length=nj)	
        #
		for(j in 1:nj){
			# Set j
			data_temp <- data_res[which(data.$j==js[j]),c(1,2,4:ncol(data_res))]
			# Prepare data matrices for ptm2
			ymat <- reshape(data_temp[,1:3], idvar = "i", timevar = "t", direction = "wide")
			ymat <- as.matrix(ymat[,-1])
			qmat <- reshape(data_temp[,c(1:2,4)], idvar = "i", timevar = "t", direction = "wide")
			qmat <- as.matrix(qmat[,-1])
			colnames(ymat) <- c()
			colnames(qmat) <- c()
			# If nameIV available: create a matrix with:
			#    rows nameIV_{i=1,l=1},nameIV_{i=2,l=1}, ... nameIV_{i=n,l=1}, nameIV_{i=1,l=2},nameIV_{i=2,l=2}, ... nameIV_{i=n,l=2}, ...
			#    columns t=1,...,T.
			if(!is.null(nameIV)){
				IV_temp <- data.[which(data.$j==js[j]),c("i","t",nameIV)]
				IVmat<-c()
				for(kiv in 1:length(nameIV)){
					IV_single <- reshape(IV_temp[,c(1:2,2+kiv)], idvar = "i", timevar = "t", direction = "wide")
					IVmat <- rbind(IVmat, IV_single)
				}
				IVmat <- IVmat[,-1]
				colnames(IVmat) <- c()
			} 	
			# Estimation
			if(Iweight==FALSE){	
				if(!is.null(nameIV)){PR <- ptm2(Y=ymat,TV=qmat,Xendo=NULL, Xexo=NULL, IV=IVmat, trimrate=trimrate, ngrid=ngrid, h0=h0, Iweight=FALSE, B=B, test.lin=test.lin)}
				if(is.null(nameIV)){PR <- ptm2(Y=ymat,TV=qmat,Xendo=NULL, Xexo=NULL,IV=NULL, trimrate=trimrate, ngrid=ngrid, h0=h0, Iweight=FALSE, B=B, test.lin=test.lin)}
			}
			if(Iweight==TRUE){	
				if(!is.null(nameIV)){PR <- ptm2(Y=ymat,TV=qmat,Xendo=NULL, Xexo=NULL,IV=IVmat, trimrate=trimrate, ngrid=ngrid, h0=h0, Iweight=TRUE, B=B, test.lin=test.lin)}
				if(is.null(nameIV)){PR <- ptm2(Y=ymat,TV=qmat,Xendo=NULL, Xexo=NULL,IV=NULL, trimrate=trimrate, ngrid=ngrid, h0=h0, Iweight=TRUE, B=B, test.lin=test.lin)}
			}
			# Save estimated parameters
			th[j,] <- c(j,PR@threshold)
			delta_c[j,]<-c(j,PR@estimates[1,])
			covmat[,,j] <- PR@cov.
			test.lin.vec[j,] <- c(j,PR@test.lin.$p.val)
			residuals.[[j]] <- PR@residuals.
		}
    param <- delta_c
	# Collect results in a list
#	if(test.lin==TRUE){
#		test.lin.res <- list(method="parametric bootstrap-based linearity test by varying the 3rd way parameter j", B=B, res=test.lin.vec) 
#		estimation <- list( threshold = th, param = delta_c, cov.=covmat, residuals.=residuals., test.lin = test.lin.res)
#	}
#	if(test.lin==FALSE){estimation <- list(threshold=th, param=delta_c, cov.=covmat, residuals.=residuals.)}
	}
	# B. Estimation: model with constants + regressors
	if(!is.null(nameX)){
		# Prepare the data
		data_res <- cbind(data.[,c("i", "t", "j", "y", "q", nameX)])
		# Create vectors of results
		th<-matrix(,nj,3)
		colnames(th) <-  c("j","coef", "p.val")
		param <- matrix(,nj,1+2*(1+2*(length(nameX))))
		test.lin.vec <- matrix(,nrow=nj, ncol=2)
		colnames(test.lin.vec) <- c("j","p.val")
		covmat <- array(NA,dim = c(2*length(nameX)+2,2*length(nameX)+2,nj))
		residuals. <- vector(mode='list', length=nj)
		for(j in 1:nj){	
			# Set j, prepare the dataset
			data_temp1 <- data_res[which(data_res$j==js[j]),]
			data_temp <- data_temp1[,c(1,2,4:ncol(data_temp1))]			
			# Create data matrices for ptm2
			ymat <- reshape(data_temp[,1:3], idvar = "i", timevar = "t", direction = "wide")
			ymat <- as.matrix(ymat[,-1])
			qmat <- reshape(data_temp[,c(1:2,4)], idvar = "i", timevar = "t", direction = "wide")
			qmat <- as.matrix(qmat[,-1])
			colnames(ymat) <- c()
			colnames(qmat) <- c()
			# Create regressor matrix:
			#    rows X_{i=1,k=1},X_{i=2,k=1}, ... IV_{i=n,k=1}, IV_{i=1,k=2},IV_{i=2,k=2}, ... IV_{i=n,k=2}, ...
			#    columns t=1,...,T.
			if(!is.null(nameXexo)){
				xmatexo <- c()
				for(k in 5:(4+length(nameXexo))){
					xexo_single <- reshape(data_temp[,c(1:2,k)], idvar = "i", timevar = "t", direction = "wide")
					xmatexo <- rbind(xmatexo,xexo_single)
				}
				xmatexo <- as.matrix(xmatexo[,-1])
				colnames(xmatexo) <- c()
			} else{xmatexo <- NULL}
			if(!is.null(nameXendo) & !is.null(nameXexo)){
				xmatendo <- c()
				for(k in (4+length(nameXexo)):ncol(data_temp)){
					xendo_single <- reshape(data_temp[,c(1:2,k)], idvar = "i", timevar = "t", direction = "wide")
					xmatendo <- rbind(xmatendo,xendo_single)
				}
				xmatendo <- as.matrix(xmatendo[,-1])
				colnames(xmatendo) <- c()
			}else if(!is.null(nameXendo) & is.null(nameXexo)){
				xmatendo <- c()
				for(k in 5:ncol(data_temp)){
					xendo_single <- reshape(data_temp[,c(1:2,k)], idvar = "i", timevar = "t", direction = "wide")
					xmatendo <- rbind(xmatendo,xendo_single)
				}
				xmatendo <- as.matrix(xmatendo[,-1])
				colnames(xmatendo) <- c()			
			}else{xmatendo <- NULL}
			# If nameIV available: create a matrix with:
			#    rows nameIV_{i=1,l=1},nameIV_{i=2,l=1}, ... nameIV_{i=n,l=1}, nameIV_{i=1,l=2},nameIV_{i=2,l=2}, ... nameIV_{i=n,l=2}, ...
			#    columns t=1,...,T.
			if(!is.null(nameIV)){
				IV_temp <- data.[which(data.$j==nj[j]),c("i","t",nameIV)]
				IVmat<-c()
				for(kiv in 1:length(nameIV)){
					IV_single <- reshape(IV_temp[,c(1:2,2+kiv)], idvar = "i", timevar = "t",  direction = "wide")
					IVmat <- rbind(IVmat, IV_single)
				}
				IVmat <- IVmat[,-1]
				colnames(IVmat) <- c()
			} 	
			# Estimation
			if(Iweight==FALSE){	
				if(!is.null(nameIV)){PR <- ptm2(Y=ymat,TV=qmat,Xendo=xmatendo,Xexo=xmatexo,IV=IVmat, trimrate=trimrate, ngrid=ngrid, h0=h0, Iweight=FALSE, B=B, test.lin=test.lin)}
				if(is.null(nameIV)){PR <- ptm2(Y=ymat,TV=qmat,Xendo=xmatendo,Xexo=xmatexo, IV=NULL, trimrate=trimrate, ngrid=ngrid, h0=h0, Iweight=FALSE,  B=B, test.lin=test.lin)}
			}
			if(Iweight==TRUE){	
				if(!is.null(nameIV)){PR <- ptm2(Y=ymat,TV=qmat,Xendo=xmatendo,Xexo=xmatexo, IV=IVmat, trimrate=trimrate, ngrid=ngrid, h0=h0, Iweight=TRUE, B=B, test.lin=test.lin)}
				if(is.null(nameIV)){PR <- ptm2(Y=ymat,TV=qmat,Xendo=xmatendo,Xexo=xmatexo, IV=NULL, trimrate=trimrate, ngrid=ngrid, h0=h0, Iweight=TRUE, B=B, test.lin=test.lin)}
			}
			# Save estimated parameters
			th[j,]<-c(j,PR@threshold)
			param[j,1] <- j
			param[j,2:ncol(param)] <- as.vector(t(PR@estimates))
			colnames(param) <- c("j", paste0( rep(rownames(PR@estimates), each=2), "." ,rep(colnames(PR@estimates),ncol(param)/ncol(PR@estimates))))
			test.lin.vec[j,] <- c(j,PR@test.lin.$p.val) 
			covmat[,,j]<-PR@cov.
			residuals.[[j]] <- PR@residuals. 
		}	
    	# Collect results in a list
#    	if(test.lin==TRUE){
#    		test.lin.res <- list(method="parametric bootstrap-based linearity test by varying the 3rd way level j", B=B, res=test.lin.vec ) 
#    		estimation <- list(threshold=th, param=param, cov.=covmat, residuals.=residuals., test.lin=test.lin.res)
#    	}
#    	if(test.lin==FALSE){estimation <- list( threshold=th, param=param, cov.=covmat, residuals.=residuals.)}		
    }#added a brace
    #
    out <- new("ptm3")
    out@threshold <- th;
    out@param <- param;
    out@cov. <- covmat;
    out@residuals. <- residuals.;
    if(test.lin==TRUE){
        test.lin.res <- list(method="parametric bootstrap-based linearity test by varying the 3rd way level j", B=B, res=test.lin.vec) 
        out@test.lin. <- test.lin.res
    };
    return(out);
}
