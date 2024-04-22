
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
setClass("ptm2",
         representation(threshold="numeric"
                        ,estimates="matrix"
                        ,cov.="matrix"
				,residuals. = "matrix"
				,test.lin.="list"
                        ),
         prototype = list(threshold=numeric()
                        ,estimates=matrix(0,0,0)
                        ,cov.=matrix(0,0,0)
				,residuals.=matrix(0,0,0)
				,test.lin.=list()
				)
         )
## ***************************************************************************************************
ptm2<-function(Y, TV=NULL, Xendo=NULL, Xexo=NULL, IV=NULL, trimrate=0.4, ngrid=100, h0=1.5, Iweight=FALSE, test.lin=TRUE, B=1000){
	#
	# Y = dependent variable matrix nxT
	# TV = transition variable matrix nxT (if not specified, the first lag of Y is taken)
	# Xendo = independent variables matrix n*kendo x T
	# Xexo = independent variables matrix n*kexo x T
	# IV = instrumental variables matrix nl xT
	# trimrate = trim rate when constructing the grid for estimating the threshold.
	# ngrid = number of grid points to estimate the threshold
	# h0 = parameter for Silverman s rule of thumb for kernel estimation
	# Iweight = 1st step weight matrix. If TRUE, the identity matrix is used. If FALSE, the 1st step weight matrix is constructed from the instrumental variables.
	# test.lin = perform linearity test?
	# B = number of bootstrap iteration
    #
	# Check
	if(!is.null(TV)){
		if(ncol(Y)!=ncol(TV) & nrow(Y)!=nrow(TV)){
			stop("Y and TV matrices of incompatible dimensions")
		}
	}
	if(!is.null(IV)){
		if(ncol(IV)!=ncol(Y) | nrow(IV)%%nrow(Y)!=0 ){
			stop("Variable Y and IV of incompatible dimensions")
		}
	}
	if(!is.null(Xendo)){
		X.reg <- TRUE
		if(ncol(Xendo)!=ncol(Y) | nrow(Xendo)%%nrow(Y)!=0 ){
			stop("Variable Y and Xendo of incompatible dimensions")
		}
	}
	if(!is.null(Xexo)){
		X.reg <- TRUE
		if(ncol(Xexo)!=ncol(Y) | nrow(Xexo)%%nrow(Y)!=0 ){
			stop("Variable Y and Xexo of incompatible dimensions")
		}
	}
	if(is.null(Xendo) & is.null(Xexo)){X.reg <- FALSE}
	if(trimrate<0 | trimrate>1){
		stop("The trim rate for grid search algorithm should be between 0 and 1, extremes excluded")
	}
	if(ngrid<=0){
		stop("The number of grid points should be greater than 0")
	}
	if(h0<=0){
		stop("Bandwidth parameter should be greater than 0")
	}
	if(Iweight!=TRUE & Iweight!=FALSE){
		stop("Iweight should be either TRUE or FALSE")
	}
	# Set the data matrices: y, x, transition variable (qmat), and IV (if available)
	if(is.null(Xendo) & !is.null(Xexo)){X <- Xexo}
	if(!is.null(Xendo) & is.null(Xexo)){X <- Xendo}
	if(!is.null(Xendo) & !is.null(Xexo)){X <- rbind(Xexo,Xendo)}
	if(is.null(Xendo) & is.null(Xexo)){X <- NULL}
	if(!is.null(TV)){
		ymat <- Y
		qmat <- TV
		if(!is.null(X)){xmat <- X}
		if(!is.null(IV)){IVmat <- IV}
		if(is.null(IV)){IVmat <- c()}
	}
	if(is.null(TV)){
		ymat <- Y[,2:ncol(Y)]
		qmat <- Y[,1:(ncol(Y)-1)]
		if(!is.null(X)){xmat <- X[,2:ncol(X)]}
		if(!is.null(IV)){IVmat <- IV[,2:ncol(IV)]}
		if(is.null(IV)){IVmat <- c()}
	}
	# Prepare the regressor matrix
	#    rows X_{i=1,k=1},X_{i=2,k=1}, ... IV_{i=n,k=1}, IV_{i=1,k=2},IV_{i=2,k=2}, ... IV_{i=n,k=2}, ...
	#    columns t=1,...,T
	const <- matrix(1,nrow(ymat),1)
	# Include constant and lag in the X
	xmat_in <- cbind(const, ymat[,1:(ncol(ymat)-1)])
	# Add X variables in input (if any)
	if(!is.null(X)){xmat <- rbind(xmat_in, xmat)}
	if(is.null(X)){xmat <- xmat_in}
	# Compute estimation parameters (data dimension)
	dimN<-nrow(ymat)
	dimt<-ncol(ymat)
	k1<- (nrow(xmat)/dimN) # Number of q and x
	# Compute number of instruments
	if(is.matrix(IVmat)){
		kinst <- ncol(IVmat)/dimt
	}else{ kinst<-0 }
	# Number of endogenous X
	if(!is.null(Xendo)){kendo <- nrow(Xendo)/dimN
	}else{kendo <- 0}
	# Number of exogenous X
	if(!is.null(Xexo)){kexo <- nrow(Xexo)/dimN
	}else{kexo <- 0}
	# Number of X
	kendon <- k1-1
	# CONSTRUCTION OF MOMENT CONDITIONS MATRIX
	MomentMat <- MomentConditions.x(ymat, k1,kendon)
	mtc <- ncol(MomentMat)
	mtr <- nrow(MomentMat)
	MomentMat<- MomentMat + cbind(matrix(0, nrow=mtr, ncol=(mtc - kexo)), matrix(2,nrow=mtr, ncol=kexo))
	# CONSTRUCT GRID
	GridTh<-GridCon(qmat, trimrate, ngrid)
	#FIRST DIFF
	ymatfd <- ymat[,2:ncol(ymat)] - ymat[,1:(ncol(ymat)-1)]
	xmatfd <- xmat[,2:ncol(xmat)] - xmat[,1:(ncol(xmat)-1)]
	# CONSTRUCTION OF Z MATRIX
	Zmat<- Zmatcon(ymat, xmat, xmatfd, MomentMat, kendon, X.reg)
	# CONSTRUCTION OF INPUT MATRICES FOR GMM 
	gGMM<-GMMvar(ymatfd, xmatfd, xmat, qmat, MomentMat, X.reg)	
	ymatrep<-gGMM$ymatrep
	qtemp1<-gGMM$qtemp1
	qtemp2<-gGMM$qtemp2
	xtemp1<-gGMM$xtemp1
	xtemp2<-gGMM$xtemp2
	xtemp3<-gGMM$xtemp3
	# 1ST-STEP WEIGHT MATRIX 
	if(Iweight==FALSE){Wn <- WCon(Zmat, MomentMat, X.reg)}
	if(Iweight==TRUE){Wn <- matrix(0, nrow(Zmat), nrow(Zmat))
		diag(Wn) <- 1 }
	g1nbar <- as.matrix(rowSums(Zmat*ymatrep)/ dimN )
	# 1ST STEP GMM, ESTIMATION OVER THE GRID
	Jresult<-c()
	# A. Model with constant + X
	if(X.reg==TRUE){
		# DROP THE y_{t-1} REGRESSOR FROM THE COMPUTATION (SEE NOTE ABOVE)
		xtemp1c <- xtemp1[,((dimN+1):ncol(xtemp1))]
		xtemp2c <- xtemp2[,((dimN+1):ncol(xtemp2))]
		xtemp3c <- xtemp3[,((dimN+1):ncol(xtemp3))]
		# GRID SEARCH ALGORITHM TO FIND THE 1ST STEP THRESHOLD
		for(rr in 1:ngrid){
			rth <- GridTh[rr]	
			GMMres<-GMMcal(g1nbar,rth,Zmat, xtemp1c, xtemp2c, xtemp3c, qtemp1, qtemp2, Wn)
			Jresult[rr]<- t(GMMres$gnbar)%*% Wn %*% GMMres$gnbar
		}
	}
	# B. Model with constants only
	if(X.reg==FALSE){
		# GRID SEARCH ALGORITHM TO FIND THE 1ST STEP THRESHOLD
		for(rr in 1:ngrid){
			rth <- GridTh[rr]	
			GMMres<-GMMcal.res(g1nbar,rth,Zmat, qtemp1, qtemp2, Wn)
			Jresult[rr]<- t(GMMres$gnbar)%*% Wn %*% GMMres$gnbar
		}
	}
	# MINIMIZE J (IF MORE THAN ONE MINIMUM, KEEP THE FIRST)
	indmin<-which.min(Jresult)
	indmin<-indmin[1]
	## 1ST-STEP THRESHOLD ESTIMATE
	rhat <- GridTh[indmin]
	## 1ST-STEP GMM ESTIMATOR WITH SELECTED r_hat
	if(X.reg==TRUE){GMMres1stth<-GMMcal(g1nbar,rhat,Zmat, xtemp1c, xtemp2c, xtemp3c, qtemp1, qtemp2, Wn)}
	if(X.reg==FALSE){GMMres1stth<-GMMcal.res(g1nbar,rhat,Zmat, qtemp1, qtemp2, Wn)}
	GMMest<-GMMres1stth$GMMest
	# RESIDUALS 
	if(X.reg==TRUE){epcollect <- epcal(rhat, xtemp1c, xtemp2c, xtemp3c, qtemp1, qtemp2, ymatrep, GMMest)}
	if(X.reg==FALSE){
		indtemp2 <- 1*(qtemp1 > rhat)
		indtemp3 <- 1*(qtemp2 > rhat)
		epcollect<- ymatrep - GMMres1stth$GMMest[1,1]*(indtemp2 - indtemp3)
	}
	# 2ND STEP WEIGTH MATRIX 
	Wn2<-W2Con(epcollect, Zmat)
	# 2ND STEP GMM, ESTIMATION OVER THE GRID
	Jresult2 <- matrix(0,nrow=ngrid, ncol=1)
	if(X.reg==TRUE){
		for(rr2 in 1:ngrid){
			rth <- GridTh[rr2]	
			GMMres<-GMMcal(g1nbar,rth,Zmat, xtemp1c, xtemp2c, xtemp3c, qtemp1, qtemp2, Wn2)
			Jresult2[rr2]<- t(GMMres$gnbar)%*% Wn2 %*% GMMres$gnbar
		}
	}
	if(X.reg==FALSE){
		for(rr2 in 1:ngrid){
			rth <- GridTh[rr2]	
			GMMres<-GMMcal.res(g1nbar,rth,Zmat, qtemp1, qtemp2, Wn2)
			Jresult2[rr2]<- t(GMMres$gnbar)%*% Wn2 %*% GMMres$gnbar
		}
	}
	# MINIMIZE J (IF MORE THAN ONE MINIMUM, KEEP THE FIRST)
	indmin2<-which.min(Jresult2)
	indmin2<-indmin2[1]
	## 2ND-STEP THRESHOLD ESTIMATE
	rhat2 <- GridTh[indmin2]
	if(X.reg==TRUE){
	## 2ND-STEP GMM ESTIMATOR WITH r_hat
		GMMres2stth<-GMMcal(g1nbar,rhat2,Zmat, xtemp1c, xtemp2c, xtemp3c, qtemp1, qtemp2, Wn2)
		GMMest2<-GMMres2stth$GMMest
	}
	if(X.reg==FALSE){
	## 2ND-STEP GMM ESTIMATOR WITH r_hat
				GMMres2stth<-GMMcal.res(g1nbar,rhat2,Zmat, qtemp1, qtemp2, Wn2)
		GMMest2<-GMMres2stth$GMMest
	}
	# 2ND-STEP RESIDUALS
	if(X.reg==TRUE){	epcollect2 <- epcal(rhat2, xtemp1c, xtemp2c, xtemp3c, qtemp1, qtemp2, ymatrep, GMMest2)}
	if(X.reg==FALSE){
		epcollect2<- ymatrep - GMMres2stth$GMMest[1,1]*(indtemp2 - indtemp3)
	}
	omegahat <- W2Con(epcollect2, Zmat)
	# COVARIANCE MATRIX 
	if(X.reg==TRUE){
		covmat <- covmatfun(GMMest2, rhat2, qmat, Zmat, xtemp1c, xtemp2c,
		 xtemp3c, qtemp1, qtemp2, h0, omegahat)
	}
	if(X.reg==FALSE){
	## COVARIANCE MATRIX 
		hband <- h0 * 1.06 * dimN^(-0.2) * sqrt(var(as.vector(t(qmat))))
		dimL<-nrow(Zmat)
		Gd <- matrix(0, nrow=dimL, ncol= 1)
		Gr <- matrix(0, nrow=dimL, ncol= 1)
	
		indtemp2 <- 1*(qtemp1 > rhat2)
		indtemp3 <- 1*(qtemp2 > rhat2)

		g2nbar <- matrix(0,nrow=dimL, ncol=1)
		g2nbar[, 1] <- rowSums(Zmat*(indtemp2 - indtemp3))/dimN
	
		Gd[, 1] <- - g2nbar[, 1]

		for (ii7 in 1:dimN){
			temprc2 = dnorm((rhat2 - qtemp2[, ii7])/hband)
			temprc4 = dnorm((rhat2 - qtemp1[, ii7])/hband) 	
			tempz = (temprc2 - temprc4) * GMMres2stth$GMMest[1,1]
			Gr = Gr + Zmat[, ii7]*tempz
		}
		Gr = Gr/(dimN*hband)
		G <- cbind(Gd, Gr)
		mattoinv<- t(G)%*%omegahat%*%G
		covmat<-pinv(mattoinv)/dimN

	}
	#BOOTSTRAP
	if(test.lin==TRUE){
		if(X.reg==TRUE){
			Zmatb <- Zmat
			ymatrepb <- ymatrep
			xtemp1b <- xtemp1
			xtemp2b <- xtemp2
			xtemp3b <- xtemp3
			qtemp1b <- qtemp1
			qtemp2b <- qtemp2
			dimL <- nrow(Zmat)
			ressupwald <- supwaldcal(ymatrepb, Zmatb, xtemp1b, xtemp2b, xtemp3b, qtemp1b, 
				 qtemp2b,  MomentMat, GridTh,  ngrid, 
				 dimN)
			supwaldori <- ressupwald$supwald
			covwaldrec <- ressupwald$covwaldrec
			waldrec <- matrix(0, nrow=B, ncol=1) # To save sup-Walds
			# Bootstrap sup-Wald stat.
			for (bb2 in 1:B) {
				supwaldb <- supwaldfb(Zmatb, xtemp1b, xtemp2b, xtemp3b, qtemp1b,
					 qtemp2b, epcollect2, Wn2, covwaldrec, GridTh, ngrid, dimN)
				waldrec[bb2] <- supwaldb
			}
			bootsp <- sum(1*(waldrec > supwaldori))/B
			resboots <- list(method="parametric bootstrap-based linearity test",B=B,p.val= bootsp)
		}
		if(X.reg==FALSE){
			Zmatb <- Zmat
			ymatrepb <- ymatrep
			qtemp1b <- qtemp1
			qtemp2b <- qtemp2
			dimL <- nrow(Zmat)
			ressupwald <- supwaldcal.con(ymatrepb, Zmatb, qtemp1b, 
				 qtemp2b, MomentMat, GridTh,  ngrid, 
				 dimN)
			supwaldori <- ressupwald$supwald
			covwaldrec <- ressupwald$covwaldrec
			waldrec <- matrix(0, nrow=B, ncol=1) # To save sup-Wald*s
			# Bootstrap sup-Wald stat.
			for (bb2 in 1:B) {
				supwaldb <- supwaldfb.con(Zmatb, qtemp1b,
					 qtemp2b, epcollect2, Wn2, covwaldrec, GridTh, ngrid, dimN)
				waldrec[bb2] <- supwaldb
			}
			bootsp <- sum(1*(waldrec > supwaldori))/B
			resboots <- list(method="parametric bootstrap-based linearity test",B=B,p.val= bootsp)
		}
	}
	p.z <- 2*(1-pnorm(abs(c(GMMres2stth$GMMest,rhat2)/sqrt(diag(covmat)))))
	params <- cbind(GMMres2stth$GMMest, p.z[1:(length(p.z)-1)])
	threshold <- c(rhat2,p.z[length(p.z)])
	names(threshold) <- c("est.coef","p.val")
	colnames(params) <- c("est.coef","p.val")
	colnames(covmat) <- c(rownames(GMMres2stth$GMMest),"threshold")
	rownames(covmat) <- c(rownames(GMMres2stth$GMMest),"threshold")
	# LIST OF RESULTS
	out <- new("ptm2")
    out@threshold <- threshold;
    out@estimates <- params
    out@cov.      <- covmat;
	out@residuals. <- epcollect2;
	if(test.lin==TRUE){out@test.lin. <- resboots };
    return(out);
}

