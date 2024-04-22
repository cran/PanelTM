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

GridCon<-function(mat, trimrate, gridnum){
	qvec<-as.vector(mat)
	gridlower<-quantile(qvec,trimrate/2)
	gridupper<-quantile(qvec,1-trimrate/2)
	gridth<-seq(gridlower, gridupper, length.out=gridnum)
return(gridth)}

## ***************************************************************************************************

MomentConditions.x<- function(ymat, k1,kendon){
	dimt <- ncol(ymat) - 1
	dimN <- nrow(ymat)
	MomentMat <- matrix(1,nrow = dimt ,ncol = (k1 + 1))
	MomentMat[,2] <- seq(0, (dimt-1),1)
	if(kendon!=0){	
		input <- seq(0,(dimt-1),1)
		MomentMat[,(k1+2-kendon)] <- as.matrix(input ,nrow = 1 ,ncol = kendon)
	}
	return(MomentMat)
}

## ***************************************************************************************************

Zmatcon <- function(ymat, xmat, xmatfd, MomentMat, kendo, X.reg){
	dimN <- nrow(ymat)
	dimT <- ncol(ymat)
	flagz<-1
	k1<-nrow(xmat)/dimN
	if(X.reg == TRUE){dimL <- sum(MomentMat) - sum(MomentMat[, 2])}
	if(X.reg == FALSE){dimL <- sum(MomentMat)}
	Zmat <- matrix(0,nrow=dimL,ncol=dimN)
	if(X.reg==TRUE){
		for(kk1 in 1:nrow(MomentMat)){
			mtemp <- MomentMat[kk1,]
			for(ll1 in 1:ncol(MomentMat)){
				if(mtemp[ll1]==0){
					next
				} else if(ll1==1){
					Zmat[flagz,]<- matrix(1,nrow=1,ncol=dimN)
					flagz <- flagz +1
				} else if(ll1 == 2){
					next
				} else if(ll1<=(k1 + 1 - kendo)){
					Zmat[flagz,]<- t(xmatfd[((ll1-2)*dimN+1) : ((ll1-1)*dimN), kk1])
					flagz <- flagz +1
				} else {
					val <- mtemp[ll1]
					Zmat[flagz : (flagz+val-1), ] = t(xmat[(((ll1-2)*dimN+1) : ((ll1-1)*dimN)), (1:val)])
					flagz <- flagz + val
				}
			}
		}
	}
	if(X.reg==FALSE){
		for(kk1 in 1:nrow(MomentMat)){
			mtemp <- MomentMat[kk1,]
			for(ll1 in 1:ncol(MomentMat)){
				if(mtemp[ll1]==0){
					next
				} else if(ll1==1){
					Zmat[flagz,]<- matrix(1,nrow=1,ncol=dimN)
					flagz <- flagz +1
				} else if(ll1 == 2){
					Zmat[flagz:(flagz+mtemp[ll1]-1), ] <- t(ymat[, (kk1-mtemp[ll1]):(kk1-1)])
					flagz = flagz + mtemp[ll1]
				}  else {
					val <- mtemp[ll1]
					Zmat[flagz : (flagz+val-1), ] = t(xmat[(((ll1-2)*dimN+1) : ((ll1-1)*dimN)), (1:val)])
					flagz <- flagz + val
				}
			}
		}
	}
	return(Zmat)
}

## ***************************************************************************************************

GMMvar<-function(ymatfd, xmatfd, xmat, qmat, MomentMat, X.reg){
	dimt<-nrow(MomentMat)+1
	if(X.reg==TRUE){dimL <- sum(MomentMat) - sum(MomentMat[, 2])}
	if(X.reg==FALSE){dimL <- sum(MomentMat)}
	dimN<-nrow(ymatfd)
	k1<-(nrow(xmat)/(dimN))	
	ymatrep<- matrix(0, nrow = dimL , ncol = dimN)
	xtemp1 = matrix(0, nrow = dimL , ncol = (dimN*k1))
	xtemp2 = matrix(0, nrow = dimL , ncol = (dimN*k1))
	xtemp3 = matrix(0, nrow = dimL , ncol = (dimN*k1))
	qtemp1<- matrix(0, nrow = dimL , ncol = dimN)
	qtemp2<- matrix(0, nrow = dimL , ncol = dimN)
	flagr1 <- 1
	for(kk in 1:(dimt-1) ) {
		if(X.reg==TRUE){mtemp <- sum(MomentMat[kk,]) - sum(MomentMat[kk,2])}
		if(X.reg==FALSE){mtemp <- sum(MomentMat[kk,])}
		if(mtemp==0){}
		if(mtemp!=0){
			krontemp1 <- kronecker(ymatfd[, kk], matrix(1,1,mtemp))
			ymatrep[flagr1:(flagr1 + mtemp - 1),] <- t(krontemp1)
			krontemp2 <- kronecker(qmat[, kk+1], matrix(1,1,mtemp))
			qtemp1[flagr1:(flagr1 + mtemp - 1),] <- t(krontemp2)
			krontemp3 <- kronecker(qmat[, kk], matrix(1,1,mtemp))
			qtemp2[flagr1:(flagr1 + mtemp - 1),] <- t(krontemp3)
            #
			krontemp4 <- kronecker(xmatfd[,kk], matrix(1,1,mtemp))
			xtemp1[flagr1:(flagr1 + mtemp - 1),] <- t(krontemp4)
			krontemp5 <- kronecker(xmat[,kk+1], matrix(1,1,mtemp))
			xtemp2[flagr1:(flagr1 + mtemp - 1),] <- t(krontemp5)
			krontemp6 <- kronecker(xmat[,kk], matrix(1,1,mtemp))
			xtemp3[flagr1:(flagr1 + mtemp - 1),] <- t(krontemp6)
            #
		  flagr1 <- flagr1 + mtemp
		}
	}
	results <- list(ymatrep = ymatrep, qtemp1=qtemp1, qtemp2=qtemp2, xtemp1=xtemp1, xtemp2=xtemp2, xtemp3=xtemp3)
	return(results)
}


## ***************************************************************************************************

WCon<-function(Zmat, MomentMat,X.reg){
	dimL<-nrow(Zmat)	
	dimN<-ncol(Zmat)
	dimt <-nrow(MomentMat)+1
	Wtemp1<-matrix(0,nrow=dimL, ncol=dimL)
	Wtemp2<-matrix(0,nrow=dimL, ncol=dimL)
	if(X.reg==TRUE){ momentvec<-rowSums(MomentMat) - MomentMat[,2]}
	if(X.reg==FALSE){ momentvec<-rowSums(MomentMat)}
	flagr2 <- 1
	flagz1 <- 1
	flagz2 <- 1
	for(kk2 in 2:dimt){
		mnum <- momentvec[kk2-1]
		if(mnum == 0){}
		if(mnum != 0){
			if(flagr2==1){
				if(!is.matrix(Zmat[1:mnum])){
					Ztemp1 <- matrix(Zmat[1:mnum,], nrow=mnum, ncol=ncol(Zmat))
				}else{
					Ztemp1 <- Zmat[1:mnum,]
				}
				Wtemp2[1:mnum, 1:mnum]<- (1/dimN)*(Ztemp1%*%t(Ztemp1))
				flagr2<- flagr2 + mnum
				flagz2<- flagz2 + mnum
			}else{ mnum1 <- momentvec[kk2-2]
				Ztemp1 <- Zmat[(flagz1:(flagz1 + mnum1 - 1)),]
				Ztemp2 <- Zmat[(flagz2:(flagz2 + mnum - 1)),]
				flagz1 <- flagz1 + mnum1
				flagz2 <- flagz2 + mnum
				Wtemp1[(flagr2 - mnum1):(flagr2 - 1), flagr2:(flagr2 + mnum - 1)] <- -(1/dimN)*(Ztemp1%*% t(Ztemp2))
				Wtemp2[flagr2:(flagr2 + mnum - 1), flagr2:(flagr2 + mnum -1)] <- (1/dimN)*(Ztemp2 %*% t(Ztemp2))
				flagr2 <- flagr2 + mnum
			}
		}
	}
	Wnout <- pinv(Wtemp1 + t(Wtemp1) + 2*Wtemp2)
	return(Wnout)
}

## ***************************************************************************************************

GMMcal<-function(g1nbar,rth,Zmat, xtemp1, xtemp2, xtemp3, qtemp1, qtemp2, Wn){
	dimN<-ncol(qtemp1)
	dimL<-nrow(qtemp1)
	k1<-ncol(xtemp2)/dimN
    #
	indtemp2 <- 1*(qtemp1 > rth)
	indtemp3 <- 1*(qtemp2 > rth)
	#
	g2nbar <- matrix(0, nrow=dimL, ncol= (2*k1 +1))
	g2nbar[, (k1+1)] <- rowSums(Zmat*(indtemp2 - indtemp3))/dimN
	#
	for(kk3 in 1:k1){
		g2nbar[,kk3] <- rowSums(Zmat*xtemp1[,(((kk3-1)*dimN+1):(kk3*dimN))])/dimN
		g2nbar[,(k1+1+kk3)] <- rowSums(Zmat*(xtemp2[,(((kk3-1)*dimN+1):(kk3*dimN))]*indtemp2 -
					xtemp3[,((kk3-1)*dimN+1):(kk3*dimN)]*indtemp3))/dimN
	}
	mattosolve<- (t(g2nbar)%*% Wn %*% g2nbar)
	GMMest <- pinv(mattosolve)%*% (t(g2nbar)%*% Wn %*% g1nbar)	
	rownames(GMMest)<-c(paste("beta.X",seq(1,k1),sep=""),"delta.c",paste("delta.X",seq(1,k1),sep=""))
	#
	gnbar <- g1nbar - g2nbar%*%GMMest
    #
	reslist = list(g2nbar=g2nbar,GMMest=GMMest, gnbar=gnbar)
	return(reslist)
}

## ***************************************************************************************************

GMMcal.res<-function(g1nbar,rth,Zmat, qtemp1, qtemp2, Wn){
	dimN<-ncol(qtemp1)
	dimL<-nrow(qtemp1)
    #
	indtemp2 <- 1*(qtemp1 > rth)
	indtemp3 <- 1*(qtemp2 > rth)
	#
	g2nbar <- matrix(0, nrow=dimL, ncol=1)
	g2nbar[, 1] <- rowSums(Zmat*(indtemp2 - indtemp3))/dimN
	#
	mattosolve<- (t(g2nbar)%*% Wn %*% g2nbar)
	GMMest <- pinv(mattosolve)%*% (t(g2nbar)%*% Wn %*% g1nbar)	
	rownames(GMMest)<-c("delta.c")
	#
	gnbar <- g1nbar - g2nbar%*%GMMest
    #
	reslist = list(g2nbar=g2nbar,GMMest=GMMest, gnbar=gnbar)
	return(reslist)
}

## ***************************************************************************************************

epcal <-function(rhat, xtemp1, xtemp2, xtemp3, qtemp1, qtemp2,ymatrep, GMMest){
	dimL <- nrow(qtemp1)
	dimN <- ncol(qtemp1)
	k1 <- ncol(xtemp2)/dimN
	epcollect <- matrix(0,nrow=dimL,ncol=dimN)
	indtemp2 <- 1*(qtemp1 > rhat)
	indtemp3 <- 1*(qtemp2 > rhat)
	indvec <- seq(0, (k1-1)*dimN, length.out= k1 )
	for(ii2 in 1:dimL){
		for(ii3 in 1:dimN){		
			ii3_rep <- matrix(ii3, k1, 1)
			A <- xtemp1[ii2, (indvec+ii3_rep)]		
			B <- indtemp2[ii2, ii3] - indtemp3[ii2, ii3]
			C <- indtemp2[ii2, ii3]*xtemp2[ii2, indvec+ii3_rep] - 	
			indtemp3[ii2, ii3]*xtemp3[ii2, indvec+ii3_rep]
			epcollect[ii2, ii3] = ymatrep[ii2, ii3] - matrix(GMMest, nrow=1)%*% matrix(c(A,B,C), ncol=1)
		}
	}	
	return(epcollect)
}

## ***************************************************************************************************

W2Con <- function(epcollect, Zmat){
	dimL<-nrow(Zmat)
	dimN<-ncol(Zmat)
	gtemp1 <- matrix(0,nrow=dimL, ncol=dimL)
	gtemp2 <- matrix(0,nrow=dimL, ncol=1)
	for (ii3 in 1:dimN) {
		gihat <- (epcollect[, ii3])*(Zmat[, ii3])
		gtemp1 <- gtemp1 + gihat%*%t(gihat)
		gtemp2 <- gtemp2 + gihat
	}
	mattoinv1 <- (1/dimN)*gtemp1 - (1/(dimN^2))*(gtemp2%*%t(gtemp2))
	Wn2 <- pinv(mattoinv1)
	return(Wn2)
}

## ***************************************************************************************************

covmatfun <- function(GMMest2, rhat2, qmat, Zmat, xtemp1, xtemp2,
xtemp3, qtemp1, qtemp2, h0, omegahat){
	dimL <- nrow(Zmat)
	dimN <- ncol(Zmat)
	# Silverman s rule of thumb + rule of Kernell
	# h0 = factor of the bandwidth
	hband <- h0 * 1.06 * dimN^(-0.2) * sqrt(var(as.vector(t(qmat))))
	k1 <- ncol(xtemp1)/dimN
    #
	Gb <- matrix(0, dimL, k1)
	Gd <- matrix(0, dimL, (k1+1))
	Gr <- matrix(0, dimL, 1)	
	#
	indtemp2 <- 1*(qtemp1 > rhat2)
	indtemp3 <- 1*(qtemp2 > rhat2)
	#
	g2nbar <- matrix(0, dimL,(2*k1+1))
	g2nbar[,(k1 + 1)] <- rowSums(Zmat*(indtemp2 - indtemp3))/dimN
	for(ii6 in 1:k1){
		tempx <- xtemp1[,(((ii6-1)*dimN+1): (ii6*dimN))]
		Gb[,ii6] <- (- rowSums(Zmat*xtemp1[, (((ii6-1)*dimN+1) : (ii6*dimN))])/dimN )
		g2nbar[,ii6] <- rowSums(Zmat*tempx)/dimN
		g2nbar[,(k1+1+ii6)] <- rowSums(Zmat*(xtemp2[,(((ii6-1)*dimN+1):(ii6*dimN))]*indtemp2
											- xtemp3[, (((ii6-1)*dimN+1):(ii6*dimN))]*indtemp3))/dimN
	}
	Gd[,1] <- (-g2nbar[, (k1 + 1)])
	Gd[, (2:(k1+1))] <- ( - g2nbar[, ((k1+2):(2*k1+1))] )
	for(ii7 in 1:dimN){
		indvec <- dimN*seq(0,(k1-1)) + ii7
		temprc1 <- cbind(matrix(1,dimL, 1), xtemp3[, indvec])
		part1<- dnorm((rhat2 - qtemp2[, ii7])/hband)
		temprc2 <- kronecker(part1,matrix(1,1,(k1+1)))
		temprc3 <- cbind(matrix(1,dimL, 1), xtemp2[, indvec])
		part2 <-  dnorm((rhat2 - qtemp1[, ii7])/hband)
		temprc4 <- kronecker(part2,matrix(1,1,(k1+1)))
		tempz <- (temprc1*temprc2 - temprc3*temprc4) %*% GMMest2[((k1+1):(2*k1+1))]
		Gr <- Gr + Zmat[, ii7]*tempz
	}
	Gr <- Gr/(dimN*hband)
	G <- cbind(Gb, Gd, Gr)
	covmat <- chol2inv(t(G)%*%omegahat%*%G)/dimN
	return(covmat)
}

## ***************************************************************************************************

supwaldcal <- function(ymatrepb, Zmatb, xtemp1b, xtemp2b, xtemp3b, qtemp1b, 
qtemp2b, MomentMat, GridTh,  ngrid, dimN){
	k1 <- ncol(xtemp1b)/dimN
	# We separate delta_hat from the GMM_estimate
	# Note that delta_hat includes constant difference
	dltlen <- k1 + 1
	covwaldrec <- matrix(0,nrow=ngrid*dltlen, ncol=dltlen)
	# To record sigma_hat for each r
    #
	# 1st-step GMM estimate 
	Wnb <- WCon(Zmatb, MomentMat, X.reg=TRUE)
	g1nbarb <- rowSums(Zmatb*ymatrepb)/dimN 
	waldvec <- matrix(0,nrow= ngrid, ncol=1)
	for(bb in 1: ngrid){
		rtt <- GridTh[bb]
		GMMresb <- GMMcal(g1nbarb,rtt,Zmatb, xtemp1b, xtemp2b, xtemp3b, qtemp1b, qtemp2b, Wnb)
		g2nbarb <- GMMresb$g2nbar
		gnbarb <- GMMresb$gnbar
		GMMestb <- GMMresb$GMMest
		epcollectb <- epcal(rtt, xtemp1b, xtemp2b, xtemp3b, qtemp1b, qtemp2b,
							ymatrepb, GMMestb)	
		Wn2b <- W2Con(epcollectb, Zmatb)
		GMMres2b <- GMMcal(g1nbarb,rtt,Zmatb, xtemp1b, xtemp2b, xtemp3b, qtemp1b, qtemp2b, Wn2b)
		g2nbar2b <- GMMresb$g2nbar
		gnbar2b <- GMMresb$gnbar
		GMMest2b <- GMMresb$GMMest
		epcollect2b <- epcal(rtt, xtemp1b, xtemp2b, xtemp3b, qtemp1b, qtemp2b,
							  ymatrepb, GMMest2b)
		omegahatb <- W2Con(epcollect2b, Zmatb)
		Vs <- covconb(rtt, Zmatb, xtemp1b, xtemp2b, xtemp3b,qtemp1b, qtemp2b, omegahatb)
		sizec <- ncol(Vs) # size of (Vs transpose times Vs)
		dlthat <- GMMest2b[(sizec-dltlen+1) : sizec]
		selm <- matrix(0,nrow=dltlen, ncol=sizec) 
		selm[, (sizec-dltlen+1):sizec] <- diag(dltlen)
		# To select deltacovmat from the whole (Vs transpose times Vs)
		covwald <- selm %*% pinv(t(Vs) %*% Vs) %*% t(selm)
		waldb <-  dimN * t(dlthat) %*% pinv(covwald) %*% dlthat
		waldvec[bb] <- waldb
		covwaldrec[((bb-1)*dltlen+1): (bb*dltlen), ] <- covwald
	}
	supwald = max(waldvec) # this is the desired sup Wald stat.
	reslist <- list(supwald=supwald,covwaldrec=covwaldrec)
	return(reslist)
}

## ***************************************************************************************************
covconb <- function(rhat2, Zmat, xtemp1, xtemp2, xtemp3, qtemp1, qtemp2, omegahat){
	#This function estimates covariance matrix for bootstrapping section
	dimL = nrow(Zmat)
	dimN = ncol(Zmat)
	#h_band = h_0 times 1.06 times N power (-0.2) times sqrt(variance(colshape(q_mat, 1))) // Not needed here
	k1 = ncol(xtemp1)/dimN
	Gb = matrix(0, nrow=dimL, ncol=k1)
	Gd = matrix(0, nrow=dimL, ncol=(k1+1))
	indtemp2 <- 1*(qtemp1 > rhat2)
	indtemp3 <- 1*(qtemp2 > rhat2)
	# This two (L times N) 0-1 matrices contain information whether q_{it} or q_{it-1} is 
	# larger than the threshold r_th
	g2nbar = matrix(0, nrow=dimL, ncol=(2*k1 + 1)) 
	g2nbar[, (k1 + 1)] = rowSums(Zmat*(indtemp2 - indtemp3))/dimN
	for(ii6 in 1:k1) {
		tempx = xtemp1[, ((ii6-1)*dimN+1):(ii6*dimN)] 
		Gb[, ii6] = -rowSums(Zmat*xtemp1[, ((ii6-1)*dimN+1):(ii6*dimN)])/dimN
		g2nbar[, ii6] = rowSums(Zmat*tempx)/dimN
		g2nbar[, (k1 + 1 + ii6)] = rowSums(Zmat*(xtemp2[, ((ii6-1)*dimN+1):(ii6*dimN)]*indtemp2 - xtemp3[, ((ii6-1)*dimN+1):(ii6*dimN)]*indtemp3))/dimN
	}
	Gd[, 1] = -g2nbar[, (k1+1)]
	Gd[, 2:(k1+1)] = -g2nbar[, (k1+2):(2*k1 + 1)]
	cholDe  <- suppressWarnings({chol(omegahat, pivot=TRUE)})
	Vs <- t(cholDe) %*% cbind(Gb, Gd)
	return(Vs)
}

## ***************************************************************************************************

supwaldfb <- function(Zmatb, xtemp1b, xtemp2b, xtemp3b, qtemp1b,qtemp2b, epcollect2, Wn2, covwaldrec, GridTh, ngrid, dimN){ 
	#
	k1 = ncol(xtemp1b)/dimN
	dimL = nrow(Wn2)
	dltlen = k1 + 1	
    #
	normvec <- qnorm(runif(dimL))
	beta. <- matrix(normvec, nrow=dimL, ncol=1) # For bootstrap weighting		
	g1nbarb <- (rowSums(Zmatb*epcollect2)*beta.)/dimN 
    #
	waldvec <- matrix(0, nrow=ngrid, ncol=1)
	for (bb in 1:ngrid){
		rtt = GridTh[bb]
		GMMcalres <- GMMcal(g1nbarb, rtt, Zmatb, xtemp1b, xtemp2b, xtemp3b, qtemp1b, qtemp2b, Wn2)
		GMMestb <- GMMcalres$GMMest
		gnbarb <- GMMcalres$gnbar 
		g2nbarb <- GMMcalres$g2nbar
        #
		sizec <- nrow(GMMestb)
		dlthat <- GMMestb[(sizec-dltlen+1):(sizec)]
	     #
		covwald <- covwaldrec[((bb-1)*dltlen+1) : (bb*dltlen), ]
		waldb <-  dimN %*% t(dlthat) %*% pinv(covwald) %*% dlthat
		waldvec[bb] <- waldb
	}
	supwaldb <- max(waldvec)
	return(supwaldb)
}

## ***************************************************************************************************

supwaldcal.con <- function(ymatrepb, Zmatb, qtemp1b, 
qtemp2b, MomentMat, GridTh,  ngrid, dimN){
	k1 <- 1
	# We separate delta_hat from the GMM_estimate
	# Note that delta_hat includes constant difference	
	dltlen <- k1 + 1
	covwaldrec <- matrix(0,nrow=ngrid*dltlen, ncol=dltlen)
	# To record sigma_hat for each r
    #
	# 1st-step GMM estimate 
	Wnb <- WCon(Zmatb, MomentMat, X.reg=FALSE)
	g1nbarb <- rowSums(Zmatb*ymatrepb)/dimN 
    #
	waldvec <- matrix(0,nrow= ngrid, ncol=1)
	for (bb in 1: ngrid) {
		rtt <- GridTh[bb]
		GMMresb <- GMMcal.res(g1nbarb,rtt,Zmatb, qtemp1b, qtemp2b, Wnb)
		g2nbarb <- GMMresb$g2nbar
		gnbarb <- GMMresb$gnbar
		GMMestb <- GMMresb$GMMest
		indtemp2 <- 1*(qtemp1b > rtt)
		indtemp3 <- 1*(qtemp2b > rtt)
		epcollectb <- ymatrepb - GMMestb[1,1]*(indtemp2 - indtemp3)
		Wn2b <- W2Con(epcollectb, Zmatb)
		GMMres2b <- GMMcal.res(g1nbarb,rtt,Zmatb, qtemp1b, qtemp2b, Wn2b)
		g2nbar2b <- GMMresb$g2nbar
		gnbar2b <- GMMresb$gnbar
		GMMest2b <- GMMresb$GMMest
		indtemp2b <- 1*(qtemp1b > rtt)
		indtemp3b <- 1*(qtemp2b > rtt)
		epcollect2b <- ymatrepb - GMMest2b[1,1]*(indtemp2b - indtemp3b)
		omegahatb <- W2Con(epcollect2b, Zmatb)
		Vs <- covconb.con(rtt, Zmatb,	qtemp1b, qtemp2b, omegahatb)
		sizec <- ncol(Vs) # size of (Vs transpose * Vs)
		dlthat <- GMMest2b[(sizec-dltlen+1) : sizec]
		selm <- as.matrix(c(1,0), nrow=dltlen, ncol=sizec)
		covwald <- selm %*% pinv(t(Vs) %*% Vs) %*% t(selm)
		waldb <-   dimN * (dlthat^2) * (covwald[1,1]^(-1))
		waldvec[bb] <- waldb
		covwaldrec[((bb-1)*dltlen+1): (bb*dltlen), ] <- covwald
	}
	supwald = max(waldvec)
	reslist <- list(supwald=supwald,covwaldrec=covwaldrec)
	return(reslist)
}

## ***************************************************************************************************

covconb.con <- function(rhat2, Zmat, qtemp1, qtemp2, omegahat){
	#This function estimates covariance matrix for bootstrapping section.
	dimL = nrow(Zmat)
	dimN = ncol(Zmat)
	Gd = matrix(0, nrow=dimL, ncol=1)
	indtemp2 <- 1*(qtemp1 > rhat2)
	indtemp3 <- 1*(qtemp2 > rhat2)
	g2nbar <- matrix(0,nrow=dimL, ncol=1)
	g2nbar[, 1] <- rowSums(Zmat*(indtemp2 - indtemp3))/dimN
	Gd[, 1] <- - g2nbar[, 1]
	cholDe  <- suppressWarnings({chol(omegahat, pivot=TRUE)})
	Vs <- t(cholDe) %*% Gd
	return(Vs)
}

## ***************************************************************************************************

supwaldfb.con <- function(Zmatb, qtemp1b,qtemp2b, epcollect2, Wn2, covwaldrec, GridTh, ngrid, dimN){ 
	k1 = 1
	dimL = nrow(Wn2)
	dltlen = k1 + 1	
	normvec <- qnorm(runif(dimL))
	beta. <- matrix(normvec, nrow=dimL, ncol=1) # For bootstrap weighting		
	g1nbarb <- (rowSums(Zmatb*epcollect2)*beta.)/dimN 
	waldvec <- matrix(0, nrow=ngrid, ncol=1)
	for (bb in 1:ngrid){
		rtt = GridTh[bb]
		GMMcalres <- GMMcal.res(g1nbarb,rtt,Zmatb, qtemp1b, qtemp2b, Wn2)
		GMMestb <- GMMcalres$GMMest
		gnbarb <- GMMcalres$gnbar 
		g2nbarb <- GMMcalres$g2nbar
		sizec <- nrow(GMMestb)
		dlthat <- GMMestb[(sizec-dltlen+1):(sizec)]
		covwald <- covwaldrec[((bb-1)*dltlen+1) : (bb*dltlen), ]	
		waldb <-  dimN * dlthat^2 * covwald[1,1]^{-1}
		waldvec[bb] <- waldb
	}
	supwaldb <- max(waldvec) 
	return(supwaldb)
}

