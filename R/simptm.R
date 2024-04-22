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
setClass("simptm",
         representation(simulation="list"),
         prototype = list(simulation=list( ))
         )
## ***************************************************************************************************
simptm <- function(n, T., J=2, CP, gamma.=c(0,0), phi_c = matrix(c(-1, 1, -0.7,1.8), nrow=2, byrow=TRUE), 
            phi_X = matrix(c(-0.2,0.2,-0.5,0.8), nrow=2, byrow=TRUE), sigmau=1, parAR=c(0.7,0.5), 
            B=200, seedstart=1){    
#
# n = number of observation i
# T. = number of times t
# CP = vector of times of regime switch. If length(CP)==1, the same CP is taken for all the js.
# J = number of j
# gamma. = vector of length J of threshold values
# phi_c = matrix Jx2 of the constant parameters of the two regimes: row 1=c(constant lower, constant upper)
# phi_X = matrix Jx2 of the regressor parameters of the two regimes: row 1=c(X lower, X upper)
# sigmau = constant to be applied to the standard gaussian error term
# B = number of replications
# parAR = vector (of length J) of the AR model parameter generating Xs; each element refers to a specific j (j=1,...,J)

	if(n<2 | T.<2){
		stop("Sample size and time series lenght should be greated than 1")
	}
	if(J<=0){
		stop("J should be greater than 0")
	}
	if(any(CP>T.-1 | CP<1)){
		stop("CP incompatible with T.")
	}
	if(length(CP)==1) {
		CP <- rep(CP,J)
	}
	if(length(CP)!=J) {
		stop("Object CP of incompatible dimensions")
	}
	if(length(gamma.)!=J) {
		stop("Length of gamma. object incompatible with J")
	}	
	if(!is.null(phi_c)) {
        if(J!=nrow(phi_c)){
		  stop("Dimension of phi_c object incompatible with J")
	   }
    }
	if(!is.null(phi_X)) {
        if(J!=nrow(phi_X)){
        		stop("Dimension of phi_X object incompatible with J")
	   }
    }
	if(is.null(phi_X) & is.null(phi_c)){
		stop("phi_c or phi_X should be defined")
	}
	if(length(parAR)!=J){
		stop("The parameter p of an AR model should be given for each level j.")
	}
	if(!is.null(phi_c) & !is.null(phi_X)){
		phi <- cbind(phi_c[,1],phi_X[,1:(ncol(phi_X)/2)],phi_c[,2],phi_X[,((ncol(phi_X)/2)+1):ncol(phi_X)])
	}
	if(is.null(phi_c) & !is.null(phi_X)){
		phi <- cbind(rep(0,J),phi_X[,1:(ncol(phi_X)/2)],rep(0,J),phi_X[,((ncol(phi_X)/2)+1):ncol(phi_X)])
	}
	if(!is.null(phi_c) & is.null(phi_X)){
		phi <- phi_c
	}
	# Output of simulation:
	simuldat <- vector("list", length= B)
	# Set seed
	seedflag <- seedstart	
	set.seed(seedflag)
    #
	for(b in 1:B){
        message(b)
		npar <- (ncol(phi))/2   # Number of parameters per regime (both cons and X)
		# Number of X
            nX <- npar - 1 # cost + X
		# Create the matrix of results. This includes i,j,t, Y and Xs.
		simulmat <- as.data.frame(matrix( , nrow=J*n*T., ncol=(4+nX) ))
		simulmat[,1] <- rep((rep(c(1:n), each=T.)), J)
		simulmat[,2] <- rep(1:T., n*J)
		simulmat[,3] <- rep(c(1:J), each=T.*n)
		colnames(simulmat)[1]<- "i"
		colnames(simulmat)[2]<- "t"
		colnames(simulmat)[3]<- "j"
		colnames(simulmat)[4]<- "Y"
		if(nX!=0){colnames(simulmat)[5:ncol(simulmat)]<-paste('X', seq(1,nX), sep="")}
		#Matrices of results of one simulation (one j)
		simulY_j <- c()
		simulX_j <- c()
	    #
		j<-1
		while(j<=J){
			phi_j   <- as.matrix(phi[j,], , )	# Select parameters of jth element of the third way
			gamma_j <- gamma.[j]			# Select the jth element of gamma.
			CP_j <- CP[j]				# Select the jth CP
			seriesj <- matrix( , n, T.)		# Y series
			if(nX>0){
				mat.Xj <- c()			# Matrix containing Xs for the j-th series
			} else { mat.Xj <- NULL }
            #
			i <- 1
			# Generate error, X and Y
			# Check Y wrt gamma and decide whether to keep or drop Y
			while(i<=n){
                # Generate error term
				err <- arima.sim(model = list(order = c(0, 0, 0)), rand.gen=rnorm, n = T.)
				# Generate Xs and Y				
				Y <- rep(0, T.)			
				tempY <- c()
                #
				if(nX == 0){
					tempY[1:CP_j] = phi_j[1] + sigmau*err[1:CP_j]
					tempY[(CP_j+1):T.] = phi_j[2] + sigmau*err[(CP_j+1):T.]				
				}
				if(nX >0){
					# Generate random Xs	(using possible different ARMA models by varying j=1,...J)
                    tempX <- as.matrix(replicate(nX,arima.sim(list(order=c(1,0,0), ar = parAR[j]), n = T.)))	
					# Generate Y and keep series according to gamma.
					tempY[1:CP_j] <- phi_j[1] + tempX[1:CP_j,] %*% as.matrix(phi_j[2:npar], ncol=1) + sigmau*err[1:CP_j]
					tempY[(CP_j+1):T.] <- phi_j[npar+1] + tempX[(CP_j+1):T.,] %*% as.matrix(phi_j[(npar+2):length(phi_j)], ncol=1) + sigmau*err[(CP_j+1):T.]				
				}
				# Check Y versus gamma.
				if( all(tempY[1:CP_j] <= gamma.[j]) & all(tempY[(CP_j+1):T.]> gamma.[j]) ){
					seriesj[i,] <- tempY
					if(nX != 0){
						mat.Xj <- rbind(mat.Xj, tempX)
					}
					i <- i + 1
				}
				if( all(tempY[1:CP_j] > gamma.[j]) & all(tempY[(CP_j+1):T.]<= gamma.[j]) ){
					seriesj[i,] <- tempY
					if(nX != 0){
						mat.Xj <- rbind(mat.Xj, tempX)
					}
					i <- i + 1
				}
			}
			#Estimation check: if unable to find solutions in estimation, discard and simulate the j-th matrices again
			#Set Xs as input of ptm2
			if(nX>0){
				xmat<-c()
				for(jj in 1:nX){
					xtemp <- matrix(mat.Xj[,jj],ncol=T.,byrow=TRUE)
					xmat<- rbind(xmat,xtemp)
				}
			}
            if(nX==0){
                xmat <- NULL
            }
			possibleError <- tryCatch(ptm2(Y=seriesj,TV=NULL,Xexo=xmat,Xendo=NULL,IV=NULL, trimrate=0.4, ngrid=100, h0=1.5, Iweight=FALSE, test.lin=FALSE), error=function(e) e)
			if(inherits(possibleError, "error")){			
			}else{
				# Save results of the i-th simulated individual of j, move to the next j
				simulY_j <- c(simulY_j,t(seriesj))
				simulX_j <- rbind(simulX_j,mat.Xj)
				j <- j+1
			}
		}
		# Save results of all the js of the b simulation
		simulmat[,4] <- simulY_j
		if(nX>0){simulmat[,5:ncol(simulmat)] <- simulX_j}
		simuldat[[b]] <- simulmat
		names(simuldat)[b] <- paste("Data matrix for B=",b, sep="")
		b <- b+1
	}
	# Show results
      out <- new("simptm")
      out@simulation <- simuldat
      return(out)
}
