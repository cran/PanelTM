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
setClass("perfm",
     	  representation(trueval="numeric",
			rb="numeric",
			nrb="numeric",
			rrmse="numeric",
			rmse="numeric"
                        ),
         prototype = list(trueval=numeric(),
				rb = numeric(),
				nrb = numeric(),
				rrmse = numeric(),
				rmse = numeric()
                        )
         )
## ***************************************************************************************************
perfm <- function(truepar,estpar){
	# truepar = value of the true parameter
	# estpar = vector of estimated parameters
    #
	# Check
	if(length(truepar)!=1){
		stop("truepar must be a unique value")
	}
	if(!is.numeric(truepar)){
		stop("truepar must be numeric")
	}
	if(!is.vector(estpar)){
		stop("estpar must be a vector of length >=1")
	}
	if(!is.numeric(estpar)){
		stop("estpar must be numeric")
	}
	# REMOVE MISSING ESTIMATES
	ep <- estpar[which(!is.na(estpar))]
	# TRUE VALUE
	trueval <- truepar
	# RELATIVE BIAS
	rb <- mean(((ep - truepar)/truepar))
	# BIAS 
	nrb <- mean((ep - truepar))
	# RRMSE
	rrmse <- sqrt(mean(((ep - truepar)/truepar)^2 ))
	# RMSE
	rmse <- sqrt(mean((ep - truepar)^2 ))
	# Show results
    out <- new("perfm")
    out@trueval <- trueval;
    out@rb <- rb;
    out@nrb <- nrb;
    out@rrmse <- rrmse;
    out@rmse <- rmse;
    #
    return(out);
}

