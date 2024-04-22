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
setClass("cpoint",
         representation(CP   = "data.frame"
                        ),
         prototype = list(CP   = data.frame()
                        )
)
## ***************************************************************************************************
cpoint <- function(data.,nameI=NA,nameT=NA,nameJ=NA,nameY=NA,thresholds){   
	# data. = one panel dataset (a matrix NxT or a data.frame with columns i,t,j,Y.)
	# thresholds = vector of estimated gamma (one for each j in the panel)
	if(!is.na(nameI) | !is.na(nameT) | !is.na(nameJ) | !is.na(nameY)){
    	#check
		if(!is.data.frame(data.)){stop("data. should be either a matrix with columns representing the times of observation and rows the individuals,
		or a data.frame. In the latter case columns nameI, nameT, nameJ, and nameY should all be specified.")}
		data.<-as.data.frame(data.)
		if(ncol(data.)<(4)){stop("data.frame structure not compatible")}	
		if(any(is.na(data.))){stop("data. should not contain missing values")}
		if(length(nameI)!=1){ stop("unique identifier nameI should be specified")}
		if(length(nameT)!=1){ stop("unique identifier nameT should be specified")}
		if(length(nameJ)!=1){ stop("unique identifier nameJ should be specified")}
		if(length(nameY)!=1){ stop("unique dependent variable nameY should be specified")}
		if(!(nameI %in% colnames(data.))){stop("unable to find the selected nameI variables in the data.")}
		if(!(nameT %in% colnames(data.))){stop("unable to find the selected nameT variables in the data.")}
		if(!(nameJ %in% colnames(data.))){stop("unable to find the selected nameJ variables in the data.")}
		if(!(nameY %in% colnames(data.))){stop("unable to find the selected nameY variable in the data.")}
        #
		# Prepare dataset
		colnames(data.)[which(colnames(data.)==nameI)] <- "i"
		colnames(data.)[which(colnames(data.)==nameJ)] <- "j"
		colnames(data.)[which(colnames(data.)==nameT)] <- "t"
		colnames(data.)[which(colnames(data.)==nameY)] <- "Y"
        #
		id.vec <- unique(data.$i)
		j.vec <- unique(data.$j)
		t_cp <- as.data.frame(matrix(,ncol=(length(j.vec)+1),nrow=length(id.vec)))
		t_cp[,1] <- id.vec
        #
		if(length(thresholds)!=length(j.vec)){stop("threshold of incompatible dimension: one threshold per each j should be specified.")}
		if(!is.numeric(thresholds)){stop("threshold should be numeric.")}
        #
		name.level <- paste("CP_","j",j.vec,sep="")
		colnames(t_cp) <- c("i",name.level)
		#
		for(j in 1:length(j.vec)){
			Ydb_sel <- data.[which(data.$j==j.vec[j]),]
			Yseries <- Ydb_sel[,c("i","t","Y")]
			threshold.est <- thresholds[j]
            #
			for(id in id.vec){
				simul_i<-Yseries[which(Yseries$i==id & Yseries$t>=3),] # DROP t=0,..,3 BECAUSE USED AS INSTRUMENTS
				colnames(simul_i)[3]<-"series"
				if(simul_i$series[1]<threshold.est){
					above<-simul_i$t[which(simul_i$series >= threshold.est)]
					subgroups.above <- split(above, cumsum(c(1L, diff(above) != 1)))
					lengthseries <- lengths(subgroups.above)
					if(sum(lengthseries)==length(lengthseries)){t_cp[id,(j+1)] <- NA}
					if(sum(lengthseries)!=length(lengthseries)){
						index.longer.series <- as.numeric(names(lengthseries[which(lengthseries==max(lengthseries))]))
						if(length(index.longer.series)!=1){
							t_cp[id,(j+1)] <- NA
						}
						if(length(index.longer.series)==1){
							longer.series <- subgroups.above[[index.longer.series]]
							t_cp[id,(j+1)] <- longer.series[1] - 1
						}
					}
				}
				if(simul_i$series[1]>=threshold.est){
					below<-simul_i$t[which(simul_i$series < threshold.est)]		
					subgroups.below <- split(below, cumsum(c(1L, diff(below) != 1)))
					lengthseries <- lengths(subgroups.below)
					if(sum(lengthseries)==length(lengthseries)){t_cp[id,(j+1)] <- NA}
					if(sum(lengthseries)!=length(lengthseries)){
						index.longer.series <- as.numeric(names(lengthseries[which(lengthseries==max(lengthseries))]))
						if(length(index.longer.series)!=1){
							t_cp[id,(j+1)] <- NA
						}
						if(length(index.longer.series)==1){
							longer.series <- subgroups.below[[index.longer.series]]
							t_cp[id,(j+1)] <- longer.series[1] - 1
						}
					}
				}
			}
		}
	}
	if(is.na(nameI) & is.na(nameT) & is.na(nameJ) & is.na(nameY)){
		#check
		if(!is.matrix(data.)){stop("data. should be either a matrix with columns representing the times of observation and rows the individuals,
		or a data.frame. In the latter case columns nameI, nameT, nameJ, and nameY should all be specified.")}
		if(length(thresholds)!=1){stop("threshold of incompatible dimension: If data. is a matrix, one threshold value should be specified.")}
		if(!is.numeric(thresholds)){stop("threshold should be numeric.")}
		if(ncol(data.)<5){stop("data. should contain at least 6 times.")}
		if(is.null(colnames(data.))){colnames(data.) <- seq(1, ncol(data.),1)}
		t_cp <- as.data.frame(matrix(,ncol=2,nrow=nrow(data.)))
		if(!is.null(rownames(data.))){t_cp[,1] <- rownames(data.)}
		if(is.null(rownames(data.))){t_cp[,1] <- seq(1,nrow(data.),1)}
		colnames(t_cp) <- c("i","CP")
		for(id in 1:nrow(data.)){
			simul_i <- data.[id,3:ncol(data.)] # DROP t=0,..,3 BECAUSE USED AS INSTRUMENTS
			simul_i <- as.data.frame(cbind(as.numeric(colnames(data.)[3:ncol(data.)]),simul_i))
			colnames(simul_i) <- c("t","series")
            threshold.est <- thresholds
				if(simul_i$series[1]<threshold.est){
					above<-simul_i$t[which(simul_i$series >= threshold.est)]
					subgroups.above <- split(above, cumsum(c(1L, diff(above) != 1)))
					lengthseries <- lengths(subgroups.above)
					if(sum(lengthseries)==length(lengthseries)){t_cp[id,2] <- NA}
					if(sum(lengthseries)!=length(lengthseries)){
						index.longer.series <- as.numeric(names(lengthseries[which(lengthseries==max(lengthseries))]))
						if(length(index.longer.series)!=1){
							t_cp[id,2] <- NA
						}
						if(length(index.longer.series)==1){
							longer.series <- subgroups.above[[index.longer.series]]
							t_cp[id,2] <- longer.series[1] - 1
						}
					}
				}
				if(simul_i$series[1]>=threshold.est){
					below<- simul_i$t[which(simul_i$series < threshold.est)]
					subgroups.below <- split(below, cumsum(c(1L, diff(below) != 1)))
					lengthseries <- lengths(subgroups.below)
					if(sum(lengthseries)==length(lengthseries)){t_cp[id,2] <- NA}
					if(sum(lengthseries)!=length(lengthseries)){
						index.longer.series <- as.numeric(names(lengthseries[which(lengthseries==max(lengthseries))]))
						if(length(index.longer.series)!=1){
							t_cp[id,2] <- NA
						}
						if(length(index.longer.series)==1){
							longer.series <- subgroups.below[[index.longer.series]]
							t_cp[id,2] <- longer.series[1] - 1
						}
					}
				}
			}
		
		}
    out <- new("cpoint")
    out@CP <- t_cp;
    return(out);
}
