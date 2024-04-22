\name{cpoint}
\alias{cpoint}
\title{Computation of the regime switch time in a (two-way or three-way) panel threshold regression model.}
\description{
A measure for computing the time at which a regime changes, i.e. the change point.
}
\usage{
cpoint(data., nameI=NA, nameT=NA, nameJ=NA, nameY=NA, thresholds)
}
\arguments{
  \item{data.}{A (two-way or three-way) panel dataset. It can be either a \eqn{N\times T} matrix of the dependent variable(s) where the row represents the statistical unit (\eqn{i})
		and the column the time (\eqn{t}) or a \code{data.frame} containing columns \code{nameI}, \code{nameT}, \code{nameJ}, and \code{nameY}, which should all be specified.}
   \item{nameI}{If data. is a \code{data.frame}, the name of the (numerical) variable that identifies the statistical units. Otherwise, it is \code{NULL}.}
  \item{nameT}{If data. is a \code{data.frame}, the name of the (numerical) time variable. Otherwise, it is \code{NULL}.}
  \item{nameJ}{If data. is a \code{data.frame}, the name of the (numerical/categorical) variable that indicates the third dimension. Otherwise, it is \code{NULL}.}
  \item{nameY}{If data. is a \code{data.frame}, the name of the (numerical) dependent variable. Otherwise, it is \code{NULL}.}
  \item{thresholds}{The vector of estimated threshold value(s). If data. is a \code{matrix}, the object should contain one value \eqn{\gamma}. If data. is a \code{data.frame}, 
		the vector should contain one threshold value per each value of j, i.e. \eqn{\gamma_j}.}
}
\details{
\code{cpoint} computes the change point on the dependent variable (\eqn{Y}) for each statistical unit i given the time series of \eqn{Y} and the estimated threshold value \eqn{\gamma_{ij}}.
If the observed value \eqn{Y} is greater (lower) than the estimated threshold parameter, \eqn{t=4} belongs to the upper (lower) regime, 
and the change point is identified as the time \eqn{t} after which the longest sequence of days in the lower (upper) regime is observed.
}
\value{
An object of S4 class "cpoint", which is a \code{matrix} reporting the change point per each statistical unit \eqn{i} and each \eqn{j}.
}
\references{
Di Lascio, F.M.L. and Perazzini, S. (202x) A three-way dynamic panel threshold regression model for change point detection in bioimpedance data. WP BEMPS <https://repec.unibz.it/bemps104.pdf>.

Di Lascio, F.M.L. and Perazzini, S. (2022) Change point detection in fruit bioimpedance using a three-way panel model. Book of short papers - SIS2022, p.1184-1189, Eds. A. Balzanella, M. Bini, C. Cavicchia, R. Verde. Pearson. ISBN: 978-88-9193-231-0.
}
\note{
The change point is detected in \eqn{t>4} since the criterion has been inspired by the (two-way or three-way) dynamic panel threshold regression model that 
requires at least 6 times of observation: four lags of the dependent and independent variables are used as instruments and two more are necessary to identify the 
regime switch (i.e., one per regime).
}
\author{
Francesca Marta Lilja Di Lascio <marta.dilascio@unibz.it>

Selene Perazzini <selene.perazzini@alumni.imtlucca.it>
}
\seealso{
	See also \code{\link{banana}}, \code{\link{ptm2}}, and \code{\link{ptm3}}.
}
\examples{
# Import data
data(banana)

## Example 1: Application to a three-way panel data

# Select data for all the three levels of the third way 
bioimp.df2 <- banana

# Estimation: three-way ptm on bioimpedance data with the lagged 
# dependent variable as transition variable
estimates2 <- ptm3(data.=bioimp.df2, nameI="i", nameT="t",  
                   nameJ="j", nameY="bioimpedance", nameTV=NULL, 
                   nameXexo="weight", nameIV=NULL, trimrate=0.4, 
                   ngrid=100, h0=1.5, Iweight=FALSE, 
                   test.lin=FALSE)

# Change point computation using a data frame (unique input choice)
CP3 <- cpoint(data.=bioimp.df2, nameI="i", nameT="t", nameJ="j", 
              nameY="bioimpedance", 
              thresholds=estimates2@"threshold"[,2])

## NOT RUN           
### Example 2: Application to a two-way panel data
#
## Select data for a specific level of the third way 
#bioimp.df <- banana[which(banana$j==1),]
#
## Prepare input matrix for ptm2
#bioimp.m <-  matrix(bioimp.df$bioimpedance, 
#                    ncol=length(unique(bioimp.df$t)), byrow=TRUE)
#
## Estimation: two-way ptm on bioimpedance data with the lagged 
## dependent variable as transition variable
#estimates <- ptm2(Y=bioimp.m, TV=NULL, Xendo=NULL, Xexo=NULL,  
#                  IV=NULL, trimrate=0.4, ngrid=100, h0=1.5, 
#                  Iweight=FALSE, test.lin=FALSE)
#
## Change point computation using a data frame as input
#CP1 <- cpoint(data.=bioimp.df, nameI="i", nameT="t", nameJ="j", 
#              nameY="bioimpedance", 
#              thresholds=estimates@"threshold"[[1]])
#    
## Change point computation using a data matrix as input
#CP2 <- cpoint(data.=bioimp.m, nameI=NA, nameT=NA, nameJ=NA, 
#              nameY=NA, thresholds=estimates@"threshold"[[1]])
##
}

