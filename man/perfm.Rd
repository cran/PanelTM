\name{perfm}
\alias{perfm}
\title{Performance measures of the (two-way or three-way) panel threshold regression model estimated on simulated data.}
\description{
Computes the relative bias, the bias, the root mean squared error, and the relative root mean squared error.
}
\usage{ perfm(truepar,estpar)
}
\arguments{
	\item{truepar}{True value of the parameter.}
	\item{estpar}{Vector of the estimated parameters (one value per each simulated dataset).}
}
\details{
    \code{perfm} computes the following performance measures: the relative bias, the bias, the root mean squared error, and the relative root mean squared error.
}
\value{
    An object of S4 class "perfm", which is a \code{list} with the following values: the true value of the parameters, the relative bias, the bias, the root relative mean squared error, the relative root relative means squared error.
}
\references{
Di Lascio, F.M.L. and Perazzini, S. (202x) A three-way dynamic panel threshold regression model for change point detection in bioimpedance data. Submitted.

Di Lascio, F.M.L. and Perazzini, S. (2022) Change point detection in fruit bioimpedance using a three-way panel model. Book of short papers - SIS2022, p.1184-1189, Eds. A. Balzanella, M. Bini, C. Cavicchia, R. Verde. Pearson. ISBN: 978-88-9193-231-0.
}
\note{
    The vector of the estimated parameters can contain missing values that are excluded from the computation of the performance mesaures.
}
\author{
Francesca Marta Lilja Di Lascio <marta.dilascio@unibz.it>

Selene Perazzini <selene.perazzini@alumni.imtlucca.it>
}
\seealso{
	See also \code{\link{ptm2}}, \code{\link{ptm3}} and \code{\link{simptm}}.
}
\examples{
# Simulation of the two-way panel threshold model 
# y_{it} = (-0.7-0.5*x_{it})1(y_{it-1}<=0) + (1.8+0.8*x_{it})1(y_{it-1}>0)
# Two datasets with 20 individuals observed over 15 times with change  
# point at time 12 and an autoregressive error term (with parameter 0.7 
# and unit error variance).

sims <- simptm(n=20, T.=15, J=1, CP=12, gamma=c(0),  
               phi_c=matrix(c(-0.7,1.8), nrow=1, byrow=TRUE), 
               phi_X=matrix(c(-0.5,0.8), nrow=1, byrow=TRUE), sigmau=1, 
               parAR=0.2, B=2, seedstart=1)

# Estimate the two simulated two-way panel threshold models

estimates1 <- ptm3(sims@simulation[[1]], nameI="i", nameT="t", nameJ="j", 
                   nameY="Y", nameXexo="X1", nameXendo=NULL, nameTV=NULL, 
                   nameIV=NULL, trimrate=0.4, ngrid=100, h0=1.5, 
                   Iweight=FALSE, test.lin=FALSE)

estimates2 <- ptm3(sims@simulation[[2]], nameI="i", nameT="t", nameJ="j", 
                   nameY="Y", nameXexo="X1", nameXendo=NULL, nameTV=NULL, 
                   nameIV=NULL, trimrate=0.4, ngrid=100, h0=1.5,
                   Iweight=FALSE, test.lin=FALSE)

# Compute performance measures for the parameter beta of the exogenous 
# regressor (for details see refs)
pm1 <- perfm(truepar=-0.5, estpar=c(estimates1@param[,2],
             estimates2@param[,2]))
pm1

# Compute performance measures for the delta corresponding to the
# difference between the constants of the model (for details see refs)
pm2 <- perfm(truepar=c(1.8+0.7), estpar=c(estimates1@param[,4], 
             estimates2@param[,4]))
pm2

# Compute performance measures for the delta corresponding to the
# difference between the regressors coefficients of the model 
# (for details see refs)
pm3 <- perfm(truepar=c(0.8+0.5), estpar=c(estimates1@param[,6],
             estimates2@param[,6]))
pm3
}
