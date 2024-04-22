\name{simptm}
\alias{simptm}
\title{Data simulation from two- or three-way panel threshold regression model.}
\description{
Random generation of data from two- or three-way panel threshold regression model with or without a time-varying regressor.
}
\usage{
simptm(n, T., J=2, CP, gamma.=c(0,0), 
        phi_c=matrix(c(-1,1,-0.7,1.8), nrow=2, byrow=TRUE), 
        phi_X=matrix(c(-0.2,0.2,-0.5,0.8), nrow=2, byrow=TRUE), 
        sigmau=1, parAR=c(0.7,0.5), B=200, seedstart=1) 
}
\arguments{
	\item{n}{number of statistical units.}
	\item{T.}{number of times.}
	\item{J}{number of third way's values.}
	\item{CP}{vector of times of regime switch (one per each \eqn{j}). If vector of length 1, the same change point is taken for all the \eqn{j}s.}
	\item{gamma.}{vector of length J of threshold values.}
	\item{phi_c}{matrix \eqn{J\times 2} of the constant parameters of the two regimes; first column for lower (upper) regime and second column for upper (lower) regime.}
	\item{phi_X}{matrix \eqn{J\times 2} of the regressor parameter of the two regimes; first column for lower (upper) regime and second column for upper (lower) regime.}
	\item{sigmau}{possible constant to be applied to the error term distributed as a Gaussian white noise.}
	\item{B}{number of datasets to be drawn.}
	\item{parAR}{vector of autoregressive parameter(s) (one for each \eqn{j}) to generate a (exogenous) time-varying regressor.}
	\item{seedstart}{number of the initial seed.} 
} 
\details{ 
\code{simptm} generates \code{B} datasets from a two- or three-way panel threshold regression model with \code{n} statistical units, \code{T.} times of observations and \code{J} values/levels for the third dimension. The data generating process is constituted by two regimes with a change point at time \code{CP} and a threshold value \eqn{y_{t-1} = \gamma_j} for each value/level \eqn{j} of the third way. 
The two regimes are defined by the constant parameters \code{phi_c} or \code{phi_c} and parameters for the time-varying regressor \code{phi_X}. In the current version of the package, it is thus possible to generate data from a model without regressors or with a time-varying regressor. The (exogenous) time-varying regressor is assumed to be distributed as an autoregressive stochastic process with the \eqn{j}-th parameter of the \code{parAR}. The error component is generated from a Gaussian white noise and can be rescaled through the constant sigmau. 
} 

\value{
An object of S4 class "simptm", which is a list of simulated data matrices.
}
\references{ 
Di Lascio, F.M.L. and Perazzini, S. (202x) A three-way dynamic panel threshold regression model for change point detection in bioimpedance data. WP BEMPS <https://repec.unibz.it/bemps104.pdf>.

Di Lascio, F.M.L. and Perazzini, S. (2022) Change point detection in fruit bioimpedance using a three-way panel model. Book of short papers - SIS2022, p.1184-1189, Eds. A. Balzanella, M. Bini, C. Cavicchia, R. Verde. Pearson. ISBN: 978-88-9193-231-0. 

Seo, M.H. and Shin, Y. (2016) Dynamic panels with threshold effect and endogeneity, Journal of Econometrics, 195(2), p.169-186. 
} 
\note{
The estimation method requires at least \eqn{T.\geq6}: four lags of the dependent and independent variables are used as instruments, and two more are necessary to identify the 
regime switch (i.e., one per regime). The output is a list of \code{B} dataframes. Each dataframe contains columns: \eqn{i}, \eqn{j}, \eqn{t}, \eqn{Y} and possibly \eqn{X}. 
} 
\author{
Francesca Marta Lilja Di Lascio <marta.dilascio@unibz.it>

Selene Perazzini <selene.perazzini@alumni.imtlucca.it>
}
\seealso{
	See also \code{\link{ptm2}} and \code{\link{ptm3}}.
}
\examples{

## NOT RUN
#
## Simulation of 10 two-way panels  
## y_{it} = (-1-0.2*x_{it})1(y_{it-1}<=0) + (1+0.2*x_{it})1(y_{it-1}>0) 
## with 50 statistical units observed for 50 times with change point at 
## time 20, autoregressive parameter 0.7, and sigmau=1. 
#
DB1 <- simptm(n=50, T.=11, J=1, CP=8, gamma.=0, phi_c=matrix(c(-1,1), 
     nrow=1, byrow=TRUE), phi_X=matrix(c(-0.2,0.2), nrow=1, byrow=TRUE), 
     sigmau=1, parAR=0.7, B=10, seedstart=1)
#head(DB1@simulation)
#str(DB1@simulation)
    
## Simulation of 10 three-way panels 
## y_{i1t} = (-1-0.2*x_{i1t})1(y_{i1t-1}<=0) + 
##           (1+0.2*x_{i1t})1(y_{i1t-1}>0)
## y_{i2t} = (-0.7-0.5*x_{i2t})1(y_{i2t-1}<=0) + 
##           (1.8+0.8*x_{i2t})1(y_{i2t-1}>0)
## with 50 statistical units, 50 times, J=2, change point corresponding to 
## time 20, with autoregressive parameter 0.7, and sigmau=1. 
#
#DB2 <- simptm(n=50, T.=20, J=2, CP=10, gamma.=c(0,0), 
#    phi_c=matrix(c(-1,1,-0.7,1.8), nrow=2, byrow=TRUE), 
#    phi_X=matrix(c(-0.2,0.2,-0.5,0.8), nrow=2, byrow=TRUE), sigmau=1, 
#    parAR=rep(0.7,2), B=5, seedstart=1)
#head(DB2@simulation)
#str(DB2@simulation)
#
##
}
  
