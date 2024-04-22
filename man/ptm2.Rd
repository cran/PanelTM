\name{ptm2}
\alias{ptm2}
\title{Two-way panel threshold regression model.}
\description{
Two-way panel threshold regression model with the 2-step GMM estimator.
}
\usage{
ptm2(Y, TV=NULL, Xendo=NULL, Xexo=NULL, IV=NULL, trimrate=0.4,  
    ngrid=100, h0=1.5, Iweight=FALSE, test.lin=TRUE, B=1000)
}
\arguments{
  \item{Y}{a \eqn{N\times T} matrix of the dependent variable where the rows represent the statistical units and the columns the time.}
  \item{TV}{a \eqn{N\times T} transition variable matrix where the rows represent the statistical units and the columns the time. If not specified, the first lag of \eqn{Y} is taken.}
  \item{Xendo}{a \eqn{Nk_{endo}\times T} matrix of the independent variables (if any). The rows contain in order \eqn{X_{1i=1}, X_{1i=2}, \dots, X_{1i=N}, X_{2i=1},X_{2i=2}, \dots, X_{2i=N}, \dots, X_{k_{endo}i=1},} \eqn{X_{k_{endo}i=2}, \dots, X_{k_{endo}i=N}}.}
  \item{Xexo}{a \eqn{Nk_{exo}\times T} matrix of the independent variables (if any). The rows contain in order \eqn{X_{1i=1}, X_{1i=2}, \dots, X_{1i=N}, X_{2i=1},X_{2i=2}, \dots, X_{2i=N}, \dots, X_{k_{exo}i=1}, } \eqn{X_{k_{exo}i=2}, \dots, X_{k_{exo}i=N}}.}
  \item{IV}{a matrix of the additional instrumental variables to include in the analysis. The rows contain in order \eqn{z_{1i=1}, z_{1i=2}, \dots, z_{1i=N}, z_{2i=1}, z_{2i=2}, \dots, z_{2i=N}, \dots, \quad} \eqn{z_{li=1}, z_{li=2}, \dots, z_{li=N}}. The lagged dependent and independent variables are automatically included in the estimation and should not be specified in IV. If not specified, only the lags of the dependent and independent variables are taken as instrumental variables.}
  \item{trimrate}{the trim rate when constructing the grid for estimating the threshold. Default is set to 0.4.}
  \item{ngrid}{the number of grid points to estimate the threshold. Default is set to 100.}
  \item{h0}{the parameter for Silverman's rule of thumb for kernel estimation. The default is set to 1.5.}
  \item{Iweight}{the 1-st step weight matrix. If \code{TRUE}, the identity matrix is used. If \code{FALSE}, the 1-st step weight matrix is constructed from the instrumental variables and the moments matrix. Default is set to \code{FALSE}.}
  \item{test.lin}{if \code{TRUE} (default), the linearity test is performed; if \code{FALSE}, not.}
  \item{B}{if \code{test.lin=TRUE}, it is the number of bootstrap iterations for the linearity test.}
}
\details{
\code{ptm2} performs the first-difference 2-step GMM estimator by Seo and Shin (2016) for panel threshold regression. The following steps are performed in the estimation algorithm:
\enumerate{
\item Computes the first differences of the dependent and independent variables.
\item Defines the instrumental variables matrix, which includes the lagged dependent and independent variables, and the additional instrumental variables in the IV matrix (if specified).
\item Constructs the 1-st step weight matrix according to the Iweight option.
\item Identifies a set of equally distant points between quantiles (\code{trimrate}/2) and 1-(\code{trimrate}/2) of the threshold variable vector. The number of points is specified in the option ngrid.
\item Computes the GMM estimator setting the threshold equal to each selected point and computes the residuals. Identifies the 1st-step estimator by minimizing the residual.
\item Constructs the 2nd-step weight matrix from the residuals associated with the 1st-step estimator.
\item Repeats steps 4. and 5. and returns the 2nd-step estimator of the 1st-differenced parameters, the threshold, and the covariance matrix.
}
}
\value{An object of S4 class "ptm2", which is a list with the following elements:
	\item{threshold}{the estimated threshold and the associated p-value.}
	\item{estimates}{the \code{matrix} of the estimated first-differenced parameters (i.e. parameter(s) estimated for the lower regime, the difference between the parameter(s) estimated for the upper and the lower regime) and the associated p-values.}
	\item{cov.}{the sample version of covariance matrix of the estimated paramaters.}
	\item{residuals.}{the \code{matrix} of the first-differenced residuals.}
	\item{test.lin.}{the \code{list} of the linearity test results containing: the name of the test; the number of bootstrap iterations carried out; the p-value of the test.}
}
\references{
Seo, M.H. and Shin, Y. (2016) Dynamic panels with threshold effect and endogeneity. Journal of Econometrics, 195(2), p.169-186.

Seo, M.H., Kim, S. and Kim Y.-J. (2019) Estimation of dynamic panel threshold model using Stata. The Stata Journal, 19(3), p.685-697.
}
\note{
The method requires at least 6 times of observation: four lags of the dependent and independent variables are used as instruments, and two more are necessary to identify the regime switch (i.e., one per regime).
} 
\author{
Francesca Marta Lilja Di Lascio <marta.dilascio@unibz.it>

Selene Perazzini <selene.perazzini@alumni.imtlucca.it>}
\seealso{
	See also \code{\link{ptm3}} and \code{\link{simptm}}.
}
\examples{
## Select data
data(banana)
data <- banana[which(banana$j==1),]

# Prepare input matrices for ptm2
bioimpedance <-  matrix(data$bioimpedance, ncol=length(unique(data$t)), 
                        byrow=TRUE)

weight <-  matrix(data$weight, ncol=length(unique(data$t)), byrow=TRUE)

## Estimation

# Model on bioimpedance y_{it} with transition variable y_{i(t-1)}
# y_{it}= (phi1_c)1(y_{it-1}<= gamma) + (phi2_c)1(y_{it-1}>gamma)
# Estimation and linear test
ptm2(Y=bioimpedance, TV=NULL, Xendo=NULL, Xexo=NULL, IV=NULL,  
     trimrate=0.4, ngrid=100, h0=1.5, Iweight=FALSE, test.lin=TRUE, 
     B=100)

## NOT RUN
## Model on bioimpedance y_{it} with transition variable y_{i(t-1)} and 
## weight as regressor x_{it}: 
## y_{it}= (phi1_c + phi1_X*x_{it})1(y_{it-1}<= gamma) + 
##         (phi2_c + phi2_X*x_{it})1(y_{it-1}>gamma)
#ptm2(Y=bioimpedance, TV=NULL, Xendo=weight, Xexo=NULL, IV=NULL,  
#     trimrate=0.4, ngrid=100, h0=1.5, Iweight=FALSE, test.lin=FALSE)
##
}
