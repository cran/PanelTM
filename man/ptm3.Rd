\name{ptm3}
\alias{ptm3}
\title{Three-way panel threshold regression model.}
\description{
Three-way panel threshold regression model and its estimation through the 2-step GMM estimator.
}
\usage{
ptm3(data., nameI, nameT, nameJ, nameY, nameTV=NULL, nameXendo=NULL, 
    nameXexo=NULL, nameIV=NULL, trimrate=0.4, ngrid=100, h0=1.5,
    Iweight=FALSE, test.lin=TRUE, B=1000)
}
\arguments{
  \item{data.}{a \code{data.frame} with the following variables (by cols): \eqn{I} (statistical units ID), \eqn{T} (time), \eqn{J} (third-level), \eqn{Y} (dependent variable), and if not null, \eqn{TV} (threshold variable), \eqn{Xendo} (endogenous regressor(s)), \eqn{Xexo} (exogenous regressor(s)), \eqn{IV} (instruments).}
  \item{nameI}{the name of the (numerical) variable that identifies the statistical units.}
  \item{nameT}{the name of the (numerical) time variable.}
  \item{nameJ}{the name of the (numerical) variable that indicates the third dimension.}
  \item{nameY}{the name of the (numerical) dependent variable.}
  \item{nameTV}{the name of the (numerical) transition variable. If not specified, the first lag of \eqn{Y} is taken.}
  \item{nameXendo}{the names of the (numerical) independent endogenous variables (if any).}
  \item{nameXexo}{the names of the (numerical) independent exogenous variables (if any).}
  \item{nameIV}{the names of the (numerical) instrumental variables (if any).}
  \item{trimrate}{the trim rate when constructing the grid for estimating the threshold. The default value is set to 0.4.}
  \item{ngrid}{the number of grid points to estimate the threshold. The default is set to 100.}
  \item{h0}{the parameter for Silverman's rule of thumb for kernel estimation. The default is set to 1.5.}
  \item{Iweight}{the 1-st step weight matrix. If \code{TRUE}, the identity matrix is used. If \code{FALSE}, the 1-st step weight matrix is constructed from the instrumental variables. The default is set to \code{FALSE}.}
  \item{test.lin}{if \code{TRUE} (default), the linearity test is performed; if \code{FALSE}, not.}
  \item{B}{if \code{test.lin=TRUE}, the number of bootstrap iterations for the linearity test.} 
}

\details{\code{ptm3} estimates a three-way dynamic panel threshold regression model using the first-difference instrumental-variable generalized method of moments (FD-IV-GMM) estimator in Di Lascio and Perazzini (2026).

Let \eqn{i = 1,\ldots,n} index individuals, \eqn{t = 1,\ldots,T} index time periods, and \eqn{j = 1,\ldots,J} index levels of the third dimension. The model is \deqn{y_{ijt} = (1,x_{ijt}')\phi_{1j} I(q_{ijt} \leq \gamma_j) + (1,x_{ijt}')\phi_{2j} I(q_{ijt} > \gamma_j) + \varepsilon_{ijt}, } where \eqn{y_{ijt}} is the dependent variable, \eqn{x_{ijt}} is a vector of time-varying regressors, \eqn{q_{ijt}} is the threshold variable, \eqn{\gamma_{j}} is the threshold parameter specific to level \eqn{j}, and \eqn{I(\cdot)} is the indicator function. The parameter vectors \eqn{\phi_{1j}} and \eqn{\phi_{2j}} contain the intercept and slope coefficients in the lower and upper regimes, respectively. The error term is \deqn{ \varepsilon_{ijt} = \mu_{i} + \lambda_{j} + \nu_{ijt}, } where \eqn{\mu_{i}} is an unobserved individual-specific effect,  \eqn{\lambda_{j}} is an unobserved effect associated with the third dimension, and \eqn{\nu_{ijt}} is an idiosyncratic error.

The model is estimated separately for each level \eqn{j}, allowing the regression coefficients and threshold parameter to vary across levels of the third dimension. The threshold variable can be specified by the user. If it is not provided, the first lag of the dependent variable is used. Regressors can be specified as exogenous or endogenous. If no regressors are supplied, the model includes an intercept only. Lagged values of the dependent variable and regressors are also used as instruments and additional instruments can be supplied through the \code{nameIV} argument.

Once fixed \eqn{j}, the steps to estimate the model for the specific level \eqn{j} coincide with the ones described in the Details section of the \code{\link{ptm2}} documentation.

The regressors in the output are always ordered with the exogenous variables first, followed by the endogenous variables. Within each group, the variables retain the same order as specified by the user in the input arguments. Therefore, the estimated coefficients and their associated p-values should be interpreted according to this ordering.

When \code{test.lin = TRUE}, the function also performs a bootstrap-based test for the presence of a threshold effect separately for each level \eqn{j}. The null hypothesis is that the regime-specific coefficient difference is zero, \eqn{H0: \delta_j = 0} against \eqn{H1: \delta_j \neq 0}. The test uses a supremum Wald statistic over the threshold grid, with the number of bootstrap replications controlled by \code{B}.

If \eqn{J = 1}, the function estimates a two-way panel threshold model as implemented in \code{\link{ptm2}}.
}

\value{An object of S4 class "ptm3", which is a list with the following elements:
	\item{threshold}{the \code{matrix} of the estimated thresholds and the associated p-values.}
	\item{param}{the \code{matrix} of the estimated first-differenced parameters of the model and the associated p-values.}
	\item{cov.}{the \code{array} containing the estimated covariance matrix per each value of the third way.}
	\item{residuals.}{the \code{list} containing the estimated first-differenced residuals per each value of the third way.}
	\item{test.lin.}{the \code{list} of the linearity test results containing: the name of the test; the number of bootstrap iterations carried out; the p-value of the test 
per each value of the third way.}
}

\references{
Di Lascio, F.M.L. and Perazzini, S. (2026) A three-way dynamic panel threshold regression model for change point detection in bioimpedance data. Stat Methods Appl 35, 849-867.

Di Lascio, F.M.L. and Perazzini, S. (2022) Change point detection in fruit bioimpedance using a three-way panel model. Book of short papers - SIS2022, p.1184-1189, Eds. A. Balzanella, M. Bini, C. Cavicchia, R. Verde. Pearson. ISBN: 978-88-9193-231-0.
}

\note{The estimation method requires at least \eqn{t\geq6}: four lags of the dependent and independent variables are used as instruments, and two more are necessary to identify the regime switch (i.e., one per regime). Note that the instant times in the nameT variable are automatically ordered increasingly.
}
\author{
Francesca Marta Lilja Di Lascio <marta.dilascio@unibz.it>

Selene Perazzini <selene.perazzini@alumni.imtlucca.it>
}
\seealso{
	See also \code{\link{ptm2}} and \code{\link{simptm}}.
}
\examples{
# Import data
data(banana)

## NOT RUN
## Model on bioimpedance y_{ijt} with transition variable y_{ij(t-1)}
## y_{ijt} = (phi1_{jc})1(y_{ijt-1}<=gamma_j) + 
##           (phi2_{jc})1(y_{ijt-1}>gamma_j)
#
#ptm3(data.=banana, nameI="i", nameT="t", nameJ="j", 
#    nameY="bioimpedance", nameTV=NULL, nameXendo=NULL, nameXexo=NULL, 
#    nameIV=NULL, trimrate=0.4, ngrid=100, h0=1.4, Iweight=FALSE, 
#    test.lin=FALSE)
##

# Model on bioimpedance y_{ijt} with transition variable y_{ij(t-1)} 
# and time-varying regressor x_{ijt}: 
# y_{it}=(phi1_{jc}+phi1_{jX}*x_{ijt})1(y_{ijt-1}<=gamma_j) + 
# (phi2_{jc}+phi2_{jX}*x_{ijt})1(y_{ijt-1}>gamma_j)

p3 <- ptm3(data.=banana, nameI="i", nameT="t", nameJ="j",nameY="bioimpedance", 
    nameTV=NULL, nameXendo="weight", nameXexo=NULL, nameIV=NULL,
    trimrate=0.4, ngrid=100, h0=1.5, Iweight=FALSE, test.lin=FALSE)
    
# Show main output 
show(p3)

}



