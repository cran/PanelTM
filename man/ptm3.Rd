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
\details{
\code{ptm3} performs the first-difference 2-step GMM estimation of the three-way panel threshold regression model in Di Lascio and Perazzini (202X).
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
Di Lascio, F.M.L. and Perazzini, S. (202x) A three-way dynamic panel threshold regression model for change point detection in bioimpedance data. WP BEMPS <https://repec.unibz.it/bemps104.pdf>.

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

ptm3(data.=banana, nameI="i", nameT="t", nameJ="j",nameY="bioimpedance", 
    nameTV=NULL, nameXendo="weight", nameXexo=NULL, nameIV=NULL,
    trimrate=0.4, ngrid=100, h0=1.5, Iweight=FALSE, test.lin=FALSE)
}



