\name{ptm2-class}
\docType{class}
\alias{ptm2-class}
\title{Class "ptm2"}
\description{A class for \code{ptm2} and its extensions}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ptm2", ...)}.
}
\section{Slots}{
     \describe{
	\item{\code{threshold}:}{Object of class \code{numeric}. The estimated threshold and the associated p-value.}
	\item{\code{estimates}:}{Object of class \code{numeric}. The vector of the parameter(s) estimated for the lower regime, the difference between the parameter(s) estimated for the upper and the lower regime, and the associated p-values.} 
	\item{\code{cov.}:}{Object of class \code{matrix}. The sample version of covariance matrix of the estimated paramaters.}
	\item{\code{residuals.}:}{Object of class \code{matrix}. A matrix of the first-differenced residuals.}
	\item{\code{test.lin.}:}{Object of class \code{list}. Result of the linearity test.}
    }
}
\section{Methods}{
No methods defined with class "ptm2" in the signature.
}
\references{
Seo, M.H. and Shin, Y. (2016) Dynamic panels with threshold effect and endogeneity. Journal of Econometrics, 195(2), p.169-186. 

Seo, M.H., Kim, S. and Kim Y.-J. (2019) Estimation of dynamic panel threshold model using Stata. The Stata Journal, 19(3), p.685-697. 
} 
\author{
Selene Perazzini <selene.perazzini@alumni.imtlucca.it>

Francesca Marta Lilja Di Lascio <marta.dilascio@unibz.it>
} 
\seealso{
    See also \code{\link{ptm2}}.
}
\examples{
showClass("ptm2")
}
