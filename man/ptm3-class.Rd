\name{ptm3-class}
\docType{class}
\alias{ptm3-class}

\title{Class "ptm3"}
\description{A class for \code{ptm3} and its extensions}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ptm3", ...)}.
}
\section{Slots}{
     \describe{
	\item{\code{threshold}:}{Object of class \code{numeric}. A vector of the estimated thresholds and the associated p-values.}
	\item{\code{param}:}{Object of class \code{matrix}. A matrix of the estimated first-differenced parameters of the model and the associated p-values.}
	\item{\code{cov.}:}{Object of class \code{array}. The sample covariance matrix of the estimated parameters per each value of the third way.}
	\item{\code{residuals.}:}{Object of class \code{list}. The estimated first-differenced residuals per each value of the third way.}
	\item{\code{test.lin.}:}{Object of class \code{list}. Result of the linearity test.}
 }
}
\section{Methods}{
No methods defined with class "ptm3" in the signature.
}
\references{
Di Lascio, F.M.L. and Perazzini, S. (202x) A three-way dynamic panel threshold regression model for change point detection in bioimpedance data. WP BEMPS <https://repec.unibz.it/bemps104.pdf>.

Di Lascio, F.M.L. and Perazzini, S. (2022) Change point detection in fruit bioimpedance using a three-way panel model. Book of short papers - SIS2022, p.1184-1189, Eds. A. Balzanella, M. Bini, C. Cavicchia, R. Verde. Pearson. ISBN: 978-88-9193-231-0.
}
\author{
Selene Perazzini <selene.perazzini@alumni.imtlucca.it>

Francesca Marta Lilja Di Lascio <marta.dilascio@unibz.it>
}
\seealso{
    See also \code{\link{ptm3}}.
}
\examples{
showClass("ptm3")
}
\keyword{classes}
