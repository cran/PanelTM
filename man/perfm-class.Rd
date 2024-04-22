\name{perfm-class}
\docType{class}
\alias{perfm-class}
\title{Class "perfm"}
\description{A class for \code{perfm} and its extensions}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("perfm", ...)}.
}
\section{Slots}{
     \describe{
	\item{\code{trueval}}{Object of class \code{numeric}. The true value of the parameter.}
	\item{\code{rb}}{Object of class \code{numeric}. The relative bias.}
	\item{\code{nrb}}{Object of class \code{numeric}. The bias.}
	\item{\code{rrmse}}{Object of class \code{numeric}. The relative root mean squared error.}
	\item{\code{rmse}}{Object of class \code{numeric}. The root mean squared error.}
 }
}
\section{Methods}{
No methods defined with class "perfm" in the signature.
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
    See also \code{\link{perfm}}.
}
\examples{
    showClass("perfm")
}
