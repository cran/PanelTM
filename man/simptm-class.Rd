\name{simptm-class}
\docType{class}
\alias{simptm-class}
\title{Class "simptm"}
\description{A class for \code{simptm} and its extensions}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("simptm", ...)}.
}
\section{Slots}{
     \describe{
	\item{\code{simulation}}{Object of class \code{"list"}. A list of dataframes of simulated data.
		The list contains \code{B} dataframes, one for each replication. Each dataframe contains columns: \eqn{i}, \eqn{j}, \eqn{t}, \eqn{Y} and possibly \eqn{X}.}
 }
}
\section{Methods}{
No methods defined with class "simptm" in the signature.
}
\references{
Di Lascio, F.M.L. and Perazzini, S. (202x) A three-way dynamic panel threshold regression model for change point detection in bioimpedance data. WP BEMPS <https://repec.unibz.it/bemps104.pdf>.

Di Lascio, F.M.L. and Perazzini, S. (2022) Change point detection in fruit bioimpedance using a three-way panel model. Book of short papers - SIS2022, p.1184-1189, Eds. A. Balzanella, M. Bini, C. Cavicchia, R. Verde. Pearson. ISBN: 978-88-9193-231-0.

Seo, M. and Shin, Y. (2016) Dynamic panels with threshold effect and endogeneity. Journal of Econometrics, 195(2), p.169-186.
}
\author{
Selene Perazzini <selene.perazzini@alumni.imtlucca.it>

Francesca Marta Lilja Di Lascio <marta.dilascio@unibz.it>
}
\seealso{
    See also \code{\link{simptm}}.
}
\examples{
    showClass("simptm")
}
