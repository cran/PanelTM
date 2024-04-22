\name{cpoint-class}
\docType{class}
\alias{cpoint-class}
\title{Class "cpoint"}
\description{A class for \code{cpoint} and its extensions}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cpoint", ...)}.
}
\section{Slots}{
     \describe{
	\item{\code{CP}:}{Object of class \code{matrix}. A matrix reporting the computed change point per each statistical unit and each level/category/value of the third way.}
 }
}
\section{Methods}{
No methods defined with class "cpoint" in the signature.
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
    See also \code{\link{ptm2}} and \code{\link{ptm3}}.
}
\examples{
showClass("cpoint")
}
