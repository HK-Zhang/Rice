\name{kronen}
\alias{kronen}
\title{ Baumkronen }
\usage{data(kronen)}
\description{
The \code{kronen} data frame has 72 rows and 2 columns.
It contains pairs of observations (var) and factor levels (class).
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{class}{ - factor with 3 levels}
    \item{var}{ - a numeric vector, the response variable}
  }
}
\details{
An example for an unbalanced one-factorial design.
}
%\source{}
%\references{}
\examples{
data(kronen)
npmc(kronen)
}
\keyword{datasets}
