\name{brain}
\alias{brain}
\title{ reaction time of brain-damaged patients }
\usage{data(brain)}
\description{
  The \code{brain} data frame has 90 rows and 2 columns.
  It contains pairs of observations (var) and factor levels (class).
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{class}{ - factor with levels
      \code{c} = control-group,
      \code{l} = left hemisphere brain-damage, 
      \code{r} = right hemisphere brain-damage
    }
    \item{var}{ - score for the reaction}
  }
}
\details{
The effects of uniliteral left hemisphere (LH) and right hemisphere (RH)
lesions on the accuracy of choice and speed of response in a four-choice
reaction time task were examined by Tartaglione et al. (1191). The
subjects consisted of 30 controls, 30 LH and 30 RH brain-damaged
patients, so this is an example for a balanced one-factorial design.
}
%\source{}
\references{
"A unified approach to Simultaneous Rank Test Procedures in
  the Unbalanced One-way Layout"
  written by Ullrich Munzel and Ludwig Hothorn.
}
\examples{
data(brain)
npmc(brain)
}
\keyword{datasets}
