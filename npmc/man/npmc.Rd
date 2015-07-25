\name{npmc}
\alias{npmc}
\alias{print.npmc}
\alias{summary.npmc}

\title{ Nonparametric Multiple Comparisons }

\description{
  npmc implements the nonparametrical multiple testprocedures
  (Behrens-Fisher- and Steel-type for the all-pairs and many-to-one situations)
  described in the paper
  "A unified approach to Simultaneous Rank Test Procedures in the 
  Unbalanced One-way Layout" 
  written by Ullrich Munzel and Ludwig Hothorn.
}

\usage{
npmc(dataset, control=NULL, df=2, alpha=0.05)

\method{summary}{npmc}(object, type="both", info=TRUE, short=TRUE, corr=FALSE, ...)
%%summary.npmc(object, type="both", info=TRUE, short=TRUE, corr=FALSE, ...)
}
\arguments{
  \item{dataset}{
    a data-frame with variables 
    'var' (the response-variable)and 
    'class' (containing the class-level).
    The data may be unsorted and unbalanced.
    The 'summary'-function supports 'name' and 'description' attributes
  }
  \item{control}{
    level of control-group for the many-to-one situation
    or NULL for the all-pairs situation (default)
  }
  \item{df}{
    determines the BF-teststatistics' asymptotic distribution function
    (0 = standard-normal approximation,
    1 = simple t-approximation,
    2 = Satterthwaite t-approximation (default))
  }
  \item{alpha}{
    Level for the (1-alpha) confidence-intervals
  }

  The following arguments are used by the 'summary'-function:
  \item{object}{
    An object of type 'npmc'
  }
  \item{type}{
    You may select either "BF" or "Steel". All other values extract
    information for both types
  }
  \item{info}{
    Prints info about the procedure and the data
  }
  \item{short}{
    Only prints the most relevant items of the test-results
  }
  \item{corr}{
    Prints the correlation-matrices
  }
  \item{\dots}{further arguments to be passed to or from methods.}
  
}
\details{
  npmc performs nonparametrical multiple testprocedures
  (Behrens-Fisher- and Steel-type for the all-pairs and many-to-one
  situations)
  and computes the simultaneous (1-alpha) confidence limits for the
  relative effects.

  The one-sided tests reject if group with smaller index has larger values
  due to the calculation of the relative effect-estimators. 
}
\value{
  A list of several other values and structures with subvalues
  \item{info}{
    group.index: an integer number to identify the group
    
    class.level: the class-level of this group
    
    nobs: the number of observations in this group
  }
  \item{corr}{
     BF: The Behrens-Fisher-type correlation-matrix
     
     Steel: The Steel-type correlation-matrix

     The matrices have attributes 'adjusted' which indicate if 
     negative eigenvalues were changed to zero (see 'Notes'-section).
   }
   \item{test}{
     BF / Steel: Two lists containing the test-results and some other
     characteristics of the multiple nonparametric Behrens-Fisher- and 
     Steel-type testprocedures
     
       cmp: names the compared groups ('a-b')

       gn: the sum of both sample-sizes

       effect: the relative effect-estimator

       variance: the variance-estimator

       std: the standard-deviation

       statistic: the test-statistic

       p-value 1s: the 1-sided p-value

       p-value 2s: the 1-sided p-value

       zero: TRUE if zero variances occured and were substituted by 0.0000001
   }
   \item{control}{
     The factor-level of the control group (NULL for
     allpairs-comparison)
   }
   \item{df.method}{
     The original df-paramter
   } 
   \item{df}{
     The estimated degrees of freedom for the multivariate t-distribution
     or 0 if the standard normal distribution was used for approximating the
     teststatistics' distribution (depending on the df-parameter).
   }
 }
 
\references{ "A unified approach to Simultaneous Rank Test Procedures in
  the Unbalanced One-way Layout"
  written by Ullrich Munzel and Ludwig Hothorn.
}

\author{ Joerg Helms }

\note{
  This function requires the 'mvtnorm' package to calculate the p-values
  for the test-statistics. If this package is not available on your system,
  the result will contain NA's as p-values.

  The functions in the 'mvtnorm' package seem to use randomized values for 
  integral calculations, so the results from npmc concerning p-values
  and confidence-intervals differ from call to call and can be
  recognized only as an approximative solution.
  
  If a correlation-matrix has negative eigenvalues, they are replaced by 0
  as recommended by the paper. The matrix's 'adjusted' attribute indicates 
  whether or not the matrix was forced to be positive semidefinite.
}

%%\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }

\examples{
data(brain)
summary(npmc(brain), type="BF")
}
\keyword{ nonparametric }
\concept{ multiple comparisons }
