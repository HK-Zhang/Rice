\name{report}
\alias{report}
\title{    Report Printing    }
\description{
  Function to present text on screen.
}
\usage{
report(msg=NULL, style=0, char="-")
}
\arguments{
  \item{msg}{ The text to be printed (any vector-type) }
  \item{style}{ The style in which to print the text. Possible values
    are 'underlined' and 'border' (and valid abbreviations).
    Any other value defaults to 'blank', i.e. the pure text-vector.
  }
  \item{char}{ A character to use for the style. }
}
%\details{}
\value{ A matrix containing the printed lines in one row. }
\author{ Joerg Helms }
\note{ report will always print a blank line first. }
\examples{ 
    report("It's a nice day!", style="border", char="*") 
}
\keyword{ print }
\keyword{ utilities }
