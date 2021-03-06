% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/methods-metabCombiner.R
\name{nonmatched}
\alias{nonmatched}
\alias{nonmatched,metabCombiner-method}
\title{Get Nonmatched Features}
\usage{
nonmatched(object, data = "x")

\S4method{nonmatched}{metabCombiner}(object, data = "x")
}
\arguments{
\item{object}{metabCombiner object}

\item{data}{dataset identifier for \code{metabCombiner} objects; if NULL,
returns full list of non-matched features}
}
\value{
Data frame of non-matched features corresponding to data argument
}
\description{
Features that lack a any counterparts in the complementary dataset may be
obtained from this method. If data is set to "x" or "y", will retrieve data
from the current X or Y dataset, respectively. If data is set to NULL, will
retrieve the list of nonmatched features.
}
\examples{
data(plasma30)
data(plasma20)

p30 <- metabData(head(plasma30,500), samples = "CHEAR")
p20 <- metabData(head(plasma20,500), samples = "Red", rtmax = 17.25)
p.comb <- metabCombiner(xdata = p30, ydata = p20, binGap = 0.005)

nnmx <- nonmatched(p.comb, data = "x")
nnmy <- nonmatched(p.comb, data = "y")

}
