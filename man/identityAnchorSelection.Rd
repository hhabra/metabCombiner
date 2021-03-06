% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectAnchors.R
\name{identityAnchorSelection}
\alias{identityAnchorSelection}
\title{Select Matching Ids as Anchors}
\usage{
identityAnchorSelection(cTable, windx, windy, useID, brackets)
}
\arguments{
\item{cTable}{data frame, contains only feature ids, mzs, rts, Qs, & labels}

\item{windx}{numeric positive retention time exclusion window in X dataset}

\item{windy}{numeric positive retention time exclusion window in Y dataset}

\item{useID}{logical. Operation proceeds if TRUE, terminates otherwise.}

\item{brackets}{If useID = TRUE, bracketed identity strings of the types
included in this argument will be ignored}
}
\value{
combinedTable with updated anchor labels
}
\description{
This is an optional helper function for \code{selectAnchors}. Uses identities
to guide selection of ordered retention time pairs. If useID option is set to
TRUE, it will select pairs of features with matching ID character strings
before proceeding with iterative anchor selection.
}
\details{
Identity anchors are allowed to violate constraints of m/z, Q, and
rtq difference tolerances, and will not be removed if they fall within a rt
exclusion window of other features. If a name appears more than once, only
the pair with the highest relative abundance is selected.
}
\seealso{
\code{\link{selectAnchors}}
}
