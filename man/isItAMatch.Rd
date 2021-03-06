% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NeedlemanWunsch.R
\name{isItAMatch}
\alias{isItAMatch}
\title{Checks if the current characters are matching}
\usage{
isItAMatch(sequenceA, sequenceB, i, j)
}
\arguments{
\item{sequenceA}{The first nucleotide sequence (positioned on the rows)}

\item{sequenceB}{The second nucleotide sequence (positioned on the columns)}

\item{i}{The index for sequenceA (rows)}

\item{j}{The index for sequenceB (columns)}
}
\value{
TRUE if the two characters are matching, FALSE otherwise
}
\description{
isItAMatch is a utility function that checks whether the i-th and j-th characters of the two sequences are matching
}
\examples{
isItAMatch("GTT","GCATT",1,2)

}
