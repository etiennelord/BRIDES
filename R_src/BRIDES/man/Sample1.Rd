\encoding{ISO-8859-2} 
\name{Sample1}
\alias{Sample1}
\docType{data}
\title{An original and an augmented sample networks}
\description{These networks contains 12 nodes for network X and 15 nodes for network Y to illustrate the application of the BRIDES package.}
\usage{data(Sample1)}
\format{
  A data frame with the following 2 variables.
  \describe{
    \item{g1}{An original \emph{igraph} network (network X)}
    \item{g2}{An augmented \emph{igraph} network (network Y)}
  }
}

\examples{
\dontrun{
	data(Sample1)
	BRIDES(g1,g2)
}
## Expected results
##
##==========================================================
##Total 66 pathways to evaluated divided into 1 groups.
##==========================================================
##Run parameters:
##Taxnames: default
##Randomize paths: FALSE
##Total nodes: 12
##Group size: 1000
##Start group: 1
##End group: 1
##Maxdistance: 100
##Maxtime: 100
##Maxnode: 100
##Maxcores: 1
##==========================================================
##       B       R       I       D       E       S       (utime  stime)
##       9       2       20      26      5       4       0.23    12.78
##====================== RESULTS ===========================
##     B      R      I      D      E      S (utime stime)
##  9.00   2.00  20.00  26.00   5.00   4.00   0.23  12.78
}
\keyword{datasets}
