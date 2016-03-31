\encoding{ISO-8859-2} 
\name{Sample_2}
\alias{Sample_2}
\docType{data}
\title{An network X of 6 nodes and an augmented network Y with 7 nodes}
\description{This file contains two artificial networks used to illustrate the application of the SDDE package.}
\usage{data(Sample_2)}
\format{
  A data frame with 6 nodes for network g1 and 7 nodes for network g2:
  \describe{
    \item{g1}{An original \emph{igraph} network}
    \item{g2}{An augmented \emph{igraph} network}
  }
}

\examples{
\dontrun{
	data(Sample_2)
	complete_network(g1,g2)
}
##          disconnected nodes shortcuts equals detours dead ends
##                           0         0      1       8         6
##          Dead ends or detour total user time system time real time
##                            0    15     0.656       0.049     1.023
# Total new nodes in g2:	1	
# Number of edges in g2:	7	
# Number of nodes in g2:	7	
# Number of nodes in g1:	6	
# ====================================
# x2	x1	Dead
# x3	x2	Dead
# x4	x2	Dead
# x5	x2	Detour
# x3	x1	Dead
# x5	x4	Detour
# x6	x2	Detour
# x4	x1	Dead
# x6	x4	Detour
# x4	x3	Dead
# x5	x1	Detour
# x5	x3	Detour
# x6	x1	Detour
# x6	x3	Detour
# x6	x5	Equal
#====================================        
}

\keyword{datasets}
