\encoding{ISO-8859-2} 
\name{networkX}
\alias{networkX}
\docType{data}
\title{An original network X with 12 nodes and 11 edges}
\description{This file contains one sample network X as described in: \cr 
			 Lord et al. (2016). BRIDES: a new fast algorithm and software for characterizing evolving similarity networks using \cr
			 breakthroughs, roadblocks, impasses, detours, equals and shortcuts. PloS one, 11(8), e0161474.}
\usage{data(networkX)}
\format{
  An igraph object containing network X
  \describe{
    \item{networkX}{An original \emph{igraph} network}
  }
}

\examples{
\dontrun{
	data(networkX)
	data(networkY)
	plot_network(networkX,networkY)
	BRIDES(networkX,networkY)
	
}
#Prunning unconnected nodes...
#==========================================================
#Total 66 pathways divided into 1 groups.
#==========================================================
#Run parameters:
#Networks                    : undirected
#Nodes in networkX           : 12 
#Nodes in networkY           : 15 
#Total added nodes (K)       : 3 
#Attributes for added nodes  : default 
#Total paths                 : 66 
#Weighted edges              : TRUE 
#Path information (trace)    : complete
#Group size                  : 1000 
#Start group                 : 1 
#End group                   : 1 
#Maxdistance                 : 100 
#Maxnode                     : 100 
#Maxcores                    : 1 
#==========================================================
#  	B	R	I	D	E	S	(utime	stime)	
#	9	2	20	26	5	4	0.07	1.27	
#====================== RESULTS ===========================
#     B      R      I      D      E      S (utime stime) 
#  9.00   2.00  20.00  26.00   5.00   4.00   0.07   1.27 
}
\keyword{dataset}
