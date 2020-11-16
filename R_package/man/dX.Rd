\encoding{ISO-8859-2} 
\name{dX}
\alias{dX}
\docType{data}
\title{An original directed network X with 4 nodes and 2 edges}
\description{This dataset is composed of a directed augmented sample network with one disconnected node}
\usage{data(networkX)}
\format{
  An igraph object containing network X
  \describe{
    \item{dX}{An augmented \emph{igraph} network}
  }
}

\examples{
\dontrun{
	data(dX)
	data(dY)
	plot_network(dX,dY, legend_position="topleft")
	BRIDES(dX,dY)
	# [Results]
	#==========================================================
	#Total 12 pathways divided into 1 groups.
	#==========================================================
	#Run parameters:
	#Networks                    : directed
	#Nodes in networkX           : 4 
	#Nodes in networkY           : 5 
	#Total added nodes (K)       : 1 
	#Attributes for added nodes  : default 
	#Total paths                 : 12 
	#Weighted edges              : FALSE 
	#Path information (trace)    : complete
	#Group size                  : 1000 
	#Start group                 : 1 
	#End group                   : 1 
	#Maxdistance                 : 100 
	#Maxnode                     : 100 
	#Maxcores                    : 1 
	#==========================================================
	#  	B	R	I	D	E	S	(utime	stime)	
	#	5	3	4	0	0	0	0.03	0.31	
	#====================== RESULTS ===========================
	#	 B      R      I      D      E      S (utime stime) 
	#    5.00   3.00   4.00   0.00   0.00   0.00   0.03   0.31 
}
}
\keyword{dataset}
