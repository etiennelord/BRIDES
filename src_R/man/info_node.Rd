\encoding{ISO-8859-2}
\name{info_node}
\alias{info_node}
\title{returns additional information regarding the nodes of networks X and Y (original and augmented)}
\description{This function returns some additional information regarding the nodes of networks X and Y such as the reachability of nodes in networks X to the augmented nodes of network Y.}
\usage{info_node(g1,g2,taxnames, maxnode, maxdistance)}
\arguments{
  \item{g1}{the original network X}
  \item{g2}{the augmented network Y with additional nodes (all the original nodes from X must be present in the augmented network Y)}
  \item{taxnames}{the taxon name of the nodes added to the original graph (default: we select all nodes that are not in g1)}
  \item{maxnode}{the maximum number of augmented nodes in network Y to take into account (default=0, no maximum number of augmented nodes to take into account. The augmented nodes are sorted by distance to the investigated node pairs by the algorithm.)} 	  
  \item{maxdistance}{the maximum search distance to the added nodes in network Y (default=0, no maximum distance for augmented nodes to take into account)} 	
}

\value{This function return information about the reachability of each node such as : \emph{reach_count}(number of reachable augmented nodes), \emph{reach_max_dist} (maximal distance to an augmented node), \emph{reach_min_dist} (minimal distance to an augmented node), \emph{reach_mean_dist} (mean distance to augmented nodes)
}

%\author{Etienne Lord, Margaux Le Cam, Eric Bapteste, Vladimir Makarenkov  and Francois-Joseph Lapointe}

\examples{
   ## Report reachability informations about the node in sample1
   \dontrun{
   data(Sample_1)
   info_node(g1,g2)
   }
   ## Results:
   ## Total of pathways to investigate: 55 
   ## Number of edges in network Y: 17 
   ## Number of nodes in network Y: 14 
   ## Number of nodes in network X: 11 
   ## Nodes distribution in network Y (first row taxa, second row count):
   ## 1  2 
   ## 11  3 
   ## $reach_count
   ## [1] 3 3 3 3 3 3 0 3 3 0 3 3 3 3 3 0 3 3 0 3 3 3 3 0 3 3 0 0 3 3 0 3 3 0 3 3 0 3
   ## [39] 3 0 3 0 3 3 0 0 3 3 0 0 0 0 3 0 0
   ##  $reach_max_dist
   ## [1] 3 3 3 3 3 1 0 3 1 0 3 3 3 1 1 0 3 1 0 3 3 3 3 0 3 3 0 0 3 1 0 3 1 0 3 1 0 3
   ## [39] 1 0 1 0 3 1 0 0 1 1 0 0 0 0 1 0 0
   ## $reach_min_dist
   ## [1] 2 1 2 2 2 2 0 2 2 0 1 1 1 2 3 0 2 3 0 1 1 1 2 0 1 2 0 0 2 2 0 2 2 0 2 2 0 2
   ## [39] 2 0 2 0 2 2 0 0 2 2 0 0 0 0 2 0 0
   ## $reach_mean_dist
   ##  [1] 2 2 2 2 2 2 0 2 2 0 2 2 2 2 2 0 2 2 0 2 2 2 2 0 2 2 0 0 2 2 0 2 2 0 2 2 0 2
   ## [39] 2 0 2 0 2 2 0 0 2 2 0 0 0 0 2 0 0   
}

%\keyword{reachability, networks}

