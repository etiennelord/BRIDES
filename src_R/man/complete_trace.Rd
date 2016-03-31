\encoding{ISO-8859-2}
\name{complete_trace}
\alias{complete_trace}
\title{return properties of a single path in two given networks (original and augmented, presented as undirected graphs) using a path analysis}
\description{This function is equivalent to the \code{\link{complete_network}} function, but it analyses a single path returning information about the: path type, the enumeration of traversed nodes, the traversed node taxa, the original path length and the augmented path length.}
\usage{
	complete_trace(g1,g2,
					  node1="default",
					  node2="default",
					  taxnames,					  
					  maxdistance=0,
					  maxtime=3600,
					  maxnode=0,
					  maxcores=1
					) 
}
\arguments{
  \item{g1}{the original network X}
  \item{g2}{the augmented network Y with additional nodes (all the original nodes from X must be present in the augmented network Y)}
  \item{taxnames}{the taxon name of the nodes added to the original graph (default: we select all nodes that are not in g1)}
   \item{node1}{the \emph{from} node (either the node name or the node number) of the path to investigate}
  \item{node2}{the \emph{to} node (either the node name or the node number) of the path to investigate}
  \item{maxnode}{the maximum number of augmented nodes in network Y to take into account (default=0, no maximum number of augmented nodes to take into account. The augmented nodes are sorted by distance to the investigated node pairs by the algorithm.)} 	  
  \item{maxtime}{the maximum search time per pathway (default=3600 seconds)} 	
  \item{maxdistance}{the maximum search distance to the added nodes in network Y (default=0, no maximum distance for augmented nodes to take into account)} 	
  \item{maxcores}{maximum number of cores to use (default=1, use 0 to use half of the total cores available )}
}

\value{This function returns a list containing the \emph{from} and \emph{to} node names, the \emph{path type} (\emph{Shortcut}, \emph{Detour}, \emph{Dead end} and \emph{Equal} paths or \emph{unreachable} node), the \emph{original path length}, the \emph{augmented path length}, the traversed \emph{path} nodes and the \emph{path_visited_taxa} (taxa of each visited node). \cr
\cr
}

%\author{Etienne Lord, Margaux Le Cam, Eric Bapteste, Vladimir Makarenkov and Francois-Joseph Lapointe}

\examples{
## Searching the sample data for the path type between node x1 and x5
data(Sample_1)
##
\dontrun{complete_trace(g1,g2,node1="x1", node2="x5")}
## Expected results:
## $from
## [1] "x1"
## $to
## [1] "x5"
## $path_type
## [1] "Detour"
## $original_path_length
## [1] 2
## $augmented_path_length
## [1] 3
## $path
## [1] "x1" "x14" "x6" "x5" 
## $path_visited_taxa
## [1] 1 2 1 1
}

%\keyword{trace, path composition}

