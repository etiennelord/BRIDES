\encoding{ISO-8859-2}
\name{sample_network}
\alias{sample_network}
\title{compares two networks using a path analysis of total pathways}
\description{This function allows the user to compare two related networks (undirected graphs) by computing the number of paths of each category for a sample of nodes from the original graph.}
\usage{
	sample_network(g1,g2,
					   size, 
					   taxnames,
					   maxdistance=0, 
					   maxtime=3600,
					   maxnode=0, 
					   verbose=FALSE, 
					   file='log.txt', 
					   maxcores=1, 
					   node1, 
					   node2,
					   sample_paths=c(),
					   old_path=c()
					) 
}
\arguments{
  \item{g1}{the original network X}
  \item{g2}{the augmented network Y with additional nodes (all the original nodes from X must be present in the augmented network Y)}
  \item{taxnames}{the taxon name of the nodes added to the original graph (default: we select all nodes that are not in g1)}
  \item{size}{the number of pathways to sample or the percent of node to sample if smaller than 1.0 (default=10)} 	  
  \item{maxnode}{the maximum number of augmented nodes in network Y to take into account (default=0, no maximum number of augmented nodes to take into account. The augmented nodes are sorted by distance to the investigated node pairs by the algorithm.)} 	  
  \item{maxtime}{the maximum search time per pathway (default=3600 seconds)} 	
  \item{maxdistance}{the maximum search distance to the added nodes in network Y (default=0, no maximum distance for augmented nodes to take into account)} 	
  \item{verbose}{flag to save into a file additionnal informations regarding the obtained paths (default=FALSE, the file name should be indicated)}
  \item{file}{filename to save the additionnal informations}
  \item{node1}{if node1 is set,we look only for paths starting from node1}
  \item{node2}{if node2 is set, we will only look at the paths starting at node1 and ending at node2}
  \item{maxcores}{maximum number of cores to use (default=1, use 0 to use half of the total cores available)}
  \item{sample_paths}{vector containing a list of node numbers to sample. Note: this override the size argument}
   \item{old_path}{vector containing a list of node numbers already sampled}
}

\value{This function returns the number of paths of each category: Shortcuts, Detours, Dead ends, Equal paths or disconnected nodes.\cr If not all augmented nodes are visited because of the maxnode, maxtime or maxdistance limits, a path can also be classified as 'Detour or Dead End'.}


%\author{Etienne Lord, Margaux Le Cam, Eric Bapteste, Vladimir Makarenkov and Francois-Joseph Lapointe}

\examples{
	## Sample 10 of the 55 pathways in Sample_1
	data(Sample_1)   
	sample_network(g1,g2, size=10)
	## Repeated sampling (5) of 10 pathways 
	\dontrun{
	old_path <- c()
	for (i in 1:5) {
       sample_paths <- sample_path(10, length(V(g1)), old_path=old_path);
       sample_network(g1,g2, sample_paths=sample_paths);
       old_path <- c(old_path,sample_paths); 
	 }
	}
}

%\keyword{ }

