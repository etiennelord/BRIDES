\encoding{ISO-8859-2}
\name{complete_restart}
\alias{complete_restart}
\title{compares two given networks (original and augmented, presented as undirected graphs) using a path analysis}
\description{This function is equivalent to the \code{\link{complete_network}} function, but it uses a fixed number of pathways to be investigated, which are divided in groups to allow restarting or better dispatching of tasks in a cluster environment.}

\usage{
	complete_restart(g1,g2,
					     taxnames,
					     resultfile='result.txt',
					     size=1000,
					     start=1,
					     end=0,
					     maxdistance=0,
					     maxtime=3600,
					     maxnode=0,
					     verbose=FALSE,
					     file='log.txt',
					     maxcores=1
					)
}




%\usage{complete_restart(g1,g2,taxnames,resultfile,size,start,end,\dots)}
\arguments{
  \item{g1}{the original network X}
  \item{g2}{the augmented network Y with additional nodes (all the original nodes from X must be present in the augmented network Y)}
  \item{taxnames}{the taxon name of the nodes added to the original graph (default: we select all nodes that are not in g1)}
  \item{resultfile}{filename to which results are saved after each group (default='result.txt')}
  \item{size}{the number of paths to investigate at each iteration of the algorithm (default=1000)}
  \item{start}{the starting group of path lengths to investigate (default=1)}
  \item{end}{the ending group of path lengths to investigate (default=0 \emph{i.e.} run until the end)}
  \item{maxnode}{the maximum number of augmented nodes in network Y to take into account (default=0, no maximum number of augmented nodes to take into account)} 	  
  \item{maxtime}{the maximum search time per pathway (default=3600 seconds)} 	
  \item{maxdistance}{the maximum search distance to the added nodes in network Y (default=0, no maximum distance for augmented nodes to take into account)} 	
  \item{verbose}{flag to save into a file additionnal informations regarding the obtained paths (default=FALSE, the file name should be indicated)}
  \item{file}{filename to save the additionnal informations}  
  \item{maxcores}{maximum number of cores to use (default=1, use 0 to use half of the total cores available )}
}

\value{This function returns the number of paths of each category : \emph{Shortcuts}, \emph{Detours}, \emph{Dead ends}, \emph{Equal} paths and unreachable nodes.\cr If some of the augmented nodes are not visited because of the limits defined by \code{maxnode}, \code{maxtime} or \code{maxdistance}, a path can also be classified as 'Detour or Dead End'. \cr
\cr
The format of the \code{resultfile} is as follow (tab-separated-value): \cr
Columns: \cr
1.  Group number \cr
2.  Unreachable nodes \cr
3.  Shortcuts \cr
4.  Egals \cr
5.  Detours \cr
6.  Dead ends \cr
7.  Undefined detour or dead end* \cr
8.  Total path evaluated \cr
9.  User time (seconds) \cr
10. System time (seconds) \cr
11. Real time (seconds) \cr
12. Unreachable total (up to this group) \cr
13. Shortcut total (up to this group) \cr
14. Egal total (up to this group) \cr
15. Detour total (up to this group) \cr
16. Dead end total (up to this group) \cr
17. Undefined total (up to this group) \cr
18. Total path evaluated (up to this group) \cr
19. Total user time (up to this group) \cr
20. Total system time (up to this group) \cr
21. Total real time (up to this group) \cr
}

%\author{Etienne Lord, Margaux Le Cam, Eric Bapteste, Vladimir Makarenkov  and Francois-Joseph Lapointe}

\examples{
## Searching the sample data (containing 11 original nodes and 3 augmented nodes)
data(Sample_1)
## Run from group number 1 to 2
complete_restart(g1,g2,size=5, start=1,end=2)
}

%\keyword{ }

