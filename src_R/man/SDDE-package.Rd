\encoding{utf-8}
\name{SDDE-package}
\alias{SDDE-package}
\alias{SDDE}
\docType{package}

\title{Shortcuts, Detours and Dead Ends (SDDE) Path Types in Genome Similarity Networks}

\description{SDDE compares the evolution of an original network X to an augmented network Y by counting the number of shortcuts, detours, dead ends, equal paths and disconnected nodes. \cr 
			 \cr
			 Four different types of paths in Y that \emph{pass by at least one augmented node} (\emph{i.e.} a node that exists in Y but not in X) can be defined as follows: \cr
			 \cr
			 1) \emph{Shortcut} is a path between two nodes, originally present in X, that is shorter in Y than in X. \cr 
			 \cr
			 2) \emph{Detour} is a path such that it connects two nodes originally present in X and is longer in Y than in X. \cr
			 \cr
			 3) \emph{Dead End} is a path that exists in X, but is impossible in Y. \cr
			 \cr
             4) \emph{Equal} is a path that has the same length in X an Y. \cr			 
			 \cr 
			 The numbers of \emph{disconnected} nodes is also reported. \cr
}

\details{
	\tabular{ll}{
		Package: \tab SDDE\cr
		Type: \tab Package\cr
		Version: \tab 1.0.1\cr
		Date: \tab 2015-11-06\cr
		License: \tab GPL-2\cr
		Maintainer: \tab Etienne Lord <m.etienne.lord@gmail.com>, \cr \tab François-Joseph Lapointe <francois-joseph.lapointe@umontreal.ca>
	}

	Function \code{\link{complete_network}} computes the number of paths in each category.
	
	Function \code{\link{complete_restart}} is similar to the \code{\link{complete_network}} function but it identifies the paths in a predefined order.

	Function \code{\link{complete_trace}} is similar to the \code{\link{complete_network}} function but allow the enumeration of nodes and taxa in a path.
	
	Function \code{\link{sample_network}} computes the number of paths in each category for a given set of nodes of the original network X.

	Function \code{\link{random_network}} creates a random network X and an augmented network Y using Erdos-Renyi or Barabási–Albert models.

	Function \code{\link{save_network}} and code{\link{plot_network}} creates rendering of network X and network Y in various file formats (svg, png and eps).

	Function \code{\link{load_network}} loads a network from a list of edges  (tab-separated) which is compatible with the \code{\link{complete_network}}, \code{\link{complete_restart}} and \code{\link{sample_network}} function.

}


\author{Etienne Lord, Margaux Le Cam, Éric Bapteste, Vladimir Makarenkov and François-Joseph Lapointe}



%\references{}

\keyword{package}
\keyword{network analysis}
\keyword{pathways}
\keyword{genome similarity networks}

\seealso{
\code{\link{complete_network}},
\code{\link{complete_restart}},
\code{\link{complete_trace}},
\code{\link{sample_network}},
\code{\link{load_network}},
\code{\link{save_network}}
\code{\link{plot_network}}
}

%\examples{}
