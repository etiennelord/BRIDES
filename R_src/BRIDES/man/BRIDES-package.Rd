\encoding{utf-8}
\name{BRIDES-package}
\alias{BRIDES-package}
\alias{BRIDES}
\docType{package}

\title{Analyse of Six Pathways Types in Networks}

\description{BRIDES compares the evolution of an original network X to an augmented network Y by counting the number of shortest paths contains at least one augmented node. \cr 
			 \cr
			 Six different types of paths in Y that \emph{pass by at least one augmented node} (\emph{i.e.} a node that exists in Y but not in X) can be defined as follows: \cr
			 \cr
			 1) \emph{Breakthrough} is a path is impossible in the network X, that become in possible in network Y. \cr 
			 \cr
			 2) \emph{Roadblock} is a path that is possible in network X, but is impossible in network Y. \cr 
			 \cr
			 3) \emph{Impasse} is a path is impossible in both networks X and Y. \cr 
			 \cr
			 4) \emph{Detour} is a path that is shorter in network X but is longer in network Y. \cr
			 \cr
			 5) \emph{Equal} is a path that has the same length in X an Y. \cr	
			 \cr
			 6) \emph{Shorcut} is a path that is longer in network X than in network Y. \cr				 
			 \cr 
			 A C++ version for bigger networks is available at : \link{https://github.com/etiennelord/BRIDES} 
			 }
\details{
	\tabular{ll}{
		Package: \tab BRIDES\cr
		Type: \tab Package\cr
		Version: \tab 0.9.0\cr
		Date: \tab 2016-03-22\cr
		License: \tab GPL-2\cr
		Maintainer: \tab Etienne Lord <m.etienne.lord@gmail.com>, \cr \tab François-Joseph Lapointe <francois-joseph.lapointe@umontreal.ca>
	}

	Function \code{\link{BRIDES}} computes the BRIDES statistics for two networks.
	
	Function \code{\link{export_network}} export to an edgelist file the specified network.

	Function \code{\link{load_network}} loads a network from a list of edges  (tab-separated) which is compatible with the \code{\link{BRIDES}} function.

}


\author{Etienne Lord, Margaux Le Cam, Éric Bapteste, Vladimir Makarenkov and François-Joseph Lapointe}



%\references{}

\keyword{package}
\keyword{network analysis}
\keyword{pathways}
\keyword{genome similarity networks}

\seealso{
\code{\link{BRIDES}},
\code{\link{export_network}}
}

%\examples{}
