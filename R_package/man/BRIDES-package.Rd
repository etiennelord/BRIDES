\encoding{utf-8}
\name{BRIDES-package}
\alias{BRIDES-package}
\docType{package}

\title{BRIDES: Characterize path types in evolving networks}

\description{This package allows the user to follow the evolution of an original network X into an augmented network Y by counting the number of Breakthroughs, Roadblocks, Impasses, Detours, Equal paths and Shortcuts (BRIDES) in the network Y. \cr 
			 \cr
			 Six different types of paths in Y that \emph{pass by at least one augmented node} (\emph{i.e.} a node that exists in Y but not in X) can be defined as follows: \cr
			 \cr
			1) \emph{Breakthrough} is a path that is impossible in network X, but is possible in network Y. \cr
			\cr
			2) \emph{Roadblock} is a path that is possible in network X, but is impossible in network Y. \cr
			\cr
			3) \emph{Impasse} is a path that is impossible in both networks, X and Y. \cr
			\cr
			4) \emph{Detour} is a path that is shorter in network X than in network Y. \cr
			\cr
			5) \emph{Equal} is a path that has the same length in networks X and Y. \cr
			\cr
			6) \emph{Shortcut} is a path that is longer in network X than in network Y.
}

\details{
	\tabular{ll}{
		Package: \tab BRIDES\cr
		Type: \tab Package\cr
		Version: \tab 1.2.-\cr
		Date: \tab 2020-09-05\cr
		License: \tab GPL-3\cr
		URL: \tab https://github.com/etiennelord/BRIDES \cr
		Maintainer: \tab Etienne Lord <m.etienne.lord@gmail.com>, \cr \tab François-Joseph Lapointe <francois-joseph.lapointe@umontreal.ca>
	}

	Function \code{\link{BRIDES}} compute the statistics between the two networks X and Y.
	
	Function \code{\link{export_network}} save to a file the selected network.
	
	Function \code{\link{load_network}} loads a network from a list of edges (tab-separated).
	
	Function \code{\link{save_network}} and \code{\link{plot_network}} creates rendering of network X and network Y in various file formats (svg, png and eps).
	
	
}


\author{Etienne Lord, Cindy Bouchard, Vladimir Makarekov and François-Joseph Lapointe}

\references{Lord, E., Le Cam, M., Bapteste, É., Méheust, R., Makarenkov, V., & Lapointe, F. J. (2016). BRIDES: a new fast algorithm and software for characterizing evolving similarity networks using breakthroughs, roadblocks, impasses, detours, equals and shortcuts. PloS one, 11(8), e0161474.}

\keyword{package}
\keyword{network analysis}
\keyword{pathways}
\keyword{evolving networks}

\seealso{
\code{\link{BRIDES}},
\code{\link{load_network}},
\code{\link{save_network}}
\code{\link{export_network}}
\code{\link{plot_network}}
}

%\examples{}
