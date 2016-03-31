\encoding{ISO-8859-2}
\name{plot_network}
\alias{plot_network}
\title{Helper function to display an illustration of a network.}
\description{This function display a representation of the networks. }
\usage{plot_network(g1,g2,layout,taxnames)}
\arguments{  
  \item{g1}{the original network X}
  \item{g2}{the augmented network Y with additional nodes (all the original nodes from X must be present in the augmented network Y)}
  \item{layout}{\emph{igraph} layout function (default=layout.kamada.kawai)}
  \item{taxnames}{the taxon name of the nodes added to the original graph. By default, we select all nodes that are not in g1. Note that in order to display all the different taxa groups, you must use 'allgroup' as the \emph{taxnames} argument}
}

%\author{Etienne Lord, Margaux Le Cam, Eric Bapteste, Vladimir Makarenkov and Francois-Joseph Lapointe}


%\keyword{Augmented network, plot function}

