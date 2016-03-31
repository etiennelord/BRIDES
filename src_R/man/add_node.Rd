\encoding{ISO-8859-2}
\name{add_node}
\alias{add_node}
\title{Add node to an existing network}
\description{This function add a node to a network.}
\usage{add_node(g,name="",taxname="2",number_of_edges=sample(2:5,1))}
\arguments{
  \item{g}{the network}
  \item{name}{the desired node name)}
  \item{taxname}{the taxon name of the nodes added to the network (default: "2")}
  \item{number_of_edges}{the number of edge connecting the new node to the network}
}

\value{This function return the new network. Note that multiple edges between two nodes will be removed.}

%\author{Etienne Lord, Margaux Le Cam, Eric Bapteste, Vladimir Makarenkov  and Francois-Joseph Lapointe}

\examples{
   ## Add a node named XY to network g1 with 6 edges
   \dontrun{
   data(Sample_1)
   add_node(g1,name="XY","2",number_of_edges=6))
   } 
}

%\keyword{reachability, networks}

