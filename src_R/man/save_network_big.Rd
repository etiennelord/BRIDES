\encoding{ISO-8859-2}
\name{save_network_big}
\alias{save_network_big}
\title{Helper function to save an illustration of a network to file with more than 50 nodes.}
\description{This function save to file a representation of the networks. Note that it will be automatically called by the save_network function.}
\usage{save_network_big(g1,g2,filename,layout,taxnames,mode,imagesize)}
\arguments{  
   \item{g1}{the original network X}
  \item{g2}{the augmented network Y with additional nodes (all the original nodes from X must be present in the augmented network Y)}
  \item{filename}{the name of the file to save the network. Note that the file extension will automatically be added}
  \item{layout}{\emph{igraph} layout function (default=layout.kamada.kawai)}
  \item{taxnames}{the taxon name of the nodes added to the original graph. By default, we select all nodes that are not in g1. Note that in order to display all the different taxa groups, you must use 'allgroup' as the \emph{taxnames} argument}
  \item{mode}{type of file to create either: 'png', 'svg', or 'eps' (default='png')}
  \item{imagesize}{size of the image in pixels (default=2500 px)} 
}

%\author{Etienne Lord, Margaux Le Cam, Eric Bapteste, Vladimir Makarenkov  and Francois-Joseph Lapointe}


%\keyword{Augmented network, save function, png, svg}

