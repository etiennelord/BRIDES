\encoding{utf-8}
\name{save_network}
\alias{save_network}
\title{is an helper function to save an illustration of a network into a file}
\description{This function saves to a file a representation of the given networks}
\usage{save_network(g1,g2,filename,layout,taxnames,mode,imagesize)}
\arguments{ 
  \item{g1}{the original network X}
  \item{g2}{the augmented network Y with additional nodes (all the original nodes from X must be present in the augmented network Y)}
  \item{filename}{the name of the file to save the network. Note that the file extension will automatically be added}
  \item{layout}{\emph{igraph} layout function (default=layout.kamada.kawai)}
  \item{taxnames}{the taxon name of the nodes added to the original graph. By default, we select all nodes that are not in g1. Note that in order to display all the different taxa groups, you must use 'allgroup' as the \emph{taxnames} argument}
  \item{mode}{type of file to create either: 'png', 'svg', or 'eps' (default='png')}
  \item{imagesize}{size of the image in pixels (default=800 px)}  
}

\examples{
	\dontrun{
		data(Sample_1)
		## Save the networks to file using the default layout (layout.kamada.kawai)
		## and size (1200px)
		save_network(g1,g2,filename="sample1",imagesize=1200)
		## Save the networks to a file using the layout.fruchterman.reingold layout, 
		## and in 'svg' format
		save_network(g1,g2,filename="s", mode="svg",layout=layout.fruchterman.reingold)
		## Save a network containing different groups
		gsample <- random_network(ngroup=4)
		save_network(gsample$g1,gsample$g2,filename="gsample", taxnames="allgroup")
	}
}

\keyword{Augmented network, save function, png, svg}

