\encoding{utf-8}
\name{save_network}
\alias{save_network}
\title{This is an helper function to save an illustration of some network(s) to a file}
\description{This function saves to a file a representation of the given networks}
\usage{save_network(g1,g2,filename,layout,attributes,mode,imagesize, 
       main, legend_position, label_filename)}
\arguments{ 
  \item{g1}{the original network X}
  \item{g2}{the augmented network Y with additional nodes (all the original nodes from X must be present in the augmented network Y)}
  \item{filename}{the name of the file to save the network. Note that the file extension will automatically be added}
  \item{layout}{\emph{igraph} layout function (default=layout.kamada.kawai)}
  \item{main}{the title of the picture}
  \item{attributes}{the taxon name of the nodes added to the original graph. 
        By default, we select all nodes that are not in g1. 
		Note that in order to display all the different attributes groups, 
		you must use 'allgroup' as the \emph{atrtibutes} argument}
  \item{mode}{type of file to create either: 'png', 'svg', or 'eps' (default='png')}
  \item{imagesize}{size of the image in pixels (default=800 px)}
  \item{legend_position}{position of the legend, default: "bottomright"}  
  \item{label_filename}{currently in development.}  
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
		save_network(gsample$g1,gsample$g2,filename="gsample", attributes="allgroup")
	}
}

\keyword{save function}
\keyword{png}
\keyword{svg}
