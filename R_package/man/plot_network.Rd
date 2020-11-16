\encoding{utf-8}
\name{plot_network}
\alias{plot_network}
\title{This is an helper function to display a augmented network}
\description{This function display a representation of the given networks}
\usage{plot_network(g1, g2 = NULL, layout = layout.kamada.kawai, taxnames
                   = "", main = "", attributes = "", legend_position =
                   "bottomright", label_filename = "")}
\arguments{ 
  \item{g1}{the original network X}
  \item{g2}{the augmented network Y}
  \item{main}{the title of the picture}
  \item{layout}{\emph{igraph} layout function (default=layout.kamada.kawai)}
  \item{attributes}{the attributes of the nodes added to the original graph.
    By default, we select all nodes that are not in g1.}
  \item{taxnames}{select nodes with belongs to a specfic group (in the V(g)$tax)
  Note that in order to display all the different groups, 
  you must use 'allgroup' as the \emph{taxnames} argument}
  \item{legend_position}{position of the legend, default: "bottomright". Set to 'NULL' to hide it.}  
  \item{label_filename}{currently in development.}
}

\examples{
	\dontrun{
		data(networkX)
		data(networkY)
		plot_network(networkX, networkY, legend_position="topleft")
	}
}

\keyword{display function}
\keyword{plot}

