\encoding{ISO-8859-2}
\name{load_network}
\alias{load_network}
\title{is a helper function to load networks from files}
\description{This function allows the user to load the network files. Note that standard \emph{igraph} functions can be also used.}
\usage{load_network(filename_or_df, 
					filename_tax_or_df, 
					edge_weight,
					directed=FALSE)}
\arguments{
   \item{filename_or_df}{the file containing an undirected list of nodes in tabular format (node1, node2, edge weight) }
  \item{filename_tax_or_df}{an additional file containing taxa information for each node}
  \item{edge_weight}{how to treat the edge weight: 'equal', 'proportional' or 'inverse' (default: equal)}
  \item{directed}{treat input data as directed graph (default: FALSE)}
}

\value{This function return generic \emph{igraph} objects suitable for use in the SDDE package.}

%\author{Etienne Lord, Margaux Le Cam, Eric Bapteste, Vladimir Makarenkov  and Francois-Joseph Lapointe}

\examples{
	## Load a network from files. 
	\dontrun{
		network=load_network('graph.txt','graph_tax.txt','equal')
		info_network(network$g1,network$g2);
	}
	# We expect the network to be a series of nodes as tab-separated values.
	#
	# Example graph.txt:
	# node1	node2  edge weight	
	# x1	x2	   	1
	# x2	x3		1
	# x3	x6		1
	# x1	x5		1
	# x5	x6		1	
	#
	# Example graph_tax.txt
	# x1	plasmid 
	# x2	plasmid
	# x3	bacteria
	# x4	plasmid
	# x5	bacteria
	# x6	virus
}

%\keyword{Load network, Tabular separated value, tsv}

