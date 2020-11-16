\encoding{ISO-8859-2} 
\name{random_network}
\alias{random_network}
\title{creates random augmented networks X and Y}
 \description{This function allows the user to create random network X and an associated augmented network Y using either the Erdos-Renyi model or the Barabasi-Albert model.}
\usage{random_network(original_node, additional_node, ngroup, edge_ratio,total_edge, type)}
\arguments{
  \item{original_node}{the number of nodes in the original network X (default=25)}
  \item{additional_node}{the number of additional node in network Y (default=5)}
  \item{ngroup}{the number of additional taxa groups in network Y (default=1)}
  \item{edge_ratio}{the edge to node ratio (default=between 1 and 5)}
  \item{total_edge}{the number of edges in the fixed model}
  \item{type}{Either 'erdos' for the Erdos-Renyi model,
              'barabasi' for the Barabasi-Albert model,
			  'watts' for the Watts-Strogatz model or
			  'fixed' fixed model allowing a defined number of edges
			  }
}

\value{
Return a data.frame containing \emph{g1}, \emph{g2}, the \emph{total_nodes} and the \emph{total_edges} numbers of network Y and the \emph{total_original_nodes} number of network X.
}

%\author{Etienne Lord, Margaux Le Cam, Eric Bapteste, Vladimir Makarenkov and Francois-Joseph Lapointe}

\examples{
	## Create a small random_network
	\dontrun{
		random_network()
	}
	## Expected result:
	#
	#  $g1
	#   IGRAPH UNW- 25 23 -- Erdos renyi (gnm) graph
	#   + attr: name (g/c), type (g/c), loops (g/l), m (g/n), name (v/c), tax
	#    (v/c), weight (e/n)
	#
	#  $g2
	#   IGRAPH UNW- 30 30 -- Erdos renyi (gnm) graph
	#   + attr: name (g/c), type (g/c), loops (g/l), m (g/n), name (v/c), tax
	#    (v/c), weight (e/n)
	#
	#   $total_nodes
	#	[1] 30
	#
	#	$total_edges
	# 	[1] 30
	#
	#   $total_original_nodes
	#   [1] 25
	#
	## Create two networks using the Erdos-Renyi model with 100 nodes in network X 
	## and 10 additional nodes in network Y of 3 types. 
	random_network(100,10,3);
	## Create a random networks of 20 +10 additional node using the  Barabsi-Albert model
	## and compute the corresponding SDDE path types.
	## l <- random_network(20,10,ngroup=1,vertex_ratio=1, type='barabasi');
	## complete_network(l$g1, l$g2);
}

