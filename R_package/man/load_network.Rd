\encoding{utf-8}
\name{load_network}
\alias{load_network}
\title{Helper function to load a network into igraph from a file}
\description{This function load/import a network to an\emph{igraph} object from a file or data.frame.
             In particular, it can handle single node and also attributes}
\usage{load_network(filename_or_df, filename_attributes_or_df, edge_weight,sep,header, directed)}
\arguments{ 
  \item{filename_or_df}{The filename or data.frame containing the list of node}
  \item{filename_attributes_or_df}{The filename or dataframe containing annotation fo the nodes}
  \item{directed}{Do we treat the network as directed, default:FALSE}
  \item{edge_weight}{How to handle edge weights: 'equal' (unweighted, 1.0), 'proportional' (load from file), 'inverse' (1/weight) }
  \item{header}{Do the files contains header, default: FALSE}
   \item{sep}{The text separator separation columns in the file, default: 'tab'}
}

\examples{
	\dontrun{
	##########################################################################
	#                    Load a directed network                             #
	##########################################################################
	# Note: see https://github.com/etiennelord/BRIDES for the sample files
	U0<-load_network("sample/U0.txt","sample/U0.attr.txt",directed=T)
	plot_network(U0, attributes="2")
	##########################################################################
	#                    Load a weigthed network                             #
	##########################################################################
	U0<-load_network("sample/U0.txt","sample/U0.attr.txt",directed=T, 
	                  edge_weight="proportional")
	plot_network(U0, attributes="2")
	
	}
}

\keyword{load function}
\keyword{import function}

