\encoding{ISO-8859-2}
\name{sample_path}
\alias{sample_path}
\title{is a helper function that creates a vector of non-repeating pathways to investigate}
\description{This internal function create a randomized vector of paired node numbers of a given length without repeat.}
\usage{sample_path(n, max_len=10, old_path,directed=FALSE)}
\arguments{ 
  \item{n}{the number of pathways}
  \item{max_len}{the total number of nodes (default=10)}
  \item{old_path}{vector of already sampled pairs of node numbers}
  \item{directed}{is it a directed or undirected network}
}
\value{This function return a randomize vector of 2 x \emph{n} node numbers.}
%\author{Etienne Lord, Margaux Le Cam, Eric Bapteste, Vladimir Makarenkov and Francois-Joseph Lapointe}
\examples{
	## Draw 5 samples of 10 pathways out of 100 nodes without repeated pathways
	\dontrun{
		old_path <- c();
		for (i in 1:5) {
		r <- sample_path(10,100, old_path);
			old_path <- c(old_path, r);
		}
	}
}

%\keyword{sampling}
