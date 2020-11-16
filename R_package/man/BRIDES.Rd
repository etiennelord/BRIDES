\encoding{utf-8}
\name{BRIDES}
\alias{BRIDES}
\title{Main function to calculate the evolving networks}
\description{This function calculates the different path types between two networks}
\usage{BRIDES(X, Y = NULL, src = "default", dest = "default",
                   attributes = "default", random = 0, maxdistance = 100,
                   maxtime = 100, maxnode = 100, maxcores = 1, outfile =
                   "", size = 1000, first = 0, last = 0, weighted = TRUE,
                   runmode = "default", path_information = TRUE, verbose
                   = TRUE, min_additional = 1, max_additional = 0,
                   max_iters = 100, mutationChance = 0.1, wt = c(3, -1,
                   -2, 1, 1, 2))}


\arguments{ 
  \item{X}{the original network X (an igraph network object)}
  \item{Y}{the augmented network Y with additional nodes (all the original nodes from X must be present in the augmented network Y)}
  \item{src}{Specify a source node (must be used with the dest argument)}
  \item{dest}{Specify a destination node (must be used with the src argument.)}
  \item{attributes}{If only one network is given, attribute which can be specified for an added nodes.\cr Must be in the "tax" attribute of the igraph network: e.g. V(Y)$tax}
  \item{first}{Specify the starting path number for the BRIDES computation if distributed}
  \item{last}{Specify the ending path number for the BRIDES computation if distributed}
  \item{min_additional}{Specify for the scenarios the minimum number of additional nodes to consider}
  \item{max_additional}{Specify for the scenarios the maximum number of additional nodes to consider}
  \item{maxdistance}{Maximum distance to an added node, default: 100}
  \item{max_iters}{Specify the number of iterations to search for the 'genetics' runmode}
  \item{maxnode}{Maximum number of random additional nodes in network Y to consider}
  \item{maxcores}{Specify the number of compute cores (threads) to use for the calculation}	
  \item{maxtime}{Specify the maximum time (in seconds) to individual path search, default: 100 seconds}
  \item{mutationChance}{Specify mutation chance for 'genetics' runmode, default=0.1}
  \item{outfile}{Specify an output file for the each computed paths between each pairs of nodes}
  \item{path_information}{Record path information (disabling this lead to faster execution), default: TRUE }
  \item{random}{If specified, will sample a number of random pathways}
  \item{runmode}{Runmode: 'default' (BRIDES statistics only)
					  \cr 'exhaustive' scenario-based search the maximum BRIDES score using the wt scoring scheme
					  \cr 'genetics' scenario-based genetics algorithm search for the best BRIDES score
					  \cr 'stepwise' scenario-based incremental search adding the best one node, than the second best to 
					  \cr maximize the BRIDES score using the wt scoring scheme 
					  \cr Note : for those scenarios, a min_additional and max_additional node must be provided
					  \cr Note^{2}: for the 'genetics' runmode, max_iters should be large enough.
					  }
  \item{size}{Specify the number of paths to consider if computation is distributed}
  \item{verbose}{Verbose to the screen some information, default: TRUE }
  \item{weighted}{Specify if we consider the edge's lengths, default: TRUE }
  \item{wt}{Scoring sheme [B R I D E S] for the runmodes 'exhaustive','genetics','incremental' }
}

\examples{
\dontrun{
	data(networkX)
	data(networkY)
	set.seed(1)
	##########################################################################
	#           Basic usage of a BRIDES network analysis in R                #
	##########################################################################	
	results<-BRIDES(networkX,networkY, outfile="outfile.txt")
	print(results)
	# [...Results...]
	#    B      R      I      D      E      S (utime stime) 
    #    9      2      20     26     5      4   0.07   1.42
	# Note: pathways and BRIDES statistics will be saved to outfile.txt
	##########################################################################
	#              Saving path information to file outfile.txt               #
	##########################################################################
	BRIDES(networkX,networkY, outfile="outfile.txt")
	##########################################################################
	# Searching for a specific path between the source and destination nodes #
	##########################################################################	
	BRIDES(networkX,networkY,src="x1",dest="x7")
	#Prunning unconnected nodes...
	#x1 -> x7 
	#original: 1 
	#augmented: 4 
	#[1] "x1"  "x14" "x6"  "x5"  "x7" 
	#[1] 0 1 0 0 0
	#$`from`
	#[1] "x1"
	#$to
	#[1] "x7"
	#$path_type
	#[1] "Detour"
	#$path_type0
	#[1] 4
	#$original_path_length
	#[1] 1
	#$augmented_path_length
	#[1] 4
	#$path
	#[1] "x1"  "x14" "x6"  "x5"  "x7" 
	#$path_visited_taxa
	#[1] 0 1 0 0 0
	##########################################################################
	#      Sampling 100 random paths in a Erdős–Rényi random network         #
	##########################################################################
	set.seed(25)
	#Create an original network X of 20 nodes with 5 additional in network Y
	g<-random_network(20,5,type="erdos") 
	BRIDES(g$g1,g$g2,random=100)
	##########################################################################
	#    Limit the distance in network Y using the maxdistance parameter     #
	##########################################################################
	BRIDES(networkX,networkY, maxdistance=1)
	# [...Results...]
	#    B      R      I      D      E      S (utime stime) 
    #    9      25     20     3      5      4   0.04   0.93  (with maxdistance=1)
    #    9      2      20     26     5      4   0.07   1.42 (without)
	#           *             *
	##########################################################################
	#      Default analysis with an attribute file and directed network      #
	##########################################################################
	# Note: see https://github.com/etiennelord/BRIDES for the sample files
	U0<-load_network("sample/U0.txt","sample/U0.attr.txt",directed=T)
	plot_network(U0, attributes="2")
	# look in vertex tax attributes for the additional nodes 
	results<- BRIDES(U0, attributes="2", outfile="directed.txt") 
	##########################################################################
	#                    Weigthed networks analysis                          #
	##########################################################################
	# Note: see https://github.com/etiennelord/BRIDES for the sample files
	U0<-load_network("sample/U0.txt","sample/U0.attr.txt",directed=T, 
	                  edge_weight="proportional")
	plot_network(U0, attributes="2")
	results<- BRIDES(U0, attributes="2", outfile="directed.txt") 
	##########################################################################
	#                    Weigthed networks analysis                          #
	##########################################################################
	# Note: if your networks already have weighted edges, 
	#       they will be taken into account.
	data(dX)
	data(dY) 
	plot_network(dX,dY,legend_position="topleft")
	BRIDES(dX,dY, outfile="directed2.txt") 
	# [ BRIDES Results ]
	#     B      R      I      D      E      S (utime stime) 
	#    5.00   3.00   4.00   0.00   0.00   0.00   0.03   0.30 
	# [ directed2.txt ]
	# src	dest	dist_x	dist_y	BRIDES	path	path_attributes
	# x2	x1	Inf	3	Breakthrough	x2,x4,y1,x1	0,0,1,0	
	# x5	x1	Inf	3	Breakthrough	x5,x4,y1,x1	0,0,1,0	
	# x4	x1	Inf	2	Breakthrough	x4,y1,x1	0,1,0	
	# x1	x2	1	Inf	Roadblock			
	# x5	x2	Inf	3	Breakthrough	x5,x4,y1,x2	0,0,1,0	
	# x4	x2	Inf	2	Breakthrough	x4,y1,x2	0,1,0	
	# x1	x5	Inf	Inf	Impasse			
	# x2	x5	Inf	Inf	Impasse			
	# x4	x5	Inf	Inf	Impasse			
	# x1	x4	2	Inf	Roadblock			
	# x2	x4	1	Inf	Roadblock			
	# x5	x4	Inf	1	Impasse			
	##########################################################################
	# Scenarios: Exhaustive search                                           #
	##########################################################################
	# Searching for the scenario favoring:Breakthough [B=1 R=0 I=0 D=0 E=0 S=0]
	# The scenario is to add 1 to 2 new nodes to network X
	# The search will be un undirected and unweighted networks
	# This search mode will evaluate ALL possible scenarios
	# Note: this search mode could take some time on bigger networks
	data(networkX)
	data(networkY)
	results<-BRIDES(networkX, networkY,runmode='exhaustive',max_additional=2, wt=c(1,0,0,0,0,0))
	# [...Results...]
	# Using scoring scheme [ B R I D E S ] : 1 0 0 0 0 0 
	# Exhaustive search for 6 iterations with min. : 1 and max.:2 additional nodes.
	# [Iteration 1/6] (score:0)
	# [Iteration 2/6] (score:9)
	# [Iteration 3/6] (score:0)
	# [Iteration 4/6] (score:9)
	# [Iteration 5/6] (score:0)
	# [Iteration 6/6] (score:9)
	# Exhaustive search done.
	#         B R I  D  E S X.utime stime. score best
	# x14     0 4 29 25 4 4    0.05   1.19     0     
	# x15     9 3 20 29 3 2    0.05   1.34     9    *
	# x13     0 5 29 32 0 0    0.04   1.30     0     
	# x14,x15 9 3 20 25 5 4    0.10   1.19     9    *
	# x14,x13 0 2 29 27 4 4    0.04   1.16     0     
	# x15,x13 9 2 20 30 3 2    0.05   1.28     9    *
	# * Best solutions found using criteria.
	##########################################################################
	# Scenarios: Genetics search                                            #
	##########################################################################
	# Searching for the scenario favoring:Breakthough [B=1 R=0 I=0 D=0 E=0 S=0]
	# The scenario is to add 1 to 2 new nodes to network X
	# The search will be un undirected and unweighted networks
	# This search mode will use a genetic algorithm to converge to a solution 
	# using artificial crossing over between the best local solutions.
	# Note: It is important to set a max_iters (iteration) in this mode
	# Note2: There is no garanty that a the maximum score will be evaluated 
	data(networkX)
	data(networkY)
	results<-BRIDES(networkX, networkY,runmode='genetics',max_additional=2, wt=c(1,0,0,0,0,0), 
	                max_iters=20)
	# [...Results...]
	# Using scoring scheme [ B R I D E S ] : 1 0 0 0 0 0 
	# Pre-training. Evaluating 3 solutions.
	# Optimization for 20 iterations with min.:1 max.:2 additional nodes.
	# [Iteration 1/20] (score:0)
	# [Iteration 2/20] (score:9)
	# [Iteration 3/20] (score:9)
	# No more suitable solutions found.
	# Genetics search done.
	#         B R  I  D E S X.utime stime. score best
	# x14     0 4 29 25 4 4    0.06   1.19     0     
	# x15     9 3 20 29 3 2    0.04   1.33     9    *
	# x13     0 5 29 32 0 0    0.03   1.31     0     
	# x14,x13 0 2 29 27 4 4    0.06   1.14     0     
	# x15,x13 9 2 20 30 3 2    0.03   1.27     9    *
	# x14,x15 9 3 20 25 5 4    0.04   1.20     9    *
	# * Best solutions found using criteria.
	##########################################################################
	# Scenarios: Stepwise search                                             #
	##########################################################################
	# Searching for the scenario favoring:Breakthough [B=1 R=0 I=0 D=0 E=0 S=0]
	# The scenario is to add 1 to 2 new nodes to network X
	# The search will be un undirected and unweighted networks
	# This search mode add iteratively the first best node, then try to add
	# the second best node, etc.
	data(networkX)
	data(networkY)
	results<-BRIDES(networkX, networkY,runmode='stepwise',max_additional=2, wt=c(1,0,0,0,0,0))
	# [...Results...]
	# Using scoring scheme [ B R I D E S ] : 1 0 0 0 0 0 
	# Stepwise search for maximum :2 additional nodes.
	# Starting adding new nodes:  1  of  2 
	# [Iteration 1/3] (in progress)
	# [Iteration 2/3] (in progress)
	# [Iteration 3/3] (in progress)
	# Starting adding new nodes:  2  of  2 
	# [Iteration 1/2] (in progress)
	# [Iteration 2/2] (in progress)
	# Stepwise search done.
	#         B R  I  D E S X.utime stime. score best
	# x14     0 4 29 25 4 4    0.03   1.16     0     
	# x15     9 3 20 29 3 2    0.07   1.34     9    *
	# x13     0 5 29 32 0 0    0.04   1.30     0     
	# x14,x15 9 3 20 25 5 4    0.04   1.19     9    *
	# x15,x13 9 2 20 30 3 2    0.05   1.28     9    *
	# * Best solutions found using criteria.
	#
	# Note that x14,x13 addition was not evaluated
	##########################################################################
	# Scenarios: Stepwise search (with 2,3) additional nodes                 #
	##########################################################################
	# The scenario is to add 1 or 2 new nodes to network X
	# Searching for the scenario using the default : [B=3 R=-1 I=-2 D=1 E=1 S=2]
	# 
	set.seed(42)
	g<-random_network(20,5,type="watts") # Create a random network
	plot_network(g$g1, g$g2)
	BRIDES(g$g1,g$g2,runmode='stepwise',max_additional=2)
	# [ Complete output ]
	# Using scoring scheme [ B R I D E S ] : 3 -1 -2 1 1 2 
	# Stepwise search for maximum :2 additional nodes.
	# Starting adding new nodes:  1  of  2 
	# [Iteration 1/5] (in progress)
	# [Iteration 2/5] (in progress)
	# [Iteration 3/5] (in progress)
	# [Iteration 4/5] (in progress)
	# [Iteration 5/5] (in progress)
	# Starting adding new nodes:  2  of  2 
	# [Iteration 1/4] (in progress)
	# [Iteration 2/4] (in progress)
	# [Iteration 3/4] (in progress)
	# [Iteration 4/4] (in progress)
	# Stepwise search done.
	#		   B  R   I  D E S X.utime stime. score best
	# x13     53 48  89  0 0 0    0.09   2.13   -67     
	# x17      0 48 142  0 0 0    0.13   1.67  -332     
	# x18      0 26 142 15 5 2    0.07   2.02  -286     
	# x23      0 48 142  0 0 0    0.11   1.66  -332     
	# x24      5 39 137  8 0 1    0.09   1.80  -288     
	# x13,x17 53 48  89  0 0 0    0.09   2.04   -67     
	# x13,x18 53 26  89 15 5 2    0.10   2.27   -21     
	# x13,x23 53 48  89  0 0 0    0.09   2.07   -67     
	# x13,x24 67 39  75  8 0 1    0.06   2.30    22    *
	#* Best solutions found using criteria.
	#
	# [ Results list object ]
	#$`population` (The evaluated scenarios)
	#  [,1] [,2] [,3] [,4] [,5]
	#1    1    0    0    0    0
	#2    0    1    0    0    0
	#3    0    0    1    0    0
	#4    0    0    0    1    0
	#5    0    0    0    0    1
	#6    1    1    0    0    0
	#7    1    0    1    0    0
	#8    1    0    0    1    0
	#9    1    0    0    0    1
	# $brides (the resulting score)
	#	   B  R   I  D E S (utime stime)
	# [1,] 53 48  89  0 0 0   0.09   2.13
	# [2,]  0 48 142  0 0 0   0.13   1.67
	# [3,]  0 26 142 15 5 2   0.07   2.02
	# [4,]  0 48 142  0 0 0   0.11   1.66
	# [5,]  5 39 137  8 0 1   0.09   1.80
	# [6,] 53 48  89  0 0 0   0.09   2.04
	# [7,] 53 26  89 15 5 2   0.10   2.27
	# [8,] 53 48  89  0 0 0   0.09   2.07
	# [9,] 67 39  75  8 0 1   0.06   2.30
	# $size (the maximum number of additional node)
	# [1] 5
	# $evaluations (the score according the the wt)
	# [1]  -67 -332 -286 -332 -288  -67  -21  -67   22
	# $additional_nodes (the combination of additionnal nodes)
	# [1] "x13" "x17" "x18" "x23" "x24" "x13,x17" "x13,x18" "x13,x23" "x13,x24"
	#
	# $final (the *BEST* scenario, there might be more than one)
	#[1] 1 0 0 0 1
	# $final_score (final scenario score)
	#[1] 22
	# $final_chromosome (the scenario position in the population)
	#[1] 9
	# $nodes (the available additional nodes)
	# [1] "x13" "x17" "x18" "x23" "x24"
	# $final_table (displayed table)
	# 		 B  R   I  D E S X.utime stime. score best
	# x13     53 48  89  0 0 0    0.09   2.13é   -67     
	# x17      0 48 142  0 0 0    0.13   1.67  -332     
	# x18      0 26 142 15 5 2    0.07   2.02  -286     
	# x23      0 48 142  0 0 0    0.11   1.66  -332     
	# x24      5 39 137  8 0 1    0.09   1.80  -288     
	# x13,x17 53 48  89  0 0 0    0.09   2.04   -67     
	# x13,x18 53 26  89 15 5 2    0.10   2.27   -21     
	# x13,x23 53 48  89  0 0 0    0.09   2.07   -67     
	# x13,x24 67 39  75  8 0 1    0.06   2.30    22    *
	#
	# $networks (a list of the *BEST* resulting networks)
	# $networks[[1]]
	# IGRAPH c9b4e9f UNW- 22 21 -- Watts-Strogatz random graph
	# + attr: name (g/c), dim (g/n), size (g/n), nei (g/n), p (g/n), 
	#   loops (g/l), multiple (g/l), name (v/c), tax (v/c),
	# | weight (e/n)
	# + edges from c9b4e9f (vertex names):
    # [...]
}
	#####################################################
	# BRIDES on a small network, single core            #
	#####################################################
	data(networkX)
	data(networkY)
	BRIDES(networkX,networkY,maxcores=1)
}

\keyword{BRIDES}
\keyword{Evolution}
\keyword{Scenario}

