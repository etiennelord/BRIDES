# BRIDES
This program allows the user to follow the evolution of an original network X into an augmented network Y by counting the number of Breakthroughs, Roadblocks, Impasses, Detours, Equal paths and Shortcuts (BRIDES) in the network Y. 
## About

* New version 2020 - Cindy Bouchard (Université de Montréal) 
* This new version (see the release section) allows scenarios simulation and search using genetics algorithm and more...  

This repository includes the R, R package and C++ source codes of BRIDES, examples of simple evolving networks and examples of real genome similarity networks which can be given as input of our program

Given an original network X and its augmented network Y (with additional nodes and edges), the BRIDES program calculates the number of Breakthroughs, Roadblocks, Impasses, Detours, Equal paths and Shortcuts characterizing the evolution of X into Y.

Six types of paths are thus possible:

1.	Breakthrough: a path that is impossible in network X, but is possible in network Y. 
2.	Roadblock: a path that is possible in network X, but is impossible in network Y.
3.	Impasse: a path that is impossible in both networks, X and Y.
4.	Detour: a path that is shorter in network X than in network Y.
5.	Equal: a path that has the same length in networks X and Y.
6.	Shortcut: a path that is longer in network X than in network Y.

The R and R package versions of BRIDES is suitable for small networks (with less than 1,000 nodes), while the C++ version can handle millions of nodes.

One novelty of the R package version (2020) is the addition of scenarios maximize a search for the best additions of nodes to the starting network. 

## Installation

1. Clone the git repository to your computer or download the corresponding zip file or download the R package version in the release section

```
git clone https://github.com/etiennelord/COMPONENT-GRAPHER.git
```

2. Follow the following steps (depending of the BRIDES version):

A) For the R standalone version

From the R_src directory, enter the R shell. Then execute: 

```
source("BRIDES_2022.R")
```

B) For the R package version

Enter the R shell. Then execute:

```
install.packages("BRIDES_1.2.1.tar.gz", repos = NULL, type = "source")
```

Alternatively, use the "Install from Package Archive File" in the R Studio package menu.

C) For the C++ version

This version requires the use of **gcc version 4.5** or higher for the C++ compiler. 
On MacOSX, the use of OpenMP will require the installation of [Homebrew](http://brew.sh/) and the execution of:

```
brew install gcc --without-multilib
```

- From the Cpp_src directory, execute:

```C
// Using a Makefile
make -f Makefile.VERSION

// Manual compilation
g++ -O3 main.cpp -o bride -fopenmp
```
where the Makefile.VERSION is either Makefile.Linux, or Makefile.MacOSX, or Makefile.MACOSX_wo_OPENMP, or Makefile.Windows . 

## Samples and Gene similarity networks

**Sample networks** can be found in both [R](https://github.com/etiennelord/BRIDES/tree/master/R_src/sample) and [C++](https://github.com/etiennelord/BRIDES/tree/master/Cpp_src/sample) source directories.

**Genome similarity networks** (70% cover, 90% similarity threshold, minimum BLAST e-value 10e-5) from 2,094,593 nucleotide sequenes of archaea, bacteria and eukaryotes can be found in the [GenomeNetwork](https://github.com/etiennelord/BRIDES/tree/master/GenomeNetwork) directory, along with BRIDES algorithm results for those network using the parameters: (MaxPathNumber=100, MaxDistance=100, MaxNode=100, please see [manual](https://github.com/etiennelord/BRIDES/blob/master/Manual/BRIDES_User_Guide.pdf) for explanations).  

**Simulation scripts** for different networks can be found in the [simulation](https://github.com/etiennelord/BRIDES/tree/master/Simulation) directory.

## Usage 

Also see the program manual in the **Manual** directory.

A) Using the R version:

```R
## Load the functions and dependencies
source("BRIDES.R") or library(BRIDES)

## Load the sample network
t0<-load_network("sample/t0.txt")
t1<-load_network("sample/t1.txt")
u0<-load_network("sample/u0.txt","sample/u0.attr.txt", directed=T)

## Execute BRIDES on networks t0 and t1
BRIDES(t0,t1)

## Execute BRIDES on directed networks U0 with an attribute file
BRIDES(u0, A="2")
```

## For the scenarios (R package version - 2020 - Cindy Bouchard)

* Scenarios: Exhaustive search 
```
# Searching for the scenario favoring: Breakthough [B=1 R=0 I=0 D=0 E=0 S=0]
# The scenario is to add 1 to 2 new nodes to network X
# The search will be un undirected and unweighted networks
# This search mode will evaluate ALL possible scenarios
# Note: this search mode could take some time on bigger networks
data(networkX)
data(networkY)
BRIDES(networkX, networkY,runmode='exhaustive',max_additional=2, wt=c(1,0,0,0,0,0))
```

* Scenarios: Genetics search 
```
# Searching for the scenario favoring:Breakthough [B=1 R=0 I=0 D=0 E=0 S=0]
# The scenario is to add 1 to 2 new nodes to network X
# The search will be un undirected and unweighted networks
# This search mode will use a genetic algorithm to converge to a solution
# using artificial crossing over between the best local solutions.
data(networkX)
data(networkY)
results<-BRIDES(networkX, networkY,runmode='genetics',max_additional=2, wt=c(1,0,0,0,0,0),
```

* Scenarios: Stepwise search 
```
# Searching for the scenario favoring:Breakthough [B=1 R=0 I=0 D=0 E=0 S=0]
# The scenario is to add 1 to 2 new nodes to network X
# The search will be un undirected and unweighted networks
# This search mode add iteratively the first best node, then try to add
# the second best node, etc.
data(networkX)
data(networkY)
results<-BRIDES(networkX, networkY,runmode='stepwise',max_additional=2, wt=c(1,0,0,0,0,0))
```

B) Using the C++ version: 

```
## Execute BRIDES on networks t0 and t1
./brides -X=sample/t0.txt -Y=sample/t1.txt

## Execute BRIDES on directed network U0 with an attribute file
./brides -X=sample/U0.txt -A=sample/U0.attr.txt -K=2 -directed
```

## Dependencies

The R version of BRIDES depends on the following R libraries:

##### [fastmap](https://cran.r-project.org/web/packages/fastmap/index.html)  
Used in the R version to help list operations.
##### [igraph](http://igraph.org/r/)
A collection of network handling tools including shortest-path analysis.   
##### [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)
For-each parallel adaptor for the 'parallel' R package.

The C++ version of the program is optimized for the use of OpenMP and dispatching tasks to different threads.
##### [OpenMP](http://openmp.org/wp/)  

