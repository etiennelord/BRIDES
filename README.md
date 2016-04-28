# BRIDES
This program allows the user to follow the evolution of an original network X into an augmented network Y by counting the number of Breakthroughs, Roadblocks, Impasses, Detours, Equal paths and Shortcuts (BRIDES) in the network Y. 
## About

This repository includes the R and C++ source codes of BRIDES, examples of simple evolving networks and examples of real genome similarity networks which can be given as input of our program

Given an original network X and its augmented network Y (with additional nodes and edges), the BRIDES program calculates the number of Breakthroughs, Roadblocks, Impasses, Detours, Equal paths and Shortcuts characterizing the evolution of X into Y.

Six types of paths are thus possible:

1.	Breakthrough: a path that is impossible in network X, but is possible in network Y. 
2.	Roadblock: a path that is possible in network X, but is impossible in network Y.
3.	Impasse: a path that is impossible in both networks, X and Y.
4.	Detour: a path that is shorter in network X than in network Y.
5.	Equal: a path that has the same length in networks X and Y.
6.	Shortcut: a path that is longer in network X than in network Y.

The R version of BRIDES is suitable for small networks (with less than 1,000 nodes), while the C++ version can handle millions of nodes.

## Installation

1. Clone the git repository to your computer or download the corresponding zip file

```
git clone https://github.com/etiennelord/COMPONENT-GRAPHER.git
```

2. Follow the following steps (depending of the BRIDES version):

A) For the R version

From the R_src directory, enter the R shell. Then execute: 

```
install.packages("SDDE")
source("BRIDES.R")
```

B) For the C++ version

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

## Usage 

Also see the program manual in the **Manual** directory.

A) Using the R version:

```R
## Load the functions and dependencies
source("BRIDES.R")

## Load the sample network
t0<-load_network("sample/t0.txt")
t1<-load_network("sample/t1.txt")
u0<-load_network("sample/u0.txt","sample/u0.attr.txt", directed=T)

## Execute BRIDES on networks t0 and t1
BRIDES(t0,t1)

## Execute BRIDES on directed networks U0 with an attribute file
BRIDES(u0, A="2")
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

##### [SDDE](https://cran.r-project.org/web/packages/SDDE/index.html)  
Wrapper to some useful graph manipulation functions of the 'igraph' library and the 'doParallel' R package.
##### [igraph](http://igraph.org/r/)
A collection of network handling tools including shortest-path analysis.   
##### [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)
For-each parallel adaptor for the 'parallel' R package.

The C++ version of the program is optimized for the use of OpenMP and dispatching tasks to different threads.
##### [OpenMP](http://openmp.org/wp/)  

