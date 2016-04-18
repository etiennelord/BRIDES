# BRIDES
Compares the evolution of an original network X to an augmented network Y by counting the number of Breakthroughs, Roadblocks, Impasses, Detours, Equal paths, and Shortcuts (BRIDES). 

## About

This is the R and C++ source code, some sample networks and some genome similarity networks for the BRIDES software.

The BRIDES software, given and  original network X  its augmented network Y (with additional nodes and edges) calculates the number of Breakthroughs, Roadblocks, Impasses, Detours, Equal paths and Shortcuts (BRIDES) characterizing the evolution of X into Y.

Six types of paths are thus possibles:

1. Breakthrough: pathway that is impossible in network X but possible in network Y. 
2. Roadblock: pathway that is possible in network X but impossible in network Y.
3. Impasse: pathway that is impossible in both networks, X and Y.
4. Detour: pathway that is shorter in network X than in network Y.
5. Equal: pathway that has the same length in networks X and Y.
6. Shortcut: pathway that is longer in network X than in network Y.

The R version is suitable for small networks (less than 1,000 nodes) while the C++ version can handle millions of nodes.

## Installation

1. clone the git repository to your computer or download the zip file

```
git clone https://github.com/etiennelord/COMPONENT-GRAPHER.git
```

2. Follow the appropriate step for the R or C++ version

A) For the R version

From the R_src directory, entre the R shell

```
install.packages("igraph")
install.packages("SDDE")
source("BRIDES.R")
```

B) For the C++ version

This version requires **gcc version 4.5** or higher. 
On MacOSX, to use OpenMP, that likely required the installation of [Homebrew](http://brew.sh/) then running the command:

```
brew install gcc --without-multilib
```

- From the Cpp_src directory

```C
// Using a Makefile
make -f Makefile.VERSION

// Manual compilation
g++ -O3 main.cpp -o bride -fopenmp
```
Where the Makefile.VERSION is either Makefile.Linux, Makefile.MacOSX, Makefile.MACOSX_wo_OPENMP or Makefile.Windows . 

## Usage 

Refer to the manual in the ** Manual ** directory.

A) Using the R version

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

B) Using the C++

```
## Execute BRIDES on networks t0 and t1
./brides -X=sample/t0.txt -Y=sample/t1.txt

## Execute BRIDES on directed networks U0 with an attribute file
./brides -X=sample/U0.txt -A=sample/U0.attr.txt -K=2 -directed
```

## Dependencies

The R version of BRIDES depends on the following R libraries:

##### [SDDE](https://cran.r-project.org/web/packages/SDDE/index.html)  
Wrapper to some function of the igraph library and doParallel package.
##### [igraph](http://igraph.org/r/)
A collection of network analysis tools including shortest-path analysis.   
##### [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)
Foreach Parallel adaptor for the 'parallel' package.

The C++ version is optimize to use OpenMP for dispatching tasks to different threads.
##### [OpenMP](http://openmp.org/wp/)  

