#Test the vanilla R version (in development)
source("BRIDES_2022.R")

X=load_network("sample_X.txt")
y=load_network("sample_Y.txt")

BRIDES(X,y,runmode = "stepwise", min_additional = 2, max_additional = 3)
BRIDES(X,y,runmode = "exhaustive", min_additional = 1, max_additional = 2)
BRIDES(X,y,runmode = "genetics", min_additional = 1, max_additional = 3)
BRIDES(X,y,runmode = "stepwise", min_additional = 1, max_additional = 3)

#Test the package version located in the Releases section of the Github
install.packages("BRIDES_1.2.1.tar.gz", repos = NULL, type = "source")
library(BRIDES)
set.seed(22)
g=random_network(20,10,type = "barabasi")
plot_network(g$g1,g$g2)
BRIDES(g$g1,g$g2,runmode = "genetics", min_additional = 1, max_additional = 3, maxcore=4, max_iters = 30) # Need a few step to ensure 

BRIDES(g$g1,g$g2,runmode = "stepwise", min_additional = 1, max_additional = 3, maxcore=4) #Relatively fast

BRIDES(g$g1,g$g2,runmode = "exhaustive", min_additional = 1, max_additional = 3) # Long
