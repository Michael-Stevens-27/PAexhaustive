library(parallel)
library(compiler)
library(microbenchmark)
library(ggplot2)
library(gtools)
library(gRbase)

#ncores <- detectCores()
#cl <- makeCluster(ncores - 1)
#a <- mapply(sum,x,y)
#b <- mcmapply(sum,x,y, mc.cores = 3)
#stopCluster(cl)


AP_allocation <- function(grid_x, grid_y)
                {
                allocation <- cbind(rep(1:((grid_x + 1)*(grid_y + 1)),each=((grid_x + 1)*(grid_y + 1))),
                      rep(1:((grid_x + 1)*(grid_y + 1)),(grid_x + 1)*(grid_y + 1)))
                allocation <- unique(t(apply(allocation, 1, sort)))
                return(allocation)
                }

#compiled_allocation <- cmpfun(AP_allocation)

#compare <- microbenchmark(AP_allocation(4,3), combn(20, 2), combinations(20, 2), times = 1000)
autoplot(compare)




###################################################################

source("http://bioconductor.org/biocLite.R")
biocLite("gRbase")
a <- t(combnPrim(10000,2))
