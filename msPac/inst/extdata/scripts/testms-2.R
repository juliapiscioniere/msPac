###
### trying out a project the size of holosim
###

library(tictoc)
library(rmetasim)
library(phyclust)
library(msPac)
source("pops2ms.R")

rows= 5
cols= rows

popn=rows*cols

pops=data.frame(pop=1:popn,row=1:popn,col=1,arrive=popn:1,source=0:(popn-1))

migmat=landscape.mig.matrix(h.dim=c(rows,cols),h=popn,s=2,
                       mig.model="distance",distance.fun=dweibull,
                       shape=1,scale=2)$R.int

mig.entries = which(migmat>0,arr.ind=T)

##hack for testing
texp=1.1 * max(pops$arrive)
nloci = 1
foundN = 1
if (dim(pops)[1]<100) popsamps=rep(10,dim(pops)[1]) else popsamps = sample(c(2,2),dim(pops)[1],T,prob=c(0.1,0.90))
popsizes = rep(1000,length(popsamps))
mutrate=1e-8
frag.length=80
island=NULL #rewrite new pops function if going to change island

#tic("sleeping")
new_tree = pops2ms_new(pops,
                     migmat,
                     texp=texp,
                     nloci = nloci,
                     foundN = foundN,
                     popsizes = popsizes,
                     popsamps = popsamps,
                     mutrate = mutrate,
                     frag.length=frag.length,
                     island=island)
output = seqgen(opts = "-mHKY -l40 -q", newick.tree = new_tree) #-u to increase tree
outFile = file("seqGenOutput_Pop25_Popsamps2.txt")
write(output_seqgen, outFile)
close(outFile)
#toc()