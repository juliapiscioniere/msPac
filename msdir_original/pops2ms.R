
# takes a pops object and a bunch of other parameters and writes an ms command line'
#ms line looks something like this...
#ms #sam #reps(loci) -t THETA 
#-I #pops #pop1 #pop2 ... #popN 
#-ma #migmat_line1 #migmat_line2 ... #migmat_lineN
#-en time popID size_scalar (new size is x*N0)
pops2ms = function(pops, migmat, texp, nloci, foundN = 1, popsizes = 1000, popsamps = 100, mutrate = 1e-8, frag.length = 80, island=NULL) 
{
  if(length(popsizes) == 1) {
		popsizevec = rep(popsizes, dim(pops)[1]) 
	} else {
		popsizevec = popsizes
	}

	if(length(popsamps) == 1) {
		popsampvec = rep(popsamps, dim(pops)[1]) 
	} else {
		popsampvec = popsamps
	}

	N0 = max(popsizevec)
	THETA = 4*N0*mutrate*frag.length
	NPOP = length(pops[,1])
	
if (is.null(island))
{                                        #MIGMAT
	migvec = signif(4*N0*as.vector(migmat), 4);  #this puts units in ms units? and 4 sig digits
	migvec = paste("-ma",paste( migvec, sep = " ", collapse = " ") )
} else migvec = paste(island)
        
	#Set all populations to their current sizes!!
	relpopsize = popsizevec/N0;
	sizeline = paste("-en 0", c(1:NPOP), relpopsize)
	sizeline = paste(sizeline, sep = " ", collapse = " ")

	#Now set the growth rates of the populations
	col_times = pops$arrive/(4*N0)
	col_times[1] = texp/(4*N0)
	g_rates = signif(-1*log(foundN/popsizevec)/col_times,4)
	growline = paste("-eg 0", c(1:NPOP), g_rates) 
	growline_new = paste(c(1:NPOP), g_rates) 
	growline_newms = paste(0, c(1:NPOP), g_rates)
	growline = paste(growline, sep = " ", collapse = " ")

	#Now the colonization history
	sorted.pops = pops[order(pops$arrive, decreasing = TRUE),]
	sorted.events = sorted.pops[(sorted.pops$arrive > 0)&(sorted.pops$source!=0),]
	times = texp - sorted.events$arrive
	times = signif(times/(4*N0),4) #4 sig digits
	jumpline = paste("-ej", times, sorted.events$pop, sorted.events$source)
	jumpline = paste(jumpline, sep = " ", collapse = " ")

  progline =  paste(sum(popsampvec), nloci)

	parline = paste("-t", THETA*10000, "-I", NPOP, paste(popsampvec, sep = " ", collapse = " "),  migvec, sizeline, jumpline,"-seeds 5004 25353 42948 ","-T")
  return(list(progline,parline))  
}

pops2ms_new = function(nsamp, nrep, pops, migmat, texp, nloci, foundN = 1, popsizes = 1000, popsamps = 100, mutrate = 1e-8, frag.length = 80, island=NULL) 
{
  if(length(popsizes) == 1) {
    popsizevec = rep(popsizes, dim(pops)[1]) 
  } else {
    popsizevec = popsizes
  }
  
  if(length(popsamps) == 1) {
    popsampvec = rep(popsamps, dim(pops)[1]) 
  } else {
    popsampvec = popsamps
  }
  
  N0 = max(popsizevec)
  THETA = 4*N0*mutrate*frag.length
  NPOP = length(pops[,1])
  
  #Set all populations to their current sizes!!
  relpopsize = popsizevec/N0;
  sizeline_matrix = matrix(nrow = NPOP, ncol = 3) #make sure this is npop
  for (i in 1:nrow(sizeline_matrix)){
    sizeline_matrix[i,] = c(0,i, relpopsize[i])
  }
  
  #Now the colonization history
  sorted.pops = pops[order(pops$arrive, decreasing = TRUE),]
  sorted.events = sorted.pops[(sorted.pops$arrive > 0)&(sorted.pops$source!=0),]
  times = texp - sorted.events$arrive
  times = signif(times/(4*N0),4) #4 sig digits
  jumpline_matrix = matrix(ncol = 3, nrow = length(times))
  
  for (i in 1:nrow(jumpline_matrix)){
    jumpline_matrix[i,] = c(times[i],sorted.events$pop[i],sorted.events$source[i])
  }
  
  progline =  paste(sum(popsampvec), nloci)
  
  t = THETA*10000
  
  popsampvec = as.integer(popsampvec)
  I = c(NPOP,popsampvec) # I for new ms

  migmat = signif(4*N0*as.vector(migmat), 4); #migmat for new ms - check the island if though!
  nsamp = as.numeric(nsamp)
  nrep = as.numeric(nrep)
  en_matrix = as.matrix(sizeline_matrix)
  ej_matrix = as.matrix(jumpline_matrix)
  
  new_tree = call_ms(nsamp, nrep, t, I, migmat, en_matrix, ej_matrix)
  return(new_tree)
  
}

call_ms = function(nsamp, nrep, t, I, migmat, en_matrix, ej_matrix){
  setwd("~/Desktop/Research/msPac_final/msPac")
  seeds = c(5004, 25353, 42948)
  migr = 0
  variable_list = c(1, 4, 0, 0, 1, 0, 1, 0, 0, 1, 2)
  devtools::load_all()
  new_tree = ms_main(nsamp, nrep, t, variable_list, I, migr, migmat, en_matrix, ej_matrix, seeds)
  return(new_tree)
}

