# msPac
This package was written based off of the MS software by Richard Hudson. It was written with the mentorship of Dr. Allan Strand. 

In order to run msPac, start an R Script, load msPac and setup this layout:<br />

library('msPac') <br />

nsam = 15 <br />
nreps = 1 <br />
t = 32 <br />
I = c(0) <br />
migration = c(0) <br />
en_matrix = matrix(c(0), ncol = 1, byrow = TRUE) <br />
ej_matrix = matrix(c(0), ncol = 1, byrow = TRUE) <br />
es_matrix = matrix(c(0), ncol = 1, byrow = TRUE) <br />
migr = 0 <br />
tree = ms_main(nsam, nreps, t, I, migr, migration, en_matrix, ej_matrix) <br />


Replace the zeroes and numbers with the ms call you would like to run. Migration will be a vector; the en, ej and es matrices need to be by row, as indicated in the matrix formation (see below for an example with everything filled out). The call will return the genetic tree. <br />

nsam = 40 <br />
nreps = 1 <br />
t = 32 <br />
I = c(4, 10, 10, 10, 10) <br />
migration = c(0.0, 1213.0, 1213.0,  986.1, 1213.0, 0.0, 986.1, 1213.0, 1213.0,  986.1, 0.0, 1213.0, 986.1, 1213.0, 1213.0, 0.0) <br />
en_matrix = matrix(c(0,1,1, 0, 2, 1, 0, 3, 1, 0, 4, 1), ncol = 3, byrow = TRUE) <br />
ej_matrix = matrix(c(0.00035, 2,1, .0006, 3, 2, .00085, 4, 3), ncol = 3, byrow = TRUE) <br />
es_matrix = matrix(c(.05, 1, 2), ncol = 3, byrow = TRUE) <br />
tree = ms_main(nsam, nreps, t, I, migr, migration, en_matrix, ej_matrix) <br />
