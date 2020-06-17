# msPac
In order to run msPac, start an R Script with this layout and set working directory to msPac (within this repository):<br />

nsam = 15 <br />
nreps = 1 <br />
t = 32 <br />
I = c(0) <br />
migration = c(0) <br />
en_matrix = matrix(c(0), ncol = 1, byrow = TRUE) <br />
ej_matrix = matrix(c(0), ncol = 1, byrow = TRUE) <br />

variable_list = c(1, 4, 0, 0, 1, 0, 1, 0, 0, 1, 2) <br />
migr = 0 <br />
seeds = c(5004, 25353, 42948) <br />
devtools::load_all() <br />
tree = ms_main(nsam, nreps, t, variable_list, I, migr, migration, en_matrix, ej_matrix, seeds) <br />


Replace the zeroes and numbers with the ms call you would like to run. Migration will be a vector; the en and ej matrices need to be by row, as indicated in the matrix formation (see below for an example with everything filled out). The variable list can be changed as well, these are default settings. The seeds can also be changed. The call will return the genetic tree. <br />

nsam = 40 <br />
nreps = 1 <br />
t = 32 <br />
I = c(4, 10, 10, 10, 10) <br />
migration = c(0.0, 1213.0, 1213.0,  986.1, 1213.0, 0.0, 986.1, 1213.0, 1213.0,  986.1, 0.0, 1213.0, 986.1, 1213.0, 1213.0, 0.0) <br />
en_matrix = matrix(c(0,1,1, 0, 2, 1, 0, 3, 1, 0, 4, 1), ncol = 3, byrow = TRUE) <br />
ej_matrix = matrix(c(0.00035, 2,1, .0006, 3, 2, .00085, 4, 3), ncol = 3, byrow = TRUE) <br />
