/*  Link in this file for random number generation using drand48() */

#include <stdio.h>
#include <stdlib.h>

         double
ran1()
{
        double drand48();
        return( drand48() );
}               

int
commandlineseed( int *seeds) /* takes seeds from random number generator in testmain to use in the program*/
{
    unsigned short seedv[3], *seed48();
    seedv[0] = seeds[0];
    seedv[1] = seeds[1];
    seedv[2] = seeds[2];
	printf("\nseeds: %d %d %d\n", seedv[0], seedv[1], seedv[2] );
	seed48(seedv);
    //printf("\nran1: %f\n", ran1());
	return(3); /* returns the number of seeds */
}

