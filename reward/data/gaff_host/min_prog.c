#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;

static MOLECULE_T *m;

static REAL_T *x,  *f, fret;

int main( argc, argv )
	int	argc;
	char	*argv[];
{
	nabout = stdout; /*default*/

	mytaskid=0; numtasks=1;
static INT_T __gdab0001__;
static INT_T __gdab0002__;
m = getpdb( "cb7_min.pdb", NULL );
readparm( m, "cb7.prmtop" );
__gdab0001__ = 3 *  *( NAB_mri( m, "natoms" ) );DA_ALLOC( x = ( REAL_T * )malloc( __gdab0001__ * ( sizeof( REAL_T ) ) ), "main", "x" );__gdab0002__ = 3 *  *( NAB_mri( m, "natoms" ) );DA_ALLOC( f = ( REAL_T * )malloc( __gdab0002__ * ( sizeof( REAL_T ) ) ), "main", "f" );setxyz_from_mol(  &m, NULL, x );
mm_options( "ntpr=5000, nsnb=99999, gb=1, gbsa=1, cut=99.0, rgbmax=99.0, diel=C" );
mme_init( m, NULL, ":::", x, NULL );



mm_options( "ntpr=100" );




nmode( x, 3 *  *( NAB_mri( m, "natoms" ) ), mme2, 0, 0, 0.000000E+00, 0.000000E+00, 0 );


	exit( 0 );
}
