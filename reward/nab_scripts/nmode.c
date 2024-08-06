#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabc.h"
#include "sff.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;
FILE* nabout;

static MOLECULE_T *m;


static REAL_T x[4000], fret;


int main( argc, argv )
	int	argc;
	char	*argv[];
{
       nabout = stdout; /*default*/

	mytaskid=0; numtasks=1;
static INT_T __it0001__;
static INT_T __it0002__;
static REAL_T __ft0001__;
static REAL_T __ft0002__;
m = getpdb( "/home/spine/DProjects/Dprograms_and_tests/Damber_gaff/gaff_script_nab/minimise/cb7/cb7_min.pdb", NULL );
readparm( m, "/home/spine/DProjects/Dprograms_and_tests/Damber_gaff/gaff_script_nab/preprocessing/cb7/cb7.prmtop" );
mm_options( "ntpr=50, nsnb=99999, gb=1, gbsa=1, cut=999.0, rgbmax=999.0, diel=C, dielc=1.0" );

mme_init( m, NULL, "::Z", x, NULL );
setxyz_from_mol(  &m, NULL, x );


conjgrad( x, ITEMP( __it0001__, 3 *  *( NAB_mri( m, "natoms" ) ) ),  &fret, mme, FTEMP( __ft0001__, 1.000000E-04 ), FTEMP( __ft0002__, 1.000000E-04 ), ITEMP( __it0002__, 200000 ) );


mm_options( "ntpr=1" );
newton( x, ITEMP( __it0001__, 3 *  *( NAB_mri( m, "natoms" ) ) ),  &fret, mme, mme2, FTEMP( __ft0001__, 1.000000E-08 ), FTEMP( __ft0002__, 0.000000E+00 ), ITEMP( __it0002__, 60 ) );


nmode( x, 3 *  *( NAB_mri( m, "natoms" ) ), mme2, 0, 0, 0.000000E+00, 0.000000E+00, 0 );

	exit( 0 );
}
