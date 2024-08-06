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
static INT_T __it0001__;
static INT_T __it0002__;
static REAL_T __ft0001__;
static REAL_T __ft0002__;
m = getpdb( "/home/spine/DProjects/DMCTS/ChemTSv2/reward/data/gaff_host/cb7_4amber.pdb", NULL );
readparm( m, "/home/spine/DProjects/DMCTS/ChemTSv2/reward/data/gaff_host/cb7.prmtop" );
__gdab0001__ = 3 *  *( NAB_mri( m, "natoms" ) );DA_ALLOC( x = ( REAL_T * )malloc( __gdab0001__ * ( sizeof( REAL_T ) ) ), "main", "x" );__gdab0002__ = 3 *  *( NAB_mri( m, "natoms" ) );DA_ALLOC( f = ( REAL_T * )malloc( __gdab0002__ * ( sizeof( REAL_T ) ) ), "main", "f" );setxyz_from_mol(  &m, NULL, x );
mm_options( "ntpr=50, nsnb=99999, gb=1, gbsa=1, cut=999.0, rgbmax=999.0, diel=C, dielc=1.0" );
mme_init( m, NULL, "::Z", x, NULL );

conjgrad( x, ITEMP( __it0001__, 3 *  *( NAB_mri( m, "natoms" ) ) ),  &fret, mme, FTEMP( __ft0001__, 1.000000E-05 ), FTEMP( __ft0002__, 1.000000E-02 ), ITEMP( __it0002__, 100000 ) );

mm_options( "ntpr=1" );
newton( x, ITEMP( __it0001__, 3 *  *( NAB_mri( m, "natoms" ) ) ),  &fret, mme, mme2, FTEMP( __ft0001__, 1.000000E-07 ), FTEMP( __ft0002__, 0.000000E+00 ), ITEMP( __it0002__, 100 ) );

putpdb( "/home/spine/DProjects/DMCTS/ChemTSv2/reward/data/gaff_host/cb7_min.pdb", m, NULL );

nmode( x, 3 *  *( NAB_mri( m, "natoms" ) ), mme2, 0, 0, 0.000000E+00, 0.000000E+00, 0 );


	exit( 0 );
}
