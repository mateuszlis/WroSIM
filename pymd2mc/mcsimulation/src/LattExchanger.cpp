#include "TriangularLattice.h"
#include "LattExchanger.h"

void LattExchanger::exchangeSites( lattIndex pos1, lattIndex pos2 ) const
{
    // fast variable values exchange trick (pos1 != pos2)
	mpLatt->getLattice()[pos1] ^= mpLatt->getLattice()[pos2];
    mpLatt->getLattice()[pos2] ^= mpLatt->getLattice()[pos1];
    mpLatt->getLattice()[pos1] ^= mpLatt->getLattice()[pos2];
}

