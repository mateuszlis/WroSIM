#include "TriangularLattice.h"
#include "LattExchanger.h"

void LattExchanger::exchangeSites( lattIndex pos1, lattIndex pos2 )
{
    // fast variable values exchange trick (pos1 != pos2)
	mpLatt->getLattice()[pos1] ^= mpLatt->getLattice()[pos2];
    mpLatt->getLattice()[pos2] ^= mpLatt->getLattice()[pos1];
    mpLatt->getLattice()[pos1] ^= mpLatt->getLattice()[pos2];
}

vectorDist LattExchanger::calcDist( lattIndex pos1, lattIndex pos2 )
{
    lattIndex startRow = pos1 / mpLatt->getRowSize();
    lattIndex startCol = pos1 % mpLatt->getRowSize();
    lattIndex endRow = pos2 / mpLatt->getRowSize();
    lattIndex endCol = pos2 % mpLatt->getRowSize();
    return vectorDist( startRow - endRow, startCol - endCol );

}

