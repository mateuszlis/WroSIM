#include "MsdEnabledLattEx.h"

bool const MsdEnabledLattEx::isNeighbor[3][3] = {
    { 0, 1, 1 },
    { 1, 1, 1 },
    { 1, 1, 0 }
};

ostream &operator<<( ostream &stream, vectorDist & dist )
{
    stream << "row=" << dist.row << " col=" << dist.col << ". ";
    return stream;
}


MsdEnabledLattEx::MsdEnabledLattEx( TriangularLattice* latt ) 
    : LattExchanger( latt ) 
{
    mTracking = new lattIndex[ latt->getLatticeSize() ];
    for ( lattIndex i = 0 ; i < latt->getLatticeSize() ; ++i )
    {
        mTracking[i] = i;
    }
    mPBCCorrection = new vectorDist[ latt->getLatticeSize() ];
}

void MsdEnabledLattEx::exchangeSites( lattIndex pos1, lattIndex pos2 ) 
{
    mTracking[pos1] ^= mTracking[pos2];
    mTracking[pos2] ^= mTracking[pos1];
    mTracking[pos1] ^= mTracking[pos2];
    LattExchanger::exchangeSites( pos1, pos2 );
    if ( !isNotPBCJump( pos1, pos2 ) )
    {
        vectorDist lDist( calcDist( pos1, pos2 ) ); // local dist
        incDist( lDist );
        mPBCCorrection[ mTracking[ pos1 ] ] += lDist;
        lDist = calcDist( pos2, pos1 );
        incDist( lDist );
        mPBCCorrection[ mTracking[ pos2 ] ] += lDist;
    }
}

double MsdEnabledLattEx::calcStat()
{
    double msd( 0 );
    for ( lattIndex i = 0 ; i < mpLatt->getLatticeSize() ; ++i )
    {
        msd += ( calcDist( i, mTracking[i] ) - mPBCCorrection[ mTracking[i] ] ).squareDisp();
    }
    msd /= mpLatt->getLatticeSize();
    return msd;
}


vectorDist MsdEnabledLattEx::calcDist( lattIndex pos1, lattIndex pos2 )
{
        lattIndex startRow = pos1 / mpLatt->getRowSize();
        lattIndex startCol = pos1 % mpLatt->getRowSize();
        lattIndex endRow = pos2 / mpLatt->getRowSize();
        lattIndex endCol = pos2 % mpLatt->getRowSize();
        return vectorDist( startRow - endRow, startCol - endCol );

}


MsdEnabledLattEx::~MsdEnabledLattEx()
{
    delete[] mPBCCorrection;
    delete[] mTracking;
}


// private member functions implementation
//

bool MsdEnabledLattEx::isNotPBCJump( lattIndex pos1, lattIndex pos2 ) 
{
    vectorDist dist( calcDist( pos1, pos2 ) );
    if ( dist.squareDisp() <= 2 )
    {
        return isNeighbor[ dist.col - 1 ][ dist.row - 1 ];
    }
    return false;
            
}

void MsdEnabledLattEx::incDist( vectorDist & pbcDist )
{
    if ( pbcDist.col > 1 )
    {
        pbcDist.col++;
    }
    if ( pbcDist.col < -1 )
    {
        pbcDist.col--;
    }
    if ( pbcDist.row > 1 )
    {
        pbcDist.row++;
    }
    if ( pbcDist.row < -1 )
    {
        pbcDist.row--;
    }
}

