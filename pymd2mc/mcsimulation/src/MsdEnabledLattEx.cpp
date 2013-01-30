#include "MsdEnabledLattEx.h"

lattIndex MsdEnabledLattEx::mLastPos( 0 );
lattIndex MsdEnabledLattEx::mLastValue( 0 );
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


MsdEnabledLattEx::MsdEnabledLattEx( TriangularLattice* latt, lattIndex *proteins, lattIndex proteinsCnt ) 
    : LattExchanger( latt, proteins, proteinsCnt ) 
{
    mTracking = new lattIndex[ latt->getLatticeSize() ];
    mWasProtein = new bool[ latt->getLatticeSize() ];
    for ( lattIndex i(0) ; i < latt->getLatticeSize() ; ++i )
    {
        mTracking[i] = i;
    }
    mPBCCorrection = new vectorDist[ latt->getLatticeSize() ];
    
    for ( lattIndex i(0) ; i < mpLatt->getLatticeSize() ; ++i )
    {
        mWasProtein[i] = false;
    }
}

void MsdEnabledLattEx::updatePBCCorrection( lattIndex site, lattIndex newSite )
{
    vectorDist localDist( calcDist( site, newSite ) );
    incDist( localDist );
    mPBCCorrection[ mTracking[ site ] ] += localDist;
}

void MsdEnabledLattEx::exchangeSites( lattIndex pos1, lattIndex pos2 ) 
{
    mTracking[pos1] ^= mTracking[pos2];
    mTracking[pos2] ^= mTracking[pos1];
    mTracking[pos1] ^= mTracking[pos2];
    LattExchanger::exchangeSites( pos1, pos2 );
    if ( !isNotPBCJump( pos1, pos2 ) )
    {
        updatePBCCorrection( pos1, pos2 );
        updatePBCCorrection( pos2, pos1 );
    }
}

void MsdEnabledLattEx::moveProtein( lattIndex site, lattIndex destination )
{
    LattExchanger::moveProtein( site, destination );
}

void MsdEnabledLattEx::moveProteinLeft( lattIndex site )
{
    LattExchanger::moveProteinLeft( site );
}
void MsdEnabledLattEx::moveProteinRight( lattIndex site )
{
    LattExchanger::moveProteinRight( site );
}
void MsdEnabledLattEx::moveProteinUp( lattIndex site )
{
    LattExchanger::moveProteinUp( site );
}
void MsdEnabledLattEx::moveProteinDown( lattIndex site )
{
    LattExchanger::moveProteinDown( site );
}

void MsdEnabledLattEx::initPushAndPop( int site )
{
   mLastPos = mpLatt->cutToPos( site );
   mLastValue = mTracking[ mLastPos ];
}
void MsdEnabledLattEx::pushAndPop( lattIndex site, int &value )
{
    site = mpLatt->cutToPos( site );
    LattExchanger::pushAndPop( site, value );
    lattIndex temp( mTracking[ site ] );
    mTracking[ site ] = mLastValue;
    mLastValue = temp;
    if ( !isNotPBCJump( mLastPos, site ) )
    {
        updatePBCCorrection( site, mLastPos );
    }
    mLastPos = site;
}
void MsdEnabledLattEx::push( lattIndex site, int & value )
{
    site = mpLatt->cutToPos( site );
    LattExchanger::push( site, value );
    mTracking[ site ] = mLastValue;
    if ( !isNotPBCJump( mLastPos, site ) )
    {
        updatePBCCorrection( site, mLastPos );
    }
    mLastPos = 0;
    mLastValue = 0;
}

void MsdEnabledLattEx::calcStat( double & msd, double & protMsd )
{
    unsigned int count( 0 ), countProt( 0 );
    msd = 0;
    protMsd = 0;
    for ( lattIndex i = 0 ; i < mpLatt->getLatticeSize() ; ++i )
    {
        if (  !mWasProtein[ mTracking[i] ] )  // lipids
        {
            msd += ( calcDist( i, mTracking[i] ) - mPBCCorrection[ mTracking[i] ] ).squareDisp();
            ++count;
        }
        else //proteins
        {
            protMsd += ( calcDist( i, mTracking[i] ) - mPBCCorrection[ mTracking[i] ] ).squareDisp();
            ++countProt;
        }
    }
    msd /= count;
    protMsd /= ( countProt );
}



MsdEnabledLattEx::~MsdEnabledLattEx()
{
    delete[] mPBCCorrection;
    delete[] mTracking;
    delete[] mWasProtein;
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

