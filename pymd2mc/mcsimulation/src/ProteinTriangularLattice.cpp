#include "ProteinTriangularLattice.h"

#include <limits>

ProteinTriangularLattice::ProteinTriangularLattice( lattIndex latticeSize
                                , lattIndex rowSize
                                , lattIndex firstTypeParticlesCnt
                                , lattIndex proteinCnt
                                , bool distributeRandomly )
    : TriangularLattice( latticeSize, rowSize, firstTypeParticlesCnt, distributeRandomly )
      , mProteinCnt( proteinCnt )
{
    mProteinSites[ 0 ] = 0;
    mProteinSites[ 1 ] = 1;
    mProteinSites[ 2 ] = -1;
    mProteinSites[ 3 ] = rowSize;
    mProteinSites[ 4 ] = -rowSize;
    mProteinSites[ 5 ] = rowSize - 1;
    mProteinSites[ 6 ] = -rowSize + 1; // TODO: this fixed the protein size

    if( distributeRandomly )
        distributeProteinsRandomly();
    else
        distributeProteins();
}

void ProteinTriangularLattice::distributeProteinsRandomly()
{
    if ( mProteinCnt )
    {
        mProteins = new lattIndex[ mProteinCnt ];
        for ( int i = 0; i < mProteinCnt; ++i )
        {
            lattIndex pos = rand() % getLatticeSize(); //FIXME: if RAND_MAX is too small this causes bad distribution
            while ( !isFree( pos ) )
            {
                    pos = rand() % getLatticeSize();
            };
            mProteins[i] = pos;
            putProtein( pos );
        }
    }
}

void ProteinTriangularLattice::distributeProteins()
{
    if ( mProteinCnt )
    {
        mProteins = new lattIndex[ mProteinCnt ];
        int counter( 0 );
        for ( lattIndex pos = 0 ; pos < getLatticeSize() && counter < mProteinCnt ; ++pos )
        {
            if ( isFree( pos ) )
            {
                putProtein( pos );
                mProteins[ counter ] = pos;
                ++counter;
            }
        }
    }
}

void ProteinTriangularLattice::moveProtein( lattIndex site, lattIndex destination )
{
    double res = ( rand() ) / static_cast< float >( RAND_MAX );
    if ( res < 0.25 ) moveProteinRight( site );
    else if ( res < 0.5 ) moveProteinLeft( site );
    else if ( res < 0.75 ) moveProteinUp( site );
    else moveProteinDown( site );
}
void ProteinTriangularLattice::moveProteinRight( lattIndex site )
{
    static const int RIGHT_NEIGH_OFFSET = 1;
    static const int RIGHT_TOP_NEIGH_OFFSET = -mRowSize + 1;
    static const int RIGHT_BOTTOM_NEIGH_OFFSET = mRowSize ;
    static const int LEFT_BOTTOM_NEIGH_OFFSET = mRowSize - 1 ;
    static const int LEFT_TOP_NEIGH_OFFSET = -mRowSize;
    static const int LEFT_NEIGH_OFFSET = -1;
    static const int MOVED_SITES_SIZE = 3;


    lattMember movedSites[ MOVED_SITES_SIZE ];
    lattMember protein( mpLattice[ site ] );
    lattIndex rightSite = site + 1; // rightest site of the protein
    
    movedSites[0] = get( rightSite + RIGHT_NEIGH_OFFSET );
    movedSites[1] = get( rightSite + RIGHT_TOP_NEIGH_OFFSET );
    movedSites[2] = get( rightSite + RIGHT_BOTTOM_NEIGH_OFFSET );

    if ( !isSpaceToMove( movedSites, MOVED_SITES_SIZE ) )
    {
        return;
    }
    // move right top site
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + LEFT_TOP_NEIGH_OFFSET, movedSites[1] );


    // move right bottom site
    pushAndPop( site + RIGHT_BOTTOM_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + RIGHT_BOTTOM_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[2] );

    // move center
    lattIndex lastMovedSite = rightSite + RIGHT_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_TOP_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_TOP_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_BOTTOM_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_BOTTOM_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_BOTTOM_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );

    // move the protein
    lattIndex newProteinCenter( cutToPos( site + RIGHT_NEIGH_OFFSET ) );
    mpLattice[ newProteinCenter ] = protein;
    for ( int i = 0 ; i < mNeighbCnt ; ++i )
    {
        mpLattice[ getNeighbIndex( newProteinCenter, i ) ] = protein;
    }
    // update array
    //FIXME: this might be efficiency drawback
    updateProteinArray( site, newProteinCenter );

}


void ProteinTriangularLattice::moveProteinLeft( lattIndex site )
{
    static const int RIGHT_NEIGH_OFFSET = 1;
    static const int RIGHT_TOP_NEIGH_OFFSET = -mRowSize + 1;
    static const int RIGHT_BOTTOM_NEIGH_OFFSET = mRowSize ;
    static const int LEFT_BOTTOM_NEIGH_OFFSET = mRowSize - 1 ;
    static const int LEFT_TOP_NEIGH_OFFSET = -mRowSize;
    static const int LEFT_NEIGH_OFFSET = -1;

    static const int MOVED_SITES_SIZE = 3;
    lattMember movedSites[ MOVED_SITES_SIZE ];
    lattMember protein( mpLattice[ site ] );
    lattIndex leftSite = site - 1; // rightest site of the protein
    
    movedSites[0] = get( leftSite + LEFT_NEIGH_OFFSET );
    movedSites[1] = get( leftSite + LEFT_TOP_NEIGH_OFFSET );
    movedSites[2] = get( leftSite + LEFT_BOTTOM_NEIGH_OFFSET );

    if ( !isSpaceToMove( movedSites, MOVED_SITES_SIZE ) )
    { 
        return;
    }

    // move left top site
    pushAndPop( site + LEFT_TOP_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + LEFT_TOP_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET, movedSites[1] );


    // move left bottom site
    pushAndPop( site + LEFT_BOTTOM_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + LEFT_BOTTOM_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + RIGHT_BOTTOM_NEIGH_OFFSET, movedSites[2] );

    // move center
    lattIndex lastMovedSite = leftSite + LEFT_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_TOP_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_TOP_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_BOTTOM_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_BOTTOM_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_BOTTOM_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );

    // move the protein
    lattIndex newProteinCenter( cutToPos( site + LEFT_NEIGH_OFFSET ) );
    mpLattice[ newProteinCenter ] = protein;
    for ( int i = 0 ; i < mNeighbCnt ; ++i )
    {
        mpLattice[ getNeighbIndex( newProteinCenter, i ) ] = protein;
    }
    // update array
    updateProteinArray( site, newProteinCenter );

}

void ProteinTriangularLattice::moveProteinUp( lattIndex site )
{
    static const int RIGHT_NEIGH_OFFSET = 1;
    static const int RIGHT_TOP_NEIGH_OFFSET = -mRowSize + 1;
    static const int RIGHT_BOTTOM_NEIGH_OFFSET = mRowSize ;
    static const int LEFT_BOTTOM_NEIGH_OFFSET = mRowSize - 1 ;
    static const int LEFT_TOP_NEIGH_OFFSET = -mRowSize;
    static const int LEFT_NEIGH_OFFSET = -1;

    static const int MOVED_SITES_SIZE = 3;
    lattMember movedSites[ MOVED_SITES_SIZE ];
    lattMember protein( mpLattice[ site ] );
    lattIndex topSite = site + LEFT_TOP_NEIGH_OFFSET; // topest site of the protein
    
    movedSites[0] = get( topSite + LEFT_TOP_NEIGH_OFFSET );
    movedSites[1] = get( topSite + LEFT_NEIGH_OFFSET );
    movedSites[2] = get( topSite + RIGHT_TOP_NEIGH_OFFSET );

    if ( !isSpaceToMove( movedSites, MOVED_SITES_SIZE ) )
    { 
        return;
    }

    // move left top site 
    pushAndPop( site + LEFT_NEIGH_OFFSET + LEFT_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + LEFT_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[1] );


    // move left bottom site
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + RIGHT_NEIGH_OFFSET, movedSites[2] );

    // move center
    lattIndex lastMovedSite = topSite + LEFT_TOP_NEIGH_OFFSET + LEFT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_BOTTOM_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_BOTTOM_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_BOTTOM_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_BOTTOM_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_TOP_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );

    // move the protein
    lattIndex newProteinCenter( cutToPos( site + LEFT_TOP_NEIGH_OFFSET ) );
    mpLattice[ newProteinCenter ] = protein;
    for ( int i = 0 ; i < mNeighbCnt ; ++i )
    {
        mpLattice[ getNeighbIndex( newProteinCenter, i ) ] = protein;
    }
    // update array
    updateProteinArray( site, newProteinCenter );

}

void ProteinTriangularLattice::moveProteinDown( lattIndex site )
{
    static const int RIGHT_NEIGH_OFFSET = 1;
    static const int RIGHT_TOP_NEIGH_OFFSET = -mRowSize + 1;
    static const int RIGHT_BOTTOM_NEIGH_OFFSET = mRowSize ;
    static const int LEFT_BOTTOM_NEIGH_OFFSET = mRowSize - 1 ;
    static const int LEFT_TOP_NEIGH_OFFSET = -mRowSize;
    static const int LEFT_NEIGH_OFFSET = -1;

    static const int MOVED_SITES_SIZE = 3;
    lattMember movedSites[ MOVED_SITES_SIZE ];
    lattMember protein( mpLattice[ site ] );
    lattIndex bottomSite = site + RIGHT_BOTTOM_NEIGH_OFFSET; // bottomest site of the protein
    
    movedSites[0] = get( bottomSite + RIGHT_BOTTOM_NEIGH_OFFSET );
    movedSites[1] = get( bottomSite + RIGHT_NEIGH_OFFSET );
    movedSites[2] = get( bottomSite + LEFT_BOTTOM_NEIGH_OFFSET );

    if ( !isSpaceToMove( movedSites, MOVED_SITES_SIZE ) )
    { 
        return;
    }

    // move left bottom site 
    pushAndPop( site + RIGHT_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + RIGHT_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET, movedSites[1] );


    // move left bottom site
    pushAndPop( site + LEFT_BOTTOM_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + LEFT_TOP_NEIGH_OFFSET + LEFT_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + LEFT_NEIGH_OFFSET, movedSites[2] );

    // move center
    lattIndex lastMovedSite = bottomSite + RIGHT_BOTTOM_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_TOP_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_TOP_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_TOP_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_TOP_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_BOTTOM_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );

    // move the protein
    lattIndex newProteinCenter( cutToPos( site + RIGHT_BOTTOM_NEIGH_OFFSET ) );
    mpLattice[ newProteinCenter ] = protein;
    for ( int i = 0 ; i < mNeighbCnt ; ++i )
    {
        mpLattice[ getNeighbIndex( newProteinCenter, i ) ] = protein;
    }
    // update array
    updateProteinArray( site, newProteinCenter );

}

void ProteinTriangularLattice::updateProteinArray( lattIndex site, lattIndex newPos )
{
    for ( int i = 0 ; i < mProteinCnt ; ++i )
    {
        if ( mProteins[i] == site )
        {
            mProteins[i] = newPos;
        }
    }
}


ProteinTriangularLattice::~ProteinTriangularLattice()
{
    if ( mProteinCnt )
        delete mProteins;
}


// helpers
//

bool ProteinTriangularLattice::isFree( lattIndex pos )
{
    bool isFree = true;
    for ( unsigned int i = 0 ; i < mProteinSize ; ++i )
    {
        isFree = isFree && ( mpLattice[ cutToPos( pos + mProteinSites[i] ) ] == LIPID_A );
    }
    return isFree;
};

void ProteinTriangularLattice::pushAndPop( int site, lattMember &value )
{
    site = cutToPos( site );
    lattMember temp = value;
    value = mpLattice[ site ];
    mpLattice[ site ] = temp;
}

void ProteinTriangularLattice::putProtein( lattIndex pos, LATTICE_FIELD_NAMES protein )
{
    for ( unsigned int i = 0 ; i < mProteinSize ; ++i )
    {
        mpLattice[ cutToPos( pos + mProteinSites[i] ) ] = protein;
    }
};

bool ProteinTriangularLattice::isSpaceToMove( lattMember movedSites[], lattIndex movedSitesSize )
{
    for ( int i = 0 ; i < movedSitesSize ; ++i )
    {
        if ( movedSites[i] == PROTEIN_A || movedSites[i] == PROTEIN_B )
            return false;
    }
    return true;
}

