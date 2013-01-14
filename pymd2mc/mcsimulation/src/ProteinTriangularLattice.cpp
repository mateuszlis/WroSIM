#include "ProteinTriangularLattice.h"

#include <limits>
#include "LattExchanger.h"

ProteinTriangularLattice::ProteinTriangularLattice( lattIndex latticeSize
                                , lattIndex rowSize
                                , lattIndex firstTypeParticlesCnt
                                , lattIndex proteinCnt
                                , bool distributeRandomly )
    : TriangularLattice( latticeSize, rowSize, firstTypeParticlesCnt, distributeRandomly )
      , mProteinCnt( proteinCnt )
      , RIGHT_NEIGH_OFFSET( 1 )
      , RIGHT_TOP_NEIGH_OFFSET( -mRowSize + 1 )
      , RIGHT_BOTTOM_NEIGH_OFFSET( mRowSize ) 
      , LEFT_BOTTOM_NEIGH_OFFSET( mRowSize - 1 )
      , LEFT_TOP_NEIGH_OFFSET( -mRowSize )
      , LEFT_NEIGH_OFFSET( -1 )
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


    // FIXME: values that are currently hardcoded - in future will be developed using optimization
    mSitesMovedMoveRight.push_back( RIGHT_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( RIGHT_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( RIGHT_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( RIGHT_TOP_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( RIGHT_TOP_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( RIGHT_BOTTOM_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( RIGHT_BOTTOM_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( RIGHT_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( mSitesMovedMoveRight.back() + LEFT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( mSitesMovedMoveRight.back() + LEFT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( mSitesMovedMoveRight.back() + LEFT_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( mSitesMovedMoveRight.back() + LEFT_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( mSitesMovedMoveRight.back() + LEFT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( mSitesMovedMoveRight.back() + LEFT_BOTTOM_NEIGH_OFFSET );
    //additional items
    mSitesMovedMoveRight.push_back( RIGHT_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( RIGHT_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveRight.push_back( RIGHT_BOTTOM_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET );

    mSitesMovedMoveLeft.push_back( LEFT_NEIGH_OFFSET + LEFT_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( LEFT_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( LEFT_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( LEFT_TOP_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( LEFT_TOP_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( LEFT_BOTTOM_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( LEFT_BOTTOM_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( LEFT_NEIGH_OFFSET + LEFT_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( mSitesMovedMoveLeft.back() + RIGHT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( mSitesMovedMoveLeft.back() + RIGHT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( mSitesMovedMoveLeft.back() + RIGHT_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( mSitesMovedMoveLeft.back() + RIGHT_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( mSitesMovedMoveLeft.back() + RIGHT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( mSitesMovedMoveLeft.back() + RIGHT_BOTTOM_NEIGH_OFFSET );
    //additional items
    mSitesMovedMoveLeft.push_back( LEFT_NEIGH_OFFSET + LEFT_NEIGH_OFFSET + LEFT_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( LEFT_NEIGH_OFFSET + LEFT_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveLeft.push_back( LEFT_BOTTOM_NEIGH_OFFSET + LEFT_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET );


    mSitesMovedMoveUp.push_back( LEFT_TOP_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( LEFT_TOP_NEIGH_OFFSET + LEFT_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( LEFT_TOP_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( LEFT_NEIGH_OFFSET + LEFT_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( LEFT_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( RIGHT_TOP_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( RIGHT_TOP_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( LEFT_TOP_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET + LEFT_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( mSitesMovedMoveUp.back() + LEFT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( mSitesMovedMoveUp.back() + LEFT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( mSitesMovedMoveUp.back() + RIGHT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( mSitesMovedMoveUp.back() + RIGHT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( mSitesMovedMoveUp.back() + RIGHT_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( mSitesMovedMoveUp.back() + RIGHT_NEIGH_OFFSET );
    //additional items
    mSitesMovedMoveUp.push_back( LEFT_TOP_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( LEFT_TOP_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveUp.push_back( RIGHT_TOP_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET );

    mSitesMovedMoveDown.push_back( RIGHT_BOTTOM_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( RIGHT_BOTTOM_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( RIGHT_BOTTOM_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( RIGHT_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( RIGHT_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( LEFT_BOTTOM_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( LEFT_TOP_NEIGH_OFFSET + LEFT_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( RIGHT_BOTTOM_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( mSitesMovedMoveDown.back() + RIGHT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( mSitesMovedMoveDown.back() + RIGHT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( mSitesMovedMoveDown.back() + LEFT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( mSitesMovedMoveDown.back() + LEFT_TOP_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( mSitesMovedMoveDown.back() + LEFT_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( mSitesMovedMoveDown.back() + LEFT_NEIGH_OFFSET );
    //additional items
    mSitesMovedMoveDown.push_back( RIGHT_BOTTOM_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( RIGHT_BOTTOM_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET );
    mSitesMovedMoveDown.push_back( LEFT_BOTTOM_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET );
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
    vectorDist dist( mpExchanger->calcDist( site, destination ) );
    // choose direction
    if ( abs( dist.col ) > abs( dist.row ) )
    {
        if ( dist.col < 0 ) 
        {
            if ( isSpaceToMove( site, mSitesMovedMoveRight ) )
                moveProteinRight( site );
            else return;
        }
        else
        {
            if ( isSpaceToMove( site, mSitesMovedMoveLeft ) )
                moveProteinLeft( site );
            else return;
        }
    }
    else
    {
        if ( dist.row > 0 ) 
        {
            if ( isSpaceToMove( site, mSitesMovedMoveUp ) )
                moveProteinUp( site );
            else return;
        }
        else 
        {
            if ( isSpaceToMove( site, mSitesMovedMoveDown ) )
                moveProteinDown( site );
            else return;
        }
    }
}
void ProteinTriangularLattice::moveProteinRight( lattIndex site )
{
    static const int MOVED_SITES_SIZE = 3;

    lattMember movedSites[ MOVED_SITES_SIZE ];
    lattMember protein( mpLattice[ site ] );
    lattIndex rightSite = site + 1; // rightest site of the protein
    
    movedSites[0] = get( rightSite + RIGHT_NEIGH_OFFSET );
    movedSites[1] = get( rightSite + RIGHT_TOP_NEIGH_OFFSET );
    movedSites[2] = get( rightSite + RIGHT_BOTTOM_NEIGH_OFFSET );

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

    static const int MOVED_SITES_SIZE = 3;
    lattMember movedSites[ MOVED_SITES_SIZE ];
    lattMember protein( mpLattice[ site ] );
    lattIndex leftSite = site - 1; // rightest site of the protein
    
    movedSites[0] = get( leftSite + LEFT_NEIGH_OFFSET );
    movedSites[1] = get( leftSite + LEFT_TOP_NEIGH_OFFSET );
    movedSites[2] = get( leftSite + LEFT_BOTTOM_NEIGH_OFFSET );

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
    static const int MOVED_SITES_SIZE = 3;
    lattMember movedSites[ MOVED_SITES_SIZE ];
    lattMember protein( mpLattice[ site ] );
    lattIndex topSite = site + LEFT_TOP_NEIGH_OFFSET; // topest site of the protein
    
    movedSites[0] = get( topSite + LEFT_TOP_NEIGH_OFFSET );
    movedSites[1] = get( topSite + LEFT_NEIGH_OFFSET );
    movedSites[2] = get( topSite + RIGHT_TOP_NEIGH_OFFSET );

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
    static const int MOVED_SITES_SIZE = 3;
    lattMember movedSites[ MOVED_SITES_SIZE ];
    lattMember protein( mpLattice[ site ] );
    lattIndex bottomSite = site + RIGHT_BOTTOM_NEIGH_OFFSET; // bottomest site of the protein
    
    movedSites[0] = get( bottomSite + RIGHT_BOTTOM_NEIGH_OFFSET );
    movedSites[1] = get( bottomSite + RIGHT_NEIGH_OFFSET );
    movedSites[2] = get( bottomSite + LEFT_BOTTOM_NEIGH_OFFSET );

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


bool ProteinTriangularLattice::isSpaceToMove( lattIndex site, const std::vector< int >& sitesMoved )
{
    std::vector< int >::const_iterator movedSite( sitesMoved.begin() )
                                     , endOfSites( sitesMoved.end() );
    for ( ; movedSite != endOfSites ; ++movedSite )
    {
        if ( get( site + ( *movedSite ) ) == PROTEIN_A 
            || get( site + ( *movedSite ) ) == PROTEIN_B )
        {
            return false;
        }
    }
    return true;

}

