#include "TriangularLattice.h"
#include "LattExchanger.h"

LattExchanger::LattExchanger( TriangularLattice* latt, lattIndex *proteins, lattIndex proteinsCnt )
            : mpLatt( latt ) 
            , RIGHT_NEIGH_OFFSET( 1 )
            , RIGHT_TOP_NEIGH_OFFSET( -latt->getRowSize() + 1 )
            , RIGHT_BOTTOM_NEIGH_OFFSET( latt->getRowSize() ) 
            , LEFT_BOTTOM_NEIGH_OFFSET( latt->getRowSize() - 1 )
            , LEFT_TOP_NEIGH_OFFSET( -latt->getRowSize() )
            , LEFT_NEIGH_OFFSET( -1 )
            , mProteins( proteins )
            , mProteinCnt( proteinsCnt )
    {
        mProteinSites[ 0 ] = 0;
        mProteinSites[ 1 ] = 1;
        mProteinSites[ 2 ] = -1;
        mProteinSites[ 3 ] = mpLatt->getRowSize();
        mProteinSites[ 4 ] = -mpLatt->getRowSize();
        mProteinSites[ 5 ] = mpLatt->getRowSize() - 1;
        mProteinSites[ 6 ] = -mpLatt->getRowSize() + 1; // TODO: this fixed the protein size

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
void LattExchanger::exchangeSites( lattIndex pos1, lattIndex pos2 )
{
    // fast variable values exchange trick (pos1 != pos2)
	mpLatt->getLattice()[pos1] ^= mpLatt->getLattice()[pos2];
    mpLatt->getLattice()[pos2] ^= mpLatt->getLattice()[pos1];
    mpLatt->getLattice()[pos1] ^= mpLatt->getLattice()[pos2];
}

void LattExchanger::moveProtein( lattIndex site, lattIndex destination )
{
    vectorDist dist( calcDist( site, destination ) );
    // choose direction
    if ( abs( dist.col ) > abs( dist.row ) )
    {
        if ( 
                ( dist.col < 0 && dist.col > -mpLatt->getRowSize() + 2 )
                || dist.col > mpLatt->getRowSize() - 4) // PBC jump
        {
            if ( isSpaceToMove( site, mSitesMovedMoveRight ) )
            {
                moveProteinRight( site );
            }
            else return;
        }
        else
        {
            if ( isSpaceToMove( site, mSitesMovedMoveLeft ) )
            {
                moveProteinLeft( site );
            }
            else return;
        }
    } 
    else
    {
    lattIndex latticeColSize( mpLatt->getLatticeSize() / mpLatt->getRowSize() );
        if ( 
              ( dist.row < 0 && dist.row > -latticeColSize + 2 )
              || dist.row > latticeColSize - 4 )  // PBC jump
        {
            if ( isSpaceToMove( site, mSitesMovedMoveDown ) )
            {
                moveProteinDown( site );
            }
            else return;
        }
        else 
        {
            if ( isSpaceToMove( site, mSitesMovedMoveUp ) )
            {
                moveProteinUp( site );
            }
            else return;
        }
    }
}

vectorDist LattExchanger::calcDist( lattIndex pos1, lattIndex pos2 )
{
    lattIndex startRow = pos1 / mpLatt->getRowSize();
    lattIndex startCol = pos1 % mpLatt->getRowSize();
    lattIndex endRow = pos2 / mpLatt->getRowSize();
    lattIndex endCol = pos2 % mpLatt->getRowSize();
    return vectorDist( startRow - endRow, startCol - endCol );

}

void LattExchanger::moveProteinRight( lattIndex site )
{
    static const int MOVED_SITES_SIZE = 3;

    int movedSites[ MOVED_SITES_SIZE ];
    lattIndex rightSite = site + 1; // rightest site of the protein
    
    movedSites[0] = mpLatt->get( rightSite + RIGHT_NEIGH_OFFSET );
    movedSites[1] = mpLatt->get( rightSite + RIGHT_TOP_NEIGH_OFFSET );
    movedSites[2] = mpLatt->get( rightSite + RIGHT_BOTTOM_NEIGH_OFFSET );

    // move right top site
    initPushAndPop( rightSite + RIGHT_TOP_NEIGH_OFFSET );
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + LEFT_TOP_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET, movedSites[1] );
    push( site + RIGHT_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET, movedSites[1] );


    // move right bottom site
    initPushAndPop( rightSite + RIGHT_BOTTOM_NEIGH_OFFSET );
    pushAndPop( site + RIGHT_BOTTOM_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + RIGHT_BOTTOM_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + RIGHT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    push( site + RIGHT_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET, movedSites[2] );

    // move center

    initPushAndPop( rightSite + RIGHT_NEIGH_OFFSET );
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
    lastMovedSite += RIGHT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_NEIGH_OFFSET;
    push( lastMovedSite, movedSites[0] );

    // move the protein
    lattIndex newProteinCenter( mpLatt->cutToPos( site + RIGHT_NEIGH_OFFSET ) );
    //( *mpLatt )[ newProteinCenter ] = protein;
    //for ( int i = 0 ; i < TriangularLattice::mNeighbCnt ; ++i )
    //{
    //    ( *mpLatt )[ mpLatt->getNeighbIndex( newProteinCenter, i ) ] = protein;
    //}
    // update array
    //FIXME: this might be efficiency drawback
    updateProteinArray( site, newProteinCenter );

}


void LattExchanger::moveProteinLeft( lattIndex site )
{

    static const int MOVED_SITES_SIZE = 3;
    int movedSites[ MOVED_SITES_SIZE ];
    lattIndex leftSite = site - 1; // rightest site of the protein
    
    movedSites[0] = mpLatt->get( leftSite + LEFT_NEIGH_OFFSET );
    movedSites[1] = mpLatt->get( leftSite + LEFT_TOP_NEIGH_OFFSET );
    movedSites[2] = mpLatt->get( leftSite + LEFT_BOTTOM_NEIGH_OFFSET );

    // move left top site
    initPushAndPop( leftSite + LEFT_TOP_NEIGH_OFFSET );
    pushAndPop( site + LEFT_TOP_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + LEFT_TOP_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + LEFT_TOP_NEIGH_OFFSET, movedSites[1] );
    push( site + LEFT_TOP_NEIGH_OFFSET + LEFT_NEIGH_OFFSET, movedSites[1] );


    // move left bottom site
    initPushAndPop( leftSite + LEFT_BOTTOM_NEIGH_OFFSET );
    pushAndPop( site + LEFT_BOTTOM_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + LEFT_BOTTOM_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + RIGHT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + LEFT_BOTTOM_NEIGH_OFFSET + LEFT_NEIGH_OFFSET, movedSites[2] );

    // move center
    initPushAndPop( leftSite + LEFT_NEIGH_OFFSET );
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
    lastMovedSite += LEFT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_NEIGH_OFFSET;
    push( lastMovedSite, movedSites[0] );

    // move the protein
    lattIndex newProteinCenter( mpLatt->cutToPos( site + LEFT_NEIGH_OFFSET ) );
    //( *mpLatt )[ newProteinCenter ] = protein;
    //for ( int i = 0 ; i < TriangularLattice::mNeighbCnt ; ++i )
    //{
    //    ( *mpLatt )[ mpLatt->getNeighbIndex( newProteinCenter, i ) ] = protein;
    //}
    //// update array
    updateProteinArray( site, newProteinCenter );

}

void LattExchanger::moveProteinUp( lattIndex site )
{
    static const int MOVED_SITES_SIZE = 3;
    int movedSites[ MOVED_SITES_SIZE ];
    lattIndex topSite = site + LEFT_TOP_NEIGH_OFFSET; // topest site of the protein
    
    movedSites[0] = mpLatt->get( topSite + LEFT_TOP_NEIGH_OFFSET );
    movedSites[1] = mpLatt->get( topSite + LEFT_NEIGH_OFFSET );
    movedSites[2] = mpLatt->get( topSite + RIGHT_TOP_NEIGH_OFFSET );

    // move left top site 
    initPushAndPop( topSite + LEFT_NEIGH_OFFSET );
    pushAndPop( site + LEFT_NEIGH_OFFSET + LEFT_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + LEFT_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + LEFT_NEIGH_OFFSET, movedSites[1] );
    push( site + LEFT_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET, movedSites[1] );


    // move left bottom site
    initPushAndPop( topSite + RIGHT_TOP_NEIGH_OFFSET );
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + RIGHT_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET, movedSites[2] );
    push( site + RIGHT_TOP_NEIGH_OFFSET + LEFT_TOP_NEIGH_OFFSET, movedSites[2] );

    // move center
    initPushAndPop( topSite + LEFT_TOP_NEIGH_OFFSET );
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
    lastMovedSite += LEFT_TOP_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_TOP_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += LEFT_TOP_NEIGH_OFFSET;
    push( lastMovedSite, movedSites[0] );

    // move the protein
    lattIndex newProteinCenter( mpLatt->cutToPos( site + LEFT_TOP_NEIGH_OFFSET ) );
    //( *mpLatt )[ newProteinCenter ] = protein;
    //for ( int i = 0 ; i < TriangularLattice::mNeighbCnt ; ++i )
    //{
    //    ( *mpLatt )[ mpLatt->getNeighbIndex( newProteinCenter, i ) ] = protein;
    //}
    // update array
    updateProteinArray( site, newProteinCenter );

}

void LattExchanger::moveProteinDown( lattIndex site )
{
    static const int MOVED_SITES_SIZE = 3;
    int movedSites[ MOVED_SITES_SIZE ];
    lattIndex bottomSite = site + RIGHT_BOTTOM_NEIGH_OFFSET; // bottomest site of the protein
    
    movedSites[0] = mpLatt->get( bottomSite + RIGHT_BOTTOM_NEIGH_OFFSET );
    movedSites[1] = mpLatt->get( bottomSite + RIGHT_NEIGH_OFFSET );
    movedSites[2] = mpLatt->get( bottomSite + LEFT_BOTTOM_NEIGH_OFFSET );

    // move left bottom site 
    initPushAndPop( bottomSite + RIGHT_NEIGH_OFFSET );
    pushAndPop( site + RIGHT_NEIGH_OFFSET + RIGHT_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + RIGHT_NEIGH_OFFSET + RIGHT_TOP_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + RIGHT_TOP_NEIGH_OFFSET, movedSites[1] );
    pushAndPop( site + RIGHT_NEIGH_OFFSET, movedSites[1] );
    push( site + RIGHT_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET, movedSites[1] );


    // move left bottom site
    initPushAndPop( bottomSite + LEFT_BOTTOM_NEIGH_OFFSET );
    pushAndPop( site + LEFT_BOTTOM_NEIGH_OFFSET + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + LEFT_TOP_NEIGH_OFFSET + LEFT_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + LEFT_NEIGH_OFFSET, movedSites[2] );
    pushAndPop( site + LEFT_BOTTOM_NEIGH_OFFSET, movedSites[1] );
    push( site + LEFT_BOTTOM_NEIGH_OFFSET + RIGHT_BOTTOM_NEIGH_OFFSET, movedSites[1] );

    // move center
    initPushAndPop( bottomSite + RIGHT_BOTTOM_NEIGH_OFFSET );
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
    lastMovedSite += RIGHT_BOTTOM_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_BOTTOM_NEIGH_OFFSET;
    pushAndPop( lastMovedSite, movedSites[0] );
    lastMovedSite += RIGHT_BOTTOM_NEIGH_OFFSET;
    push( lastMovedSite, movedSites[0] );

    // move the protein
    lattIndex newProteinCenter( mpLatt->cutToPos( site + RIGHT_BOTTOM_NEIGH_OFFSET ) );
    //( *mpLatt )[ newProteinCenter ] = protein;
    //for ( int i = 0 ; i < TriangularLattice::mNeighbCnt ; ++i )
    //{
    //    ( *mpLatt )[ mpLatt->getNeighbIndex( newProteinCenter, i ) ] = protein;
    //}
    // update array
    updateProteinArray( site, newProteinCenter );

}

void LattExchanger::updateProteinArray( lattIndex site, lattIndex newPos )
{
    for ( int i = 0 ; i < mProteinCnt ; ++i )
    {
        if ( mProteins[i] == site )
        {
            mProteins[i] = newPos;
        }
    }
}

void LattExchanger::pushAndPop( int site, int &value )
{
    site = mpLatt->cutToPos( site );
    lattMember temp = value;
    value = mpLatt->getLattice()[ site ];
    mpLatt->getLattice()[ site ] = temp;
}

void LattExchanger::initPushAndPop( int /*site*/ )
{
}

void LattExchanger::push( int site, int & value )
{
    site = mpLatt->cutToPos( site );
    mpLatt->getLattice()[ site ] = value;
}
bool LattExchanger::isSpaceToMove( lattIndex site, const std::vector< int >& sitesMoved ) const
{
    std::vector< int >::const_iterator movedSite( sitesMoved.begin() )
                                     , endOfSites( sitesMoved.end() );
    for ( ; movedSite != endOfSites ; ++movedSite )
    {
        if ( mpLatt->get( site + ( *movedSite ) ) == PROTEIN_A 
            || mpLatt->get( site + ( *movedSite ) ) == PROTEIN_B )
        {
            return false;
        }
    }
    return true;

}
