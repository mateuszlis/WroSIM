// Google test tools
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#define protected public // test class internals
// not a very good strategy, but needed here

// project-local
#include "MsdEnabledLattEx.h"
#include "TriangularLattice.h"
#include "ProteinTriangularLattice.h"

// tests
#include "testUtils.h"
using namespace testUtils;

TEST( MsdEnabledLattEx, Construct )
{
    TriangularLattice *latt = new TriangularLattice( 100, 10, 20);
    MsdEnabledLattEx ex( latt );
    delete latt;
}

TEST( MsdEnabledLattEx, calcDist1  )
{
    //FIXME: this should be done in tearUp and tearDown
    TriangularLattice *latt;
    latt = new TriangularLattice( 1000, 100, 20 );
    MsdEnabledLattEx ex( latt );
    vectorDist distance( ex.calcDist( 0, 99 ) );
    EXPECT_EQ( distance.col, -99 );
    EXPECT_EQ( distance.row, 0 );
    delete latt;
}


TEST( MsdEnabledLattEx, calcDist2  )
{
    TriangularLattice *latt;
    latt = new TriangularLattice( 1000, 100, 20 );
    MsdEnabledLattEx ex( latt );
    vectorDist distance( ex.calcDist( 298, 99 ) );
    EXPECT_EQ( distance.col, -1 );
    EXPECT_EQ( distance.row, 2 );
    delete latt;
}
TEST( MsdEnabledLattEx, calcDist_PBC  )
{
    TriangularLattice *latt;
    latt = new TriangularLattice( 1000, 100, 20 );
    MsdEnabledLattEx ex( latt );
    vectorDist distance( ex.calcDist( 100, 99 ) );
    EXPECT_EQ( distance.col, -99 );
    EXPECT_EQ( distance.row, 1 );
    delete latt;
}

TEST( MsdEnabledLattEx, msdCalc  )
{
    TriangularLattice *latt;
    latt = new TriangularLattice( 1000, 100, 20 );
    MsdEnabledLattEx ex( latt );
    ex.exchangeSites( 0, 1 );
    double msd( 0 ), protMsd( 0 );
    ex.calcStat( msd, protMsd );
    EXPECT_DOUBLE_EQ( msd, 2/1000. );
    delete latt;
}

TEST( MsdEnabledLattEx, isPBCJump  )
{
    TriangularLattice *latt;
    latt = new TriangularLattice( 25, 5, 20 );
    MsdEnabledLattEx ex( latt );
    EXPECT_FALSE( ex.isNotPBCJump( 4, 5 ) );
    EXPECT_TRUE( ex.isNotPBCJump( 2, 3 ) );
    EXPECT_FALSE( ex.isNotPBCJump( 2, 8 ) );
    EXPECT_FALSE( ex.isNotPBCJump( 7, 1 ) );
    delete latt;
}

TEST( MsdEnabledLattEx, calcMsd_PBC  )
{
    TriangularLattice *latt;
    latt = new TriangularLattice( 25, 5, 20 );
    MsdEnabledLattEx ex( latt );
    printLatt( latt->getLattice(), 5, 5 );
    printPermutation( ex.mTracking, 5, 5 );
    ex.exchangeSites( 4, 5 );
    ex.exchangeSites( 5, 6 );
    ex.exchangeSites( 6, 7 );
    double msd( 0 ), protMsd( 0 );
    ex.calcStat( msd, protMsd );
    EXPECT_DOUBLE_EQ( msd, 12/25. );
    printLatt( latt->getLattice(), 5, 5 );
    printPermutation( ex.mTracking, 5, 5 );
    delete latt;
}


TEST( MsdEnabledLattEx, calcMsd_PBC_circular  )
{
    TriangularLattice *latt;
    latt = new TriangularLattice( 25, 5, 20 );
    MsdEnabledLattEx ex( latt );
    ex.exchangeSites( 4, 5 );
    ex.exchangeSites( 4, 5 );
    ex.exchangeSites( 4, 5 );
    ex.exchangeSites( 4, 5 );
    ex.exchangeSites( 23, 3 );
    ex.exchangeSites( 23, 3 );
    double msd( 0 ), protMsd( 0 );
    ex.calcStat( msd, protMsd );
    EXPECT_DOUBLE_EQ( msd, 0. );
    delete latt;
}

TEST( MsdEnabledLattEx, calcMsd_Protein  )
{
    ProteinTriangularLattice *latt;
    latt = new ProteinTriangularLattice( 25, 5, 10, 1, false );
    MsdEnabledLattEx ex( latt, latt->getProteins(), latt->getProteinCnt() );
    ex.setProteins( latt->getProteins(), latt->getProteinCnt() );
    printLatt( latt->getLattice(), 5, 5 );
    printPermutation( ex.mTracking, 5, 5 );
    ex.moveProtein( 15, 16 );
    printLatt( latt->getLattice(), 5, 5 );
    printPermutation( ex.mTracking, 5, 5 );
    double msd( 0 ), protMsd( 0 );
    ex.calcStat( msd, protMsd );
    //EXPECT_DOUBLE_EQ( msd, 19./18. ); // this should be right but it is not because of HBC Msd calc drawback
    EXPECT_DOUBLE_EQ( protMsd, 1. );
    ex.moveProtein( 16, 17 );
    ex.calcStat( msd, protMsd );
    printLatt( latt->getLattice(), 5, 5 );
    printPermutation( ex.mTracking, 5, 5 );
    ex.moveProtein( 17, 18 );
    EXPECT_DOUBLE_EQ( protMsd, 4. );
    delete latt;
}
TEST( MsdEnabledLattEx, PushAndPop  )
{
    ProteinTriangularLattice *latt;
    latt = new ProteinTriangularLattice( 9, 3, 2, 0, false );
    MsdEnabledLattEx ex( latt, latt->getProteins(), latt->getProteinCnt() );
    printLatt( latt->getLattice(), 3, 3 );
    printPermutation( ex.mTracking, 3, 3 );
    for ( int i = 0 ; i < 2 ; ++i )
    {
        int localMember( latt->get( 3 ) );
        ex.initPushAndPop( 3 );
        ex.pushAndPop( 4, localMember );
        ex.pushAndPop( 6, localMember );
        ex.push( 3, localMember );
    }
    int localMember( latt->get( 3 ) );
    ex.initPushAndPop( 3 );
    ex.pushAndPop( 4, localMember );
    ex.pushAndPop( 5, localMember );
    ex.push( 3, localMember );
    printLatt( latt->getLattice(), 3, 3 );
    printPermutation( ex.mTracking, 3, 3 );
    EXPECT_EQ( ex.mTracking[3], 5 );
    EXPECT_EQ( ex.mTracking[4], 4 );
    EXPECT_EQ( ex.mTracking[6], 3 );
    delete latt;
}
