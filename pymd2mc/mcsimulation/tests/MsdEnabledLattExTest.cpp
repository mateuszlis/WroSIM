// Google test tools
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#define protected public // test class internals
// not a very good strategy, but needed here

// project-local
#include "MsdEnabledLattEx.h"
#include "TriangularLattice.h"

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
    EXPECT_DOUBLE_EQ( ex.calcStat(), 2/1000. );
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
    ex.exchangeSites( 4, 5 );
    ex.exchangeSites( 5, 6 );
    ex.exchangeSites( 6, 7 );
    EXPECT_DOUBLE_EQ( ex.calcStat(), 12/25. );
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
    EXPECT_DOUBLE_EQ( ex.calcStat(), 0. );
    delete latt;
}
