// Google test tools
#include "gtest/gtest.h"
#include "gmock/gmock.h"

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
    ex.exchangeSites( 4, 6 );
    EXPECT_DOUBLE_EQ( ex.calcStat(), 0.01 );
    delete latt;
}
