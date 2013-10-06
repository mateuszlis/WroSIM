/*
 * TriangularLatticeTest.cpp
 *
 *  Created on: 24-05-2011
 *      Author: lisu
 */

// std
#include <iostream>
#include <assert.h>
#include <tr1/memory>

// stl
#include <set>

// Google test tools
#include "gtest/gtest.h"
#include "gmock/gmock.h"

// project-local
#include "TriangularLattice.h"

// testing
#include "testUtils.h"
using namespace testUtils;

TEST( TriangularLattice, Construct )
{
    TriangularLattice *latt = new TriangularLattice( 100, 10, 20);
    EXPECT_EQ( calcSum( latt ), 20 * LIPID_B );
    delete latt;
}

TEST( TriangularLattice, ConstructException )
{
    try
    {
        TriangularLattice *latt = new TriangularLattice( 100, 10, 101 );
        // will never happen
        printLatt( latt->getLattice(), 10, 10 );
    }
    catch ( InputParametersException& e )
    {
        SUCCEED();
        return;
    }
    FAIL();
}

TEST( TriangularLatticeTest, GetNeighborIndex )
{
    TriangularLattice *latt = new TriangularLattice( 9, 3, 3 );
    for( int currentSite = 0 ; currentSite < 9 ; ++currentSite )
    {
        std::set< int > neighbors = getNeighborsOf( currentSite, latt );
        EXPECT_EQ( neighbors.size(), 6u ); // check that there are no duplicate entries
        EXPECT_TRUE( checkIfIsNeighborOf( currentSite, neighbors, latt ) ); // check if neighbor relation is reflexive
    }
    delete latt;
}



TEST( TriangularLatticeTest, CalcNeighbors )
{
    TriangularLattice *latt = new TriangularLattice( 100, 10, 1);
    long sum = 0;
    for ( int i = 0 ; i < 100 ; ++i )
    {
        sum += latt->simNeighbCount( i );
    }
    EXPECT_EQ( sum, 6 * 100 - 12 );
    EXPECT_EQ( calcSum( latt ), LIPID_B );
    delete latt;
}

TEST( TriangularLatticeTest, ExchangeSites )
{
    TriangularLattice *latt = new TriangularLattice( 100, 10, 20);

    for ( int i = 0; i < 50000; i++ )
    {
        int pos1 = rand() % latt->getLatticeSize();
        int pos2 = rand() % latt->getLatticeSize();
        if ( pos1 != pos2 )
            latt->exchangeSites( pos1, pos2);
    }
    EXPECT_EQ( calcSum(latt), 20 * LIPID_B );
    delete latt;
}

TEST( TriangularLatticeTest, ClusterAnalysis )
{
    const int LATT_SIZE = 35;
    std::tr1::shared_ptr< TriangularLattice > latt( new TriangularLattice( LATT_SIZE, 5, 1) );
    TriangularLattice::clustersMap map;
    latt->calculateClusters( map );

    EXPECT_EQ( map[ 1 ], 1 );
    for( int i = 0 ; i < LATT_SIZE ; ++i )
    {
        if( latt->getLattice()[i] == LIPID_B )
        {
            latt->getLattice()[ ( i + 1 ) % LATT_SIZE ] = LIPID_B;
            break;
        }
    }
    map.clear();
    latt->calculateClusters( map );
    EXPECT_EQ( map[ 2 ], 1 );
    for( int i = 0 ; i < LATT_SIZE ; ++i )
    {
        if( latt->getLattice()[i] == LIPID_B )
        {
            latt->getLattice()[ ( i + 14 ) % LATT_SIZE ] = LIPID_B;
            break;
        }
    }
    map.clear();
    latt->calculateClusters( map );
    printLatt( latt->getLattice(), 5, 7 );
    EXPECT_EQ( map[ 2 ], 1 );
    EXPECT_EQ( map[ 1 ], 1 );
    latt.reset( new TriangularLattice( LATT_SIZE, 5, 34 ) );
    map.clear();
    latt->calculateClusters( map );
    printLatt( latt->getLattice(), 5, 7 );
    EXPECT_EQ( map[ 34 ], 1 );
}

