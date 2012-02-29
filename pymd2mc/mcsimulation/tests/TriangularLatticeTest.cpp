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

// helpers
void printLatt( lattMember *latt, int rowSize, int rowsCount );
int calcSum( TriangularLattice* );
std::set< int > getNeighborsOf( int site, TriangularLattice* latt );
bool checkIfIsNeighborOf( const int& neighbor, const std::set< int >& elements, TriangularLattice* latt );


TEST( TriangularLattice, Construct )
{
    TriangularLattice *latt = new TriangularLattice( 100, 10, 20);
    EXPECT_EQ( calcSum( latt ), 20 );
    delete latt;
}

TEST( TriangularLatticeTest, GetNeighborIndex )
{
    TriangularLattice *latt = new TriangularLattice( 9, 3, 3 );
    for( int currentSite = 0 ; currentSite < 9 ; ++currentSite )
    {
        std::set< int > neighbors = getNeighborsOf( currentSite, latt );
        EXPECT_EQ( neighbors.size(), 6 ); // check that there are no duplicate entries
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
    EXPECT_EQ( calcSum( latt ), 1);
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
    EXPECT_EQ( calcSum(latt), 20 );
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
        if( latt->getLattice()[i] == 1 )
        {
            latt->getLattice()[ ( i + 1 ) % LATT_SIZE ] = 1;
            break;
        }
    }
    map.clear();
    latt->calculateClusters( map );
    EXPECT_EQ( map[ 2 ], 1 );
    for( int i = 0 ; i < LATT_SIZE ; ++i )
    {
        if( latt->getLattice()[i] == 1 )
        {
            latt->getLattice()[ ( i + 14 ) % LATT_SIZE ] = 1;
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
// helpers

void printLatt( lattMember *latt, int rowSize, int rowsCount )
{
    for ( int i = 0 ; i < rowsCount ; ++i )
    {
        cout << setw( 2 * ( i + 1 ) ) << "  " ;
        for( int j = 0 ; j < rowSize; ++j )
        {   
            cout << setw(2) << static_cast< int >( latt[ i * rowSize + j ] ) << "  ";
        }
        cout << endl;
    }
}
int calcSum( TriangularLattice *latt )
{
    int sum = 0;
    for ( int i = 0; i < latt->getLatticeSize(); i++ )
    {
        sum += ( *latt )[i];
    }
    return sum;
}

std::set< int > getNeighborsOf( int site, TriangularLattice* latt )
{
    std::set< int > neighbors;
    for( int neighNum = 0 ; neighNum < 6 ; ++neighNum )
    {
        neighbors.insert( latt->getNeighbIndex( site, neighNum ) );
    }
    return neighbors;
}

/**
 * @brief for each element from elements, checks if neighbor is its neighbor on
 * the lattice latt
 */
bool checkIfIsNeighborOf( const int& neighbor, const std::set< int >& elements, TriangularLattice* latt )
{
    for ( std::set< int >::const_iterator it = elements.begin() 
            ; it != elements.end() 
            ; ++it )
    {
        std::set< int > elementsNeighbors = getNeighborsOf( *it, latt );
        if ( elementsNeighbors.find( neighbor ) == elementsNeighbors.end() )
        {
            return false;
        }
    }
    return true;
}

