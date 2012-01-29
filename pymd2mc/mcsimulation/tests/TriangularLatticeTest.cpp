/*
 * TriangularLatticeTest.cpp
 *
 *  Created on: 24-05-2011
 *      Author: lisu
 */

// std
#include <iostream>
#include <assert.h>

// stl
#include <set>

// project-local
#include "TriangularLattice.h"

void testTriangularLattice();
void testExchange();
void testCalcNeigh();
void testGetNeighbIndex();

// helpers
int calcSum( TriangularLattice* );
std::set< int > getNeighborsOf( int site, TriangularLattice* latt );
bool checkIfIsNeighborOf( const int& neighbor, const std::set< int >& elements, TriangularLattice* latt );

int main()
{
    testTriangularLattice();
    testGetNeighbIndex();
    testExchange();
    testCalcNeigh();

    return 0;
}

void testTriangularLattice()
{
    TriangularLattice *latt = new TriangularLattice( 100, 10, 20);
    assert(calcSum(latt) == 20);
    delete latt;
}

void testGetNeighbIndex()
{
    TriangularLattice *latt = new TriangularLattice( 9, 3, 3 );
    for( int currentSite = 0 ; currentSite < 9 ; ++currentSite )
    {
        std::set< int > neighbors = getNeighborsOf( currentSite, latt );
        assert( neighbors.size() == 6 ); // check that there are no duplicate entries
        assert( checkIfIsNeighborOf( currentSite, neighbors, latt ) ); // check if neighbor relation is reflexive
    }    
    delete latt;
}



void testCalcNeigh()
{
    TriangularLattice *latt = new TriangularLattice( 100, 10, 1);
    long sum = 0;
    for ( int i = 0 ; i < 100 ; ++i )
    {
        sum += latt->simNeighbCount( i );
    }
    assert( sum == 6 * 100 - 12 );
    assert(calcSum(latt) == 1);
    delete latt;
}

void testExchange()
{
    TriangularLattice *latt = new TriangularLattice( 100, 10, 20);

    for ( int i = 0; i < 50000; i++ )
    {
        int pos1 = rand() % latt->getLatticeSize();
        int pos2 = rand() % latt->getLatticeSize();
        if ( pos1 != pos2 )
            latt->exchangeSites( pos1, pos2);
    }
    assert(calcSum(latt) == 20);
    delete latt;
}

// helpers

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

