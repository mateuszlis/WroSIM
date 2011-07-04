/*
 * TriangularLatticeTest.cpp
 *
 *  Created on: 24-05-2011
 *      Author: lisu
 */

#include <iostream>
using namespace std;

#include <assert.h>

#include "TriangularLattice.h"

void testTriangularLattice();
void testExchange();
int calcSum( TriangularLattice* );
int main()
{
    srand( ( unsigned )time( 0));
    testTriangularLattice();
    testExchange();

    return 0;
}

void testTriangularLattice()
{
    TriangularLattice *latt = new TriangularLattice( 100, 10, 20);
    assert(calcSum(latt) == 20);
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
int calcSum( TriangularLattice *latt )
{
    int sum = 0;
    for ( int i = 0; i < latt->getLatticeSize(); i++ )
    {
        cout << ( *latt )[i];
        sum += ( *latt )[i];
    }
    return sum;
}
