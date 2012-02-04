/*
 * TriangularLattice.h
 *
 *  Created on: 23-05-2011
 *      Author: lisu
 */

#pragma once
#include <iostream>
#include <string>
#include <iomanip>
#include <cstdlib>

using namespace std;

class TriangularLattice;
ostream &operator<<( ostream &stream, TriangularLattice &latt );

class TriangularLattice
{
    private:
        int *mpLattice;
        int mLatticeSize;
        int mRowSize;
        static const int mNeighbCnt = 6;
        int mNeighb[mNeighbCnt];
        int mNeighbLeft[mNeighbCnt];
        int mNeighbRight[mNeighbCnt];

        void clearArr();
        void distributeParticles( int firstTypeParticlesCnt );

    public:
        TriangularLattice( string filename );
        TriangularLattice( int latticeSize, int rowSize, int firstTypeParticlesCnt );

        int operator[]( int index ) const;

        int getLatticeSize() const;
        int getRowSize() const;

        int simNeighbCount( int pos );
        void exchangeSites( int pos1, int pos2 );

        int getNeighbIndex( int pos, int neighborNum ) const;
        int getNeighborsCnt() const;

        virtual ~TriangularLattice();

        friend ostream &operator<<( ostream &stream, TriangularLattice &latt );
};

