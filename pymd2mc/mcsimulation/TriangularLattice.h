/*
 * TriangularLattice.h
 *
 *  Created on: 23-05-2011
 *      Author: lisu
 */

#pragma once
// std
#include <iostream>
#include <string>
#include <iomanip>
#include <cstdlib>


// stl
#include <map>
#include <list>

// project local
#include "types.h"

using namespace std;

class TriangularLattice;
ostream &operator<<( ostream &stream, TriangularLattice &latt );

class TriangularLattice
{
    public:
        typedef ::lattMember lattMember; // TODO: not implemented everywhere
        typedef std::map< unsigned int, unsigned int >  clustersMap;
        typedef map< unsigned int, unsigned int > doneMap;

        TriangularLattice( string filename );
        TriangularLattice( int latticeSize, int rowSize, int firstTypeParticlesCnt );

        lattMember operator[]( int index ) const;

        int getLatticeSize() const;
        int getRowSize() const;
        lattMember* getLattice() const { return mpLattice; };

        int simNeighbCount( int pos );
        void exchangeSites( int pos1, int pos2 );

        int getNeighbIndex( int pos, int neighborNum ) const;
        int getNeighborsCnt() const;

        void calculateClusters( clustersMap& map );
        virtual ~TriangularLattice();

        friend ostream &operator<<( ostream &stream, TriangularLattice &latt );
    private:
        lattMember *mpLattice;
        int mLatticeSize;
        int mRowSize;
        static const int mNeighbCnt = 6;
        int mNeighb[mNeighbCnt];
        int mNeighbLeft[mNeighbCnt];
        int mNeighbRight[mNeighbCnt];

        void clearArr();
        void distributeParticles( int firstTypeParticlesCnt );
        void pushNeighborsToQueue( std::list< lattMember > & queue, unsigned int siteInd );

};

