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
#include "InputParametersException.h"
#include "types.h"

using namespace std;

class TriangularLattice;
ostream &operator<<( ostream &stream, TriangularLattice &latt );

class TriangularLattice
{
    public:
        typedef ::lattMember lattMember; // TODO: not implemented everywhere
        typedef int lattIndex;
        typedef std::map< lattIndex, lattIndex >  clustersMap; // clusterSize ---> cluster count
        typedef map< lattIndex, bool > doneMap;

        TriangularLattice( string filename );

        TriangularLattice( lattIndex latticeSize
                         , lattIndex rowSize
                         , lattIndex firstTypeParticlesCnt
                         , bool distributeRandomly = true );

        lattMember operator[]( lattIndex index ) const;

        lattIndex getLatticeSize() const;
        lattIndex getRowSize() const;
        lattMember* getLattice() const { return mpLattice; };

        lattIndex simNeighbCount( lattIndex pos );
        void exchangeSites( lattIndex pos1, lattIndex pos2 );

        lattIndex getNeighbIndex( lattIndex pos, int neighborNum ) const;
        unsigned int getNeighborsCnt() const;

        void calculateClusters( clustersMap& map );
        virtual ~TriangularLattice();

        friend ostream &operator<<( ostream &stream, TriangularLattice &latt );
    private:
        lattMember *mpLattice;
        lattIndex mLatticeSize;
        lattIndex mRowSize;
        static const int mNeighbCnt = 6;
        lattIndex mNeighb[mNeighbCnt];

        void clearArr();
        void distributeParticlesRandomly( lattIndex firstTypeParticlesCnt );
        void distributeParticles( lattIndex firstTypeParticlesCnt );
        void pushNeighborsToQueue( std::list< lattIndex > & queue, lattIndex siteInd );
		static bool gotDifferentNeighbors(list<int> neighLabels, int currentLabel);
		static int findAncestor(int currentLabel, TriangularLattice::clustersMap& map);

};

