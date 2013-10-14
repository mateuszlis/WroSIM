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

class LattExchanger;
class TriangularLattice;
ostream &operator<<( ostream &stream, TriangularLattice &latt );

class TriangularLattice
{
    public: // typedefs
        typedef ::lattMember lattMember; // TODO: not implemented everywhere
        typedef int lattIndex;
        typedef std::map< lattIndex, lattIndex >  clustersMap; // clusterSize ---> cluster count
        typedef map< lattIndex, bool > doneMap;

    public: // functions
        TriangularLattice( string filename ); //TODO: not yet implemented

        TriangularLattice( lattIndex latticeSize
                         , lattIndex rowSize
                         , lattIndex firstTypeParticlesCnt
                         , bool distributeRandomly = true );

        TriangularLattice( const TriangularLattice& );

        lattMember &operator[]( lattIndex index ) const;
        lattMember get( int index ) const; // applies helical boundary conditions

        lattIndex getLatticeSize() const;
        lattIndex getRowSize() const;
        lattMember* getLattice() const { return mpLattice; };
        lattIndex getNonLipidMembersCount() const { return mNonLipidMembersCount; };

        lattIndex simNeighbCount( lattIndex pos );
        void exchangeSites( lattIndex pos1, lattIndex pos2 );

        lattIndex getNeighbIndex( lattIndex pos, int neighborNum ) const;
        int getNeighborsCnt() const;

        void calculateClusters( clustersMap& map );
        void setExchanger( LattExchanger* );
        LattExchanger* getExchanger();
        virtual lattIndex cutToPos( int pos ) const;
        static const int mNeighbCnt = 6;
        virtual ~TriangularLattice();

        friend ostream &operator<<( ostream &stream, TriangularLattice &latt );

    protected:
        lattMember *mpLattice;
        lattIndex mLatticeSize;
        lattIndex mRowSize;
        LattExchanger* mpExchanger;
        bool selfLattExchanger;
        lattIndex mNeighb[mNeighbCnt];
        lattIndex mNonLipidMembersCount;

        virtual void clearArr();
        virtual void distributeParticlesRandomly( lattIndex firstTypeParticlesCnt );
        virtual void distributeParticles( lattIndex firstTypeParticlesCnt );
        virtual bool isFree( lattIndex pos ) const;
        void pushNeighborsToQueue( std::list< lattIndex > & queue, lattIndex siteInd );
        void clearExchanger();
		static bool gotDifferentNeighbors(list<int> neighLabels, int currentLabel);
		static int findAncestor(int currentLabel, TriangularLattice::clustersMap& map);

};

