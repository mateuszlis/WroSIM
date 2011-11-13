/*
 * ClustersAnalysis.h
 *
 *  Created on: 13-11-2011
 *      Author: lisu
 */

#pragma once

#include <iostream>
#include <string>
#include <cstdlib>
#include <map>
#include <vector>
#include "Atom.h"

using namespace std;

class ClustersAnalyzer
{
    public: 
        ClustersAnalyzer( int frameNumThr, double distThr ):
            frameNumThr(frameNumThr), distThr(distThr)
        {}

        void registerAtom( int atomInd, vector< Distance > sortedDistances, int currFrameNum );
        bool isClustered( int atomInd, int currFrameNum )
        {}

        bool isInMixedCluster( int atomInd, int currFrameNum )
        {}
        friend ostream& operator<<( ostream& s, ClustersAnalyzer& an )
        {
            NeighbMem::iterator it = an.neighborPairs.begin();
            for ( ; it != an.neighborPairs.end() ; ++it )
            {
                s << "[" << it->first.first << "," << it->first.second << "] : " << it->second.first << ", " << it->second.second << endl;
            }
            return s;
        }

        typedef pair< int, int > Pair;
        typedef map< Pair, Pair >  NeighbMem;
        NeighbMem neighborPairs; // std::pair< atomInd1, atomInd2 > => std::pair< numOfFrames, lastFrameNum >

    private:
        int frameNumThr;
        double distThr;
};



