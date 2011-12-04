#include "ClustersAnalyzer.h"

void ClustersAnalyzer::registerAtom( int atomInd, vector< Distance > sortedDistances, int currFrameNum, Atom* atoms )
{
    unsigned int i( 0 );
    string symbol;
    bool isMixed = false;
    if ( atoms ) 
        symbol = atoms[atomInd].resname + atoms[atomInd].name;
    while ( sortedDistances.size() > i && sortedDistances[i].d < distThr )
    {
        NeighbMem::iterator it = neighborPairs.find( Pair( atomInd, sortedDistances[i].at2Ind ) );
        if ( it != neighborPairs.end() )
        {
            Pair &neighPair = it->second;
            if ( neighPair.second < currFrameNum - 1 )
            {
                neighPair.first = 1;
            }
            else
            {
                neighPair.first++;
                if ( neighPair.first > frameNumThr )
                { // at this point we know, that we have close neighbors for long enough
                    // so we use that to mark it as a clustered atom in current
                    // frame and if neighbor is different, set clustertype to
                    // mixed.
                    isClusteredMap[atomInd] = pair< bool, int >( true, currFrameNum );
                    if ( atoms && ( symbol != atoms[it->first.second].resname + atoms[it->first.second].name ) )
                        isMixed = true;
                }
            }
            neighPair.second = currFrameNum;
        }
        else
        {
            neighborPairs[ Pair( atomInd, sortedDistances[i].at2Ind ) ] = Pair( 1, currFrameNum );
        }
        ++i;
    }
    if ( isClustered( atomInd, currFrameNum ) )
    {
        isMixedMap[atomInd].first = isMixed;
        isMixedMap[atomInd].second = currFrameNum;
    }
}

bool ClustersAnalyzer::isClustered( int atomInd, int currFrameNum )
{
    mapAtomNumClust::iterator it = isClusteredMap.find( atomInd );
    if ( it != isClusteredMap.end() )
    {
        if ( ( it->second.second ) == currFrameNum && it->second.first )
            return true;
    }
    return false;
}
bool ClustersAnalyzer::isInMixedCluster( int atomInd, int currFrameNum )
{
    mapAtomNumClust::iterator it = isMixedMap.find( atomInd );
    if ( it != isMixedMap.end() )
    {
        if ( ( it->second.second ) == currFrameNum && it->second.first )
            return true;
    }
    return false;
}
