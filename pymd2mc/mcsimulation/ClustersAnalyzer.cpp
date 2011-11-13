#include "ClustersAnalyzer.h"

void ClustersAnalyzer::registerAtom( int atomInd, vector< Distance > sortedDistances, int currFrameNum )
{
    int i( 0 );
    while ( sortedDistances[i].d < distThr )
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
            }
            neighPair.second = currFrameNum;
        }
        else
        {
            neighborPairs[ Pair( atomInd, sortedDistances[i].at2Ind ) ] = Pair( 1, currFrameNum );
        }
        ++i;
    }
}

