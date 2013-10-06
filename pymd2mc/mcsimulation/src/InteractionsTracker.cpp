#include "InteractionsTracker.h"

#include "TriangularLattice.h"
#include <cmath>

using std::sqrt;

InteractionsTracker::InteractionsTracker(
        TriangularLattice* pLatt
        , lattIndex* proteins
        , size_t proteinsCount )
: mpLatt( pLatt )
    , mCols( pLatt->getRowSize() )
    , mRows( pLatt->getLatticeSize() / mCols )
    , mpProtInt( new double[ mRows * mCols ] )
    , mIntValues( // we use this to "map" interactions from single 
        {{ 2, 1 } // protein
        , { -2, 1 }
        , { 3, 0.5 }
        , { -3, 0.5 }
        , { 4, 0.333333}
        , { -4, 0.333333}
        , { mCols + 1, 1 }
        , { mCols + 2, 0.5 }
        , { mCols + 3, 0.333333 }
        , { mCols - 2, 1 }
        , { mCols - 3, 0.5 }
        , { mCols - 4, 0.333333 }
        , { -mCols - 1, 1 }
        , { -mCols - 2, 0.5 }
        , { -mCols - 3, 0.333333 }
        , { -mCols + 2, 1 }
        , { -mCols + 3, 0.5 }
        , { -mCols + 4, 0.333333 }
        , { 2 * mCols, 1 }
        , { 2 * mCols + 1, 0.5 }
        , { 2 * mCols + 2, 0.333333}
        , { 2 * mCols - 1, 1 }
        , { 2 * mCols - 2, 1 }
        , { 2 * mCols - 3, 0.5}
        , { 2 * mCols - 4, 0.333333}
        , { 3 * mCols, 0.5}
        , { 3 * mCols - 1, 0.5}
        , { 3 * mCols - 2, 0.333333}
        , { 3 * mCols + 1, 0.333333}
        , { 4 * mCols, 0.333333 }
        , { 4 * mCols - 1, 0.333333 }
        , { -2 * mCols, 1 }
        , { -2 * mCols - 1, 0.5 }
        , { -2 * mCols - 2, 0.333333}
        , { -2 * mCols + 1, 1 }
        , { -2 * mCols + 2, 1 }
        , { -2 * mCols + 3, 0.5}
        , { -2 * mCols + 4, 0.333333}
        , { -3 * mCols, 0.5 }
        , { -3 * mCols + 1, 0.5 }
        , { -3 * mCols - 1, 0.333333 }
        , { -3 * mCols + 2, 0.333333 }
        , { -4 * mCols, 0.333333 }
        , { -4 * mCols + 1, 0.333333 }
        })
{
    for ( lattIndex i( 0 ) ; i < mRows * mCols ; ++i )
    {
        mpProtInt[i] = 0;
    }

    for ( size_t prot( 0 ) ; prot < proteinsCount ; ++prot )
    {
        for ( auto intVal : mIntValues )
        {
            mpProtInt[ mpLatt->cutToPos( intVal.first + proteins[ prot ] ) ] += intVal.second; 
        }
    }

}

void InteractionsTracker::registerProteinMove( lattIndex site
                                            , lattIndex newSite )
{
    for ( auto intVal : mIntValues )
    {
        mpProtInt[ mpLatt->cutToPos( intVal.first + site ) ] -= intVal.second; 
        mpProtInt[ mpLatt->cutToPos( intVal.first + newSite ) ] += intVal.second; 
    }
}

double InteractionsTracker::getInteraction( lattIndex site )
{
    return mpProtInt[ site ];
}
