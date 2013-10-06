#pragma once
// stl
#include <map>
#include <vector>
#include <memory>

// project local
#include "types.h"

using namespace std;


class TriangularLattice;

/**
 * @brief
 *
 **/
class InteractionsTracker
{
    public: // fields

    public: // functions

        /**
         * @brief Ctor - creates structure for protein interaction tracking
         *
         **/
        InteractionsTracker( TriangularLattice* latt, lattIndex* proteins, size_t proteincCount );

        /**
         * @brief Function needs to be called at every move of protein
         *
         **/
        virtual void registerProteinMove( lattIndex site, lattIndex newSite );

        /**
         * @brief returns protein interaction on given site
         **/
        double getInteraction( lattIndex site );


    protected: // typedefs
        typedef std::pair< lattIndex, double > Interaction; /// site->interaction
        typedef std::vector< Interaction > IntVec;
    protected: // functions

    protected: // fields
        TriangularLattice* mpLatt;
        const lattIndex mCols;
        const lattIndex mRows;
        std::unique_ptr< double[] > mpProtInt; // prot interaction tracking

        // this below is constant used as a "map" of interactions for single protein
        const std::vector< std::pair< lattIndex, double > > mIntValues;
}; // class InteractionsTracker
