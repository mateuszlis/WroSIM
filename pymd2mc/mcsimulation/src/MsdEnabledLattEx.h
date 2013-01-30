#include "LattExchanger.h"

#include <math.h>
using namespace std;


/**
 * @brief Class responsible for tracking everything thata happens in the lattice. It is able to calculate Mean Square Displacement based on that information.
 *
 **/
class MsdEnabledLattEx : public LattExchanger
{
    public: // functions
        /**
         * @brief Ctor - assigns tracked lattice (not possible to change). Allocates memory for tracking
         *
         * @param latt - lattice that will be tracked
         **/
        MsdEnabledLattEx( TriangularLattice* latt, lattIndex * proteins = NULL, lattIndex proteinCnt = 0 );

        /**
         * @brief Since currently it is responsibility of Metropolis class to output statistics, this function is used to determine, whether we have Msd stat enable or not. Should be organized in different way.
         *
         **/
        virtual bool hasMsd() const { return true; }

        /**
         * @brief This function is used by Lattice to facilitate tracking of lattice movement
         **/
        virtual void exchangeSites( lattIndex pos1, lattIndex pos2 );

        /**
         * @brief This function is used by Lattice to facilitate tracking of Protein movement
         **/
        virtual void moveProtein( lattIndex site, lattIndex destination );

        /**
         * @brief Calculates Mean Square Displacement
         **/
        virtual void calcStat( double & msd, double & protMsd );

        virtual void setProteins( lattIndex *proteins, lattIndex proteinCnt )
        {
            LattExchanger::setProteins( proteins, proteinCnt );
            for ( lattIndex i(0) ; i < mpLatt->getLatticeSize() ; ++i )
            {
                mWasProtein[i] = isProtein( mpLatt->get(i) );
            }
        }


        virtual ~MsdEnabledLattEx();

    protected: // fields

        static lattIndex mLastPos; ///< for pushAndPop status for PBC
        static lattIndex mLastValue; ///< for pushAndPop status for PBC
        lattIndex* mTracking; ///< track current position of every lipid
        bool* mWasProtein; ///< track initial status
        vectorDist* mPBCCorrection; ///< track PBC jumps
        static const bool isNeighbor[3][3]; ///< used to determine if sites are direct neighbors

    protected: // functions

        bool isNotPBCJump( lattIndex pos1, lattIndex pos2 );
        void incDist( vectorDist & pbcDist );
        /**
         * @brief Moves protein right
         **/
        void moveProteinRight( lattIndex site );

        /**
         * @brief Moves protein left
         **/
        void moveProteinLeft( lattIndex site );

        /**
         * @brief Moves protein to the Top 
         **/
        void moveProteinUp( lattIndex site );
        /**
         * @brief Moves protein down 
         **/
        void moveProteinDown( lattIndex site );
        virtual void initPushAndPop( int site );
        virtual void pushAndPop( lattIndex site, int &value );
        virtual void push( lattIndex site, int & value );

        /**
         * @brief When particle jumps through PBC, it traces it using internal variable
         **/
        virtual void updatePBCCorrection( lattIndex site, lattIndex newSite );

}; // class MsdEnabledLattEx


