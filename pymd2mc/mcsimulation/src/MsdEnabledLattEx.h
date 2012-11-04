#include "LattExchanger.h"

#include <math.h>
using namespace std;
/**
 * @brief Struct used to represent rows and cols distance between lattice sites
 *
 **/
struct vectorDist
{
    vectorDist() : row( 0 ), col( 0 ) {};
    vectorDist( lattIndex argRow, lattIndex argCol ) : row( argRow ), col( argCol ) {};
    vectorDist operator+( vectorDist added ) { return vectorDist( row + added.row, col + added.col ); }
    vectorDist operator-( vectorDist added ) { return vectorDist( row - added.row, col - added.col ); }
    vectorDist operator+=( vectorDist added ) 
    { 
        col += added.col;
        row += added.row;
        return *this;
    }
    lattIndex squareDisp() { return row*row + col*col; }
    lattIndex row;
    lattIndex col;
}; // struct vectorDist

/**
 * @brief Outputs vector class contents. Used for debug purposes
 **/
ostream &operator<<( ostream &stream, vectorDist & dist );

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
        MsdEnabledLattEx( TriangularLattice* latt );

        /**
         * @brief Since currently it is responsibility of Metropolis class to output statistics, this function is used to determine, whether we have Msd stat enable or not. Should be organized in different way.
         *
         **/
        virtual bool hasMsd() { return true; }

        /**
         * @brief This function is used by Lattice to facilitate tracking of lattice movement
         **/
        virtual void exchangeSites( lattIndex pos1, lattIndex pos2 );

        /**
         * @brief Calculates Mean Square Displacement
         **/
        virtual double calcStat();

        /**
         * @brief Calculates distance between two lattice sites in the means of rows and columns
         *
         * @return number of rows and cols between sites pos1 and pos2
         *
         **/
        vectorDist calcDist( lattIndex pos1, lattIndex pos2 );

        virtual ~MsdEnabledLattEx();

    protected: // fields

        lattIndex* mTracking; ///< track current position of every lipid
        vectorDist* mPBCCorrection; ///< track PBC jumps
        static const bool isNeighbor[3][3]; ///< used to determine if sites are direct neighbors

    protected: // functions

        bool isNotPBCJump( lattIndex pos1, lattIndex pos2 );
        void incDist( vectorDist & pbcDist );

}; // class MsdEnabledLattEx


