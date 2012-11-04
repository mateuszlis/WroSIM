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
    lattIndex squareDisp() { return row*row + col*col; }
    lattIndex row;
    lattIndex col;
};

/**
 * TODO: document
 *
 **/
class MsdEnabledLattEx : public LattExchanger
{
    public: // typedefs

    public: // functions
        //FIXME: move implementations to cpp file!
        /**
         * TODO: document
         **/
        MsdEnabledLattEx( TriangularLattice* latt ) 
            : LattExchanger( latt ) 
        {
            mTracking = new lattIndex[ latt->getLatticeSize() ];
            for ( lattIndex i = 0 ; i < latt->getLatticeSize() ; ++i )
            {
                mTracking[i] = i;
            }
            mPBCCorrection = new vectorDist[ latt->getLatticeSize() ];
        }
        /**
         * TODO: document
         **/
        virtual bool hasMsd() { return true; }

        /**
         * TODO: document
         **/
        virtual void exchangeSites( lattIndex pos1, lattIndex pos2 ) 
        {
            mTracking[pos1] ^= mTracking[pos2];
            mTracking[pos2] ^= mTracking[pos1];
            mTracking[pos1] ^= mTracking[pos2];
            LattExchanger::exchangeSites( pos1, pos2 );
            if ( isPBCJump( pos1, pos2 ) )
            {
                mPBCCorrection[ mTracking[ pos1 ] ].col = 0;
                mPBCCorrection[ mTracking[ pos2 ] ].row = 0;
            }
        }


        /**
         * TODO: document
         **/
        virtual double calcStat()
        {
            double msd( 0 );
            for ( lattIndex i = 0 ; i < mpLatt->getLatticeSize() ; ++i )
            {
                msd += calcDist( i, mTracking[i] ).squareDisp();
            }
            msd /= mpLatt->getLatticeSize();
            return msd;
        }

        /**
         * TODO: document
         *
         **/
        vectorDist calcDist( lattIndex pos1, lattIndex pos2 )
        {
                lattIndex startRow = pos1 / mpLatt->getRowSize();
                lattIndex startCol = pos1 % mpLatt->getRowSize();
                lattIndex endRow = pos2 / mpLatt->getRowSize();
                lattIndex endCol = pos2 % mpLatt->getRowSize();
                return vectorDist( startRow - endRow, startCol - endCol );

        }

        virtual ~MsdEnabledLattEx()
        {
            delete[] mPBCCorrection;
            delete[] mTracking;
        }

    protected: // fields
        lattIndex* mTracking;

        vectorDist* mPBCCorrection;

    protected: // functions

        virtual bool isPBCJump( lattIndex pos1, lattIndex pos2 ) 
        {
            return ( abs( pos1 - pos2 ) == 1 || abs( pos1 - pos2 ) == mpLatt->getRowSize() || pos1 - pos2 == -mpLatt->getRowSize() + 1 || pos1 - pos2 == mpLatt->getRowSize() - 1 );
        }

};


