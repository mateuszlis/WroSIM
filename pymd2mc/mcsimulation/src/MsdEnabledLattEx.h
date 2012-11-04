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
    lattIndex squareDisp() { return row*row + col*col; }
    lattIndex row;
    lattIndex col;
};
ostream &operator<<( ostream &stream, vectorDist & dist );

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
            if ( !isNotPBCJump( pos1, pos2 ) )
            {
                vectorDist lDist( calcDist( pos1, pos2 ) );
                mPBCCorrection[ mTracking[ pos1 ] ] = calcDist( pos1, pos2 );
                incDist( mPBCCorrection[ mTracking[ pos1 ] ] );
                mPBCCorrection[ mTracking[ pos2 ] ] = calcDist( pos2, pos1 );
                incDist( mPBCCorrection[ mTracking[ pos2 ] ] );
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
                msd += ( calcDist( i, mTracking[i] ) - mPBCCorrection[ mTracking[i] ] ).squareDisp();
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

        static const bool isNeighbor[3][3];

    protected: // functions

        bool isNotPBCJump( lattIndex pos1, lattIndex pos2 ) 
        {
            vectorDist dist( calcDist( pos1, pos2 ) );
            if ( dist.squareDisp() <= 2 )
            {
                return isNeighbor[ dist.col - 1 ][ dist.row - 1 ];
            }
            return false;
                    
        }

        void incDist( vectorDist & pbcDist )
        {
            if ( pbcDist.col > 1 )
            {
                pbcDist.col++;
            }
            if ( pbcDist.col < -1 )
            {
                pbcDist.col--;
            }
            if ( pbcDist.row > 1 )
            {
                pbcDist.row++;
            }
            if ( pbcDist.row < -1 )
            {
                pbcDist.row--;
            }
        }

}; // class MsdEnabledLattEx


