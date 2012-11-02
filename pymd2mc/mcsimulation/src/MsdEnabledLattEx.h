#include "LattExchanger.h"

#include <math.h>
using namespace std;
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
            mPBCCorrection = new vector[ latt->getLatticeSize() ];
        }
        /**
         * TODO: document
         **/
        virtual bool hasMsd() { return true; }

        /**
         * TODO: document
         **/
        virtual void exchangeSites( lattIndex pos1, lattIndex pos2 ) const
        {
            mTracking[pos1] ^= mTracking[pos2];
            mTracking[pos2] ^= mTracking[pos1];
            mTracking[pos1] ^= mTracking[pos2];
            LattExchanger::exchangeSites( pos1, pos2 );
            if ( isPBCJump( pos1, pos2 ) )
            {
                mPBCCorrection[ mpLatt[ pos1 ] ].col = 0;
                mPBCCorrection[ mpLatt[ pos2 ] ].row = 0;
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
                lattIndex startRow = i / mpLatt->getRowSize();
                lattIndex startCol = i % mpLatt->getRowSize();
                lattIndex currRow = mTracking[i] / mpLatt->getRowSize();
                lattIndex currCol = mTracking[i] % mpLatt->getRowSize();
                msd += ( ( startRow - currRow ) * ( startRow - currRow )  + ( startCol - currCol ) * ( startCol - currCol ) );
            }
            msd /= mpLatt->getLatticeSize();
            return msd;
        }

        virtual ~MsdEnabledLattEx()
        {
            delete mPBCCorrection;
            delete mTracking;
        }

    protected: // fields
        lattIndex* mTracking;

        struct vector
        {
            vector() : row( 0 ), col( 0 ) {};
            lattIndex row;
            lattIndex col;
        };
        vector* mPBCCorrection;

    protected: // functions

        virtual bool isPBCJump( lattIndex pos1, lattIndex pos2 )
        {
            return ( abs( pos1 - pos2 ) == 1 || abs( pos1 - pos2 ) == mpLatt->getRowSize() || pos1 - pos2 == -mpLatt->getRowSize() + 1 || pos1 - pos2 == mpLatt->getRowSize() - 1 );
        }

};


