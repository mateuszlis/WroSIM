#include "LattExchanger.h"

#include <math.h>
using namespace std;
class MsdEnabledLattEx : public LattExchanger
{
    public: // typedefs

    public: // functions
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
        }

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

        ~MsdEnabledLattEx()
        {
            delete mTracking;
        }

    protected: // fields
        lattIndex* mTracking;

};


