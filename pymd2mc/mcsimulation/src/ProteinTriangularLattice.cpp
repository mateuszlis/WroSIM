#include "ProteinTriangularLattice.h"

#include "LattExchanger.h"
#include <limits>

ProteinTriangularLattice::ProteinTriangularLattice( lattIndex latticeSize
                                , lattIndex rowSize
                                , lattIndex firstTypeParticlesCnt
                                , lattIndex proteinCnt
                                , bool distributeRandomly )
    : TriangularLattice( latticeSize, rowSize, firstTypeParticlesCnt, distributeRandomly )
      , mProteinCnt( proteinCnt )
{
    if( distributeRandomly )
        distributeProteinsRandomly();
    else
        distributeProteins();
    mpExchanger->setProteins( mProteins, mProteinCnt );
    TriangularLattice::mNonLipidMembersCount = 7 * proteinCnt;
}

void ProteinTriangularLattice::distributeProteinsRandomly()
{
    if ( mProteinCnt )
    {
        mProteins = new lattIndex[ mProteinCnt ];
        for ( int i = 0; i < mProteinCnt; ++i )
        {
            lattIndex pos = rand() % getLatticeSize(); //FIXME: if RAND_MAX is too small this causes bad distribution
            while ( !isFree( pos ) )
            {
                    pos = rand() % getLatticeSize();
            };
            mProteins[i] = pos;
            putProtein( pos );
        }
    }
}

void ProteinTriangularLattice::distributeProteins()
{
    if ( mProteinCnt )
    {
        mProteins = new lattIndex[ mProteinCnt ];
        int counter( 0 );
        for ( lattIndex pos = 0 ; pos < getLatticeSize() && counter < mProteinCnt ; ++pos )
        {
            if ( isFree( pos ) )
            {
                putProtein( pos );
                mProteins[ counter ] = pos;
                ++counter;
            }
        }
    }
}

void ProteinTriangularLattice::moveProtein( lattIndex site, lattIndex destination )
{
    mpExchanger->moveProtein( site, destination );
}

ProteinTriangularLattice::~ProteinTriangularLattice()
{
    if ( mProteinCnt )
        delete mProteins;
}


// helpers
//

bool ProteinTriangularLattice::isFree( lattIndex pos )
{
    static const int safeDistance( 9 ); // squared!
    bool isFree = true;
    // check if only lipids A will be overwritten
    for ( unsigned int i = 0 ; i < mProteinSize ; ++i )
    {
        isFree = isFree && ( mpLattice[ cutToPos( pos + mpExchanger->mProteinSites[i] ) ] == LIPID_A );
    }

    // check for collisions
    for ( int site = 0 ; site < mLatticeSize ; ++site )
    {
        if ( ( mpExchanger->calcDist( pos, site ) ).squareDisp() <= safeDistance )
        {
            isFree = isFree && !isProtein( mpLattice[ site ] );
        }
    }
    return isFree;
};


void ProteinTriangularLattice::putProtein( lattIndex pos, LATTICE_FIELD_NAMES protein )
{
    for ( unsigned int i = 0 ; i < mProteinSize ; ++i )
    {
        mpLattice[ cutToPos( pos + mpExchanger->mProteinSites[i] ) ] = protein;
    }
};



