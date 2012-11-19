#include "KawasakiProteins.h"

KawasakiProteins::KawasakiProteins( ProteinTriangularLattice* latt, double omegaAB, int T, int equilibSteps )
    : KawasakiSimulation( latt, omegaAB, T, equilibSteps )
{
    mpLatt = latt;
}

KawasakiProteins::~KawasakiProteins()
{}


void KawasakiProteins::run( int steps )
{
    for ( int i = 0; i < steps; i++ )
    {
        for ( int j = 0; j < mStepSize; j++ )
        {
            metropolisStep();
        }

        proteinStep();

        if ( analysisStep( i ) && mIsSetFrameStream )
        {
            ( *mpFrameStream ) << ( *mpLatt ); //print frame to output
        }
        performAnalysis( i );

    }
    finishAnalysis( steps );

}

void KawasakiProteins::proteinStep()
{
    if ( mpLatt->getProteinCnt() )
    {
        ProteinTriangularLattice::lattIndex protein( rand() % mpLatt->getProteinCnt() );
        mpLatt->moveProtein( mpLatt->getProteins()[protein], mpLatt->getProteins()[protein] + 1 ); //FIXME +1 is just a quick thing
    }

};


void KawasakiProteins::metropolisStep()
{
    int pos1 = 0;
    int pos2 = 0;
    if ( mpSampler != NULL )
    {
        mpSampler( mpLatt, pos1, pos2); //sets pos1 and pos2 by reference
    }
    else
    {
        cout << "Sampler not set" << endl;
        exit( 0 );
    }

    if ( ( *mpLatt )[pos1] != ( *mpLatt )[pos2] && isLipid( ( *mpLatt )[pos2] ) && isLipid( ( *mpLatt )[pos1] )  )
    {
        double p = prob( calcEnergyDiff( pos1, pos2));
        double acceptance = rand() / ( float( RAND_MAX ) + 1 );
        if ( p >= 1 or p > ( acceptance ) )
        {
            mpLatt->exchangeSites( pos1, pos2);
        }
    }
}

bool KawasakiProteins::isLipid( lattMember site )
{
    return ( site == LIPID_A || site == LIPID_B );
}

double KawasakiProteins::calcEnergyDiff( int pos1, int pos2 )
{
    static const double mProtein_A = 300;
    static const double mProtein_B = -200;
    int s1Diff = 6 - mpLatt->simNeighbCount( pos1 );
    int s2Diff = 6 - mpLatt->simNeighbCount( pos2 );
    
    int diff1 = s1Diff + s2Diff;
    int diff2 = 14 - ( s1Diff + s2Diff ) ;
    //mpLatt->getProteins()
    return ( diff2 - diff1 ) * mOmegaAB;
}
