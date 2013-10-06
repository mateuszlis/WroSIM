#include "KawasakiProteins.h"
#include "LattExchanger.h"

KawasakiProteins::KawasakiProteins( ProteinTriangularLattice* latt
                        , double omegaAB
                        , double omegaAC
                        , double omegaBC
                        , int T
                        , int equilibSteps
                        , int proteinStepFreq
                        , int proteinStepSize )
    : KawasakiSimulation( latt, omegaAB, T, equilibSteps )
    , mOmegaAC( omegaAC )
    , mOmegaBC( omegaBC )
    , mProteinStepFreq( proteinStepFreq )
    , mProteinStepSize( proteinStepSize )
{
    mpLatt = latt;
}

KawasakiProteins::~KawasakiProteins()
{}


void KawasakiProteins::run( int steps )
{
    for ( int i = 0; i < steps; i++ )
    {
        if ( analysisStep( i ) && mIsSetFrameStream )
        {
            ( *mpFrameStream ) << ( *mpLatt ); //print frame to output
        }
        performAnalysis( i );

        for ( int j = 0; j < mStepSize; ++j )
        {
            metropolisStep();
        }

        if ( i % mProteinStepFreq == 0 )
        {
            for ( int j = 0 ; j < mProteinStepSize ; ++j )
            {
                proteinStep();
            }
        }


    }
    finishAnalysis( steps );

}

void KawasakiProteins::proteinStep()
{
    if ( mpLatt->getProteinCnt() )
    {
        double energyBefore( calcEnergy() );
        ProteinTriangularLattice::lattIndex proteinIndex( rand() % mpLatt->getProteinCnt() );
        lattIndex protein( mpLatt->getProteins()[ proteinIndex ] );
        lattIndex directionSite( rand() % mpLatt->getNeighborsCnt() );
        mpLatt->moveProtein( protein
                           , mpLatt->getNeighbIndex( protein, directionSite  ) );
        double energyAfter( calcEnergy() );
        double p = prob( energyAfter - energyBefore );
        double acceptance = rand() / ( float( RAND_MAX ) + 1 );
        if ( !( p >= 1 or p > ( acceptance ) ) )
        {
           mpLatt->moveProtein( protein
                              , mpLatt->getNeighbIndex( protein, directionSite + 3 ) ); //revoke
        }
    }

};


void KawasakiProteins::metropolisStep()
{
    int pos1 = 0;
    int pos2 = 0;
    if ( unlikely( mpSampler != NULL ) )
    {
        mpSampler( mpLatt, pos1, pos2); //sets pos1 and pos2 by reference
    }
    else
    {
        cout << "Sampler not set" << endl;
        exit( 0 );
    }

    if ( likely( ( *mpLatt )[pos1] != ( *mpLatt )[pos2]
            && isLipid( ( *mpLatt )[pos2] )
            && isLipid( ( *mpLatt )[pos1] ) )  )
    {
        double p = prob( calcEnergyDiff( pos1, pos2) );
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
    int s1Diff = mpLatt->calcOtherLipidNeighbors( pos1 );
    int s2Diff = mpLatt->calcOtherLipidNeighbors( pos2 );
    double s1_protInter = mpLatt->getExchanger()->getProteinInteraction( pos1 );
    double s2_protInter = mpLatt->getExchanger()->getProteinInteraction( pos2 );
    int s1protANeighb = mpLatt->calcNeighbors( pos1, PROTEIN_A );
    int s2protANeighb = mpLatt->calcNeighbors( pos2, PROTEIN_A );

    int s1_after( 7 - ( s1Diff + s1protANeighb ) );
    int s2_after( 7 - ( s2Diff + s2protANeighb ) );
    //int s1prot_after( s2protANeighb );
    //int s2prot_after( s2protANeighb );


    if ( ( *mpLatt )[ pos1 ] == LIPID_A )
    {
        double result = ( s1_after + s2_after - ( s1Diff + s2Diff ) ) * mOmegaAB + ( s2_protInter - s1_protInter ) * mOmegaAC + ( s1_protInter - s2_protInter ) * mOmegaBC;
        return result;
    }
    else
    {
        double result = ( s1_after + s2_after - ( s1Diff + s2Diff ) ) * mOmegaAB + ( s2_protInter - s1_protInter ) * mOmegaBC + ( s1_protInter - s2_protInter ) * mOmegaAC;
        return result;
    }

}

double KawasakiProteins::calcEnergy()
{
    unsigned int ABcount( 0 ), ACcount( 0 ), BCcount( 0 );

    for ( lattIndex i( 0 ) ; i < mpLatt->getLatticeSize() ; ++i )
    {
        for ( lattIndex neighbIndex( 0 ) ; ( neighbIndex < mpLatt->getNeighborsCnt() ) ; ++neighbIndex )
        {
            lattIndex neighb( mpLatt->getNeighbIndex( i, neighbIndex ) );
            if ( ( *mpLatt )[i] != ( *mpLatt )[ neighb ] )
            {
                if ( ( ( *mpLatt )[i] == LIPID_A && ( *mpLatt )[ neighb ] == LIPID_B )
                        || ( ( *mpLatt )[i] == LIPID_B && ( *mpLatt )[ neighb ] == LIPID_A ) )
                {
                    ABcount++;
                }
                if ( ( ( *mpLatt )[i] == LIPID_A && ( *mpLatt )[ neighb ] == PROTEIN_A )
                        || ( ( *mpLatt )[i] == PROTEIN_A && ( *mpLatt )[ neighb ] == LIPID_A ) )
                {
                    ACcount++;
                }
                if ( ( ( *mpLatt )[i] == LIPID_B && ( *mpLatt )[ neighb ] == PROTEIN_A )
                        || ( ( *mpLatt )[i] == PROTEIN_A && ( *mpLatt )[ neighb ] == LIPID_B ) )
                {
                    BCcount++;
                }
            }

        }

    }
    ABcount /= 2; // each pair was counted twice
    BCcount /= 2;
    ACcount /= 2;
    return ABcount * mOmegaAB + ACcount * mOmegaAC + BCcount * mOmegaBC;
}
