#include "KawasakiSimulation.h"

//Non object functions
void almeidaSampler( TriangularLattice *latt, int &pos1, int &pos2 )
{
    pos1 = rand() % latt->getLatticeSize();
    pos2 = rand() % latt->getLatticeSize();
}
void kawasakiSampler( TriangularLattice *latt, int &pos1, int &pos2 )
{
    pos1 = rand() % latt->getLatticeSize();
    pos2 = latt->getNeighbIndex( pos1, rand() % latt->getNeighborsCnt());
}

/**
 * @brief This is experimental feature that is here to test ideas connected with
 * GPU implementation of Kawasaki sampler.
 */
void massiveParallelKawasakiSampler( TriangularLattice *latt, int &pos1, int &pos2 )
{
    static long step = 0;
    static int start = 0;
    if( step == latt->getLatticeSize() / 7 || step == 0 ) // As every step involves 7 sites
        // we perform latticeSize/7 steps to involve all sites in lattice
    {
        pos1 = rand() % 7; // chose initial point for sampling (note that there are only 7 possibilities
        start = pos1;
        step = 0;
    }
    else
    {
        pos1 = ( 7 * step + start ) % latt->getLatticeSize();
    }
    pos2 = latt->getNeighbIndex( pos1, rand() % latt->getNeighborsCnt()); // randomly chose second exchange site
    step++;
}

// KawasakiSimulation public functions

KawasakiSimulation::KawasakiSimulation( TriangularLattice* latt, double omegaAB, int T, int equilibSteps )
    : Metropolis( latt, omegaAB, T, equilibSteps )
{}

KawasakiSimulation::~KawasakiSimulation()
{}


void KawasakiSimulation::setSampler( Sampler s )
{
    if ( s == Kawasaki )
    {
        mpSampler = kawasakiSampler;
    }
    if ( s == Almeida )
    {
        mpSampler = almeidaSampler;
    }
    if ( s == MassiveParallelKawasaki )
    {
        cout << "Massive parallel Kawasaki Sampler test!" << endl;
        mpSampler = massiveParallelKawasakiSampler;
    }
}
void KawasakiSimulation::run( int steps )
{
    for ( int i = 0; i < steps; i++ )
    {
        for ( int j = 0; j < mStepSize; j++ )
        {
            metropolisStep();
        }

        if ( analysisStep( i ) && mIsSetFrameStream )
        {
            ( *mpFrameStream ) << ( *mpLatt ); //print frame to output
        }
        performAnalysis( i );

    }
    finishAnalysis( steps );

}

void KawasakiSimulation::metropolisStep()
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

    if ( ( *mpLatt )[pos1] != ( *mpLatt )[pos2] )
    {
        double p = prob( calcEnergyDiff( pos1, pos2));
        double acceptance = rand() / ( float( RAND_MAX ) + 1 );
        if ( p >= 1 or p > ( acceptance ) )
        {
            mpLatt->exchangeSites( pos1, pos2);
        }
    }
}
