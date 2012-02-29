/*
 * Metropolis.cpp
 *
 *  Created on: 20-06-2011
 *      Author: lisu
 */

#include "Metropolis.h"
#include "math.h"

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
//Metropolis public functions

Metropolis::Metropolis( TriangularLattice* latt, double omegaAB, int T ) 
    : mOmegaAB( omegaAB )
      , mOutputFreq( 100 )
      , mpNeighOutputFile( NULL )
      , mpFNFOutputFile( NULL )
      , mpFrameStream( NULL )
      , mpStatusStream( NULL )
      , mpClusterStream( NULL )
      , mIsSetFrameStream( false )
      , mIsSetNeighOutputFile( false )
      , mIsSetStatusStream( false )
      , mpLatt( latt )
      , mpHistArr( NULL )
      , T( T )

{

    mStepSize = ( latt->getLatticeSize() );
    setSampler( Kawasaki);

}

void Metropolis::setOutputFreq( int freq )
{
    mOutputFreq = freq;
}

void Metropolis::setNeighOutput( ostream &neighOutputFile )
{
    mIsSetNeighOutputFile = true;
    mpNeighOutputFile = &neighOutputFile;
}

void Metropolis::setOutput( ostream &frameStream )
{
    mIsSetFrameStream = true;
    mpFrameStream = &frameStream;
}
void Metropolis::setStatus( ostream &statusStream )
{
    mIsSetStatusStream = true;
    mpStatusStream = &statusStream;
}
void Metropolis::setFNFStream( ostream &fnfStream )
{
    mpFNFOutputFile = &fnfStream;
}
void Metropolis::setClusterStream( ostream &clusterStream )
{
    mpClusterStream = &clusterStream;
}

void Metropolis::setOmegaAB( double omegaAB )
{
    mOmegaAB = omegaAB;
}

void Metropolis::run( int steps )
{
    long long neighHist[7] = { 0, 0, 0, 0, 0, 0, 0 };
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
        if ( analysisStep( i ) && mIsSetNeighOutputFile )
        {
            createNeighHist( neighHist );
        }
        if ( analysisStep( i ) && mpFNFOutputFile != NULL )
        {
            ( *mpFNFOutputFile ) << setw( 10 ) << i << "\t" << calcFirstNeighboursFract() << endl;
        }
        if ( analysisStep( i ) && mpClusterStream != NULL )
        {
            TriangularLattice::clustersMap map;
            mpLatt->calculateClusters( map );
            int sum = 0;
            for( TriangularLattice::clustersMap::const_iterator it = map.begin() ; it != map.end() ; ++it )
            {
                sum += ( *it ).second * ( *it ).first;
                ( *mpClusterStream ) << ( *it ).first << "\t" << ( *it ).second << std::endl;
            }
            cout << sum << " " << map.size()<<  std::endl;
        }
        if ( mIsSetStatusStream )
        {
            ( *mpStatusStream ) << "\r" << i; //print status message
        }

    }
    if ( mIsSetNeighOutputFile )
        for ( int i = 0; i < 7; i++ )
        {
            double freq = static_cast< double > ( neighHist[i] ) / ( mpLatt->getLatticeSize() * ( steps - EQUIB_STEPS ) );
            ( *mpNeighOutputFile ) << i 
                << " " 
                <<  freq 
                << "\t" << ( neighHist[ i ] ) 
                << endl;
        }

}

void Metropolis::setSampler( Sampler s )
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

//private functions

double Metropolis::prob( double dG )
{
    return exp( -dG / ( R * T ));
}

double Metropolis::calcEnergyDiff( int pos1, int pos2 )
{
    int s1Diff = 6 - mpLatt->simNeighbCount( pos1 );
    int s2Diff = 6 - mpLatt->simNeighbCount( pos2 );
    
    int diff1 = s1Diff + s2Diff;
    int diff2 = 14 - ( s1Diff + s2Diff ) ;
    return ( diff2 - diff1 ) * mOmegaAB;
}

void Metropolis::metropolisStep()
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
        double acceptance = rand() / ( float( RAND_MAX) + 1 );
        if ( p >= 1 or p > ( acceptance ) )
        {
            //cout << "DONE MOVE " << pos1 << " " << pos2 << endl;
            mpLatt->exchangeSites( pos1, pos2);
        }
    }
}
void Metropolis::createNeighHist( long long *histArr )
{
    for ( int i = 0; i < mpLatt->getLatticeSize(); i++ )
    {
        histArr[ mpLatt->simNeighbCount( i ) ] += 1;
    }
}
double Metropolis::calcFirstNeighboursFract() 
{
    int simFirstNeighbCount = 0;
    for ( int i = 0 ; i < mpLatt->getLatticeSize() ; i++ )
    {
        int pos = mpLatt->getNeighbIndex( i, rand() % mpLatt->getNeighborsCnt());
        if ( ( *mpLatt )[pos] == ( *mpLatt )[i] )
        {
            simFirstNeighbCount++;
        }
    }
    return static_cast< double >( simFirstNeighbCount ) / mpLatt->getLatticeSize();
}



Metropolis::~Metropolis()
{
    // TODO Auto-generated destructor stub
}

