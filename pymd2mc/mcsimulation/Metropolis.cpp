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

//Metropolis public functions

Metropolis::Metropolis( TriangularLattice* latt, double omegaAB ) :
    mOmegaAB( omegaAB), mOutputFreq( 100), mIsSetFrameStream( false), mIsSetNeighOutputFile( false),
            mIsSetStatusStream( false), mpLatt( latt)
{

    mStepSize = ( latt->getRowSize() );
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

        if ( i % mOutputFreq == 0 && mIsSetFrameStream )
        {
            ( *mpFrameStream ) << ( *mpLatt ); //print frame to output
        }
        if ( mIsSetNeighOutputFile && i > EQUIB_STEPS )
        {
            createNeighHist( neighHist);
        }
        if ( mIsSetStatusStream )
        {
            ( *mpStatusStream ) << "\r" << i; //print status message
        }

    }
    if ( mIsSetNeighOutputFile )
        for ( int i = 0; i < 7; i++ )
        {
            ( *mpNeighOutputFile ) << i << " " << neighHist[i] / ( steps - EQUIB_STEPS ) << endl;
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
}

//private functions

double Metropolis::prob( double dG )
{
    return exp( -dG / ( R * T ));
}

double Metropolis::calcEnergyDiff( int pos1, int pos2 )
{
    int diff1 = ( 6 - mpLatt->simNeighbCount( pos1) ) + ( 6 - mpLatt->simNeighbCount( pos2) );
    mpLatt->exchangeSites( pos1, pos2);
    int diff2 = ( 6 - mpLatt->simNeighbCount( pos1) ) + ( 6 - mpLatt->simNeighbCount( pos2) );
    mpLatt->exchangeSites( pos1, pos2);
    return ( diff2 - diff1 ) * mOmegaAB;
}

void Metropolis::metropolisStep()
{
    int pos1 = 0;
    int pos2 = 0;
    mpSampler( mpLatt, pos1, pos2); //sets pos1 and pos2 by reference

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
        if ( ( *mpLatt )[i] )
            histArr[mpLatt->simNeighbCount( i)] += 1;
    }
}
Metropolis::~Metropolis()
{
    // TODO Auto-generated destructor stub
}

