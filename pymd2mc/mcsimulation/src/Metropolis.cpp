/*
 * Metropolis.cpp
 *
 *  Created on: 20-06-2011
 *      Author: lisu
 */

#include "Metropolis.h"
#include "math.h"

//Metropolis public functions

Metropolis::Metropolis( TriangularLattice* latt, double omegaAB, int T, int equilibSteps) 
    : mEquilibSteps( equilibSteps )
      , mT( T )
      , mOutputFreq( 100 )
      , mpNeighOutputFile( NULL )
      , mpFNFOutputFile( NULL )
      , mpFrameStream( NULL )
      , mpStatusStream( NULL )
      , mpClusterStream( NULL )
      , mpMsdOutputFile( NULL )
      , mIsSetFrameStream( false )
      , mIsSetNeighOutputFile( false )
      , mIsSetStatusStream( false )
      , mIsSetClusterStream( false )
      , mpLatt( latt )
      , mpHistArr( NULL )
      , mOmegaAB( omegaAB )

{
    mStepSize = ( latt->getLatticeSize() );
    for ( unsigned int i = 0 ; i < LATTICE_NEIGHB_COUNT ; ++i )
    {
        mNeighHist[ i ] = 0;
    }
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

void Metropolis::setMsdOutput( ostream &MsdOutput )
{
    mpMsdOutputFile = &MsdOutput;
}

void Metropolis::setOmegaAB( double omegaAB )
{
    mOmegaAB = omegaAB;
}


//private functions

void Metropolis::performAnalysis( int stepNum )
{
    if ( analysisStep( stepNum ) && mIsSetNeighOutputFile )
    {
        createNeighHist( mNeighHist );
    }
    if ( analysisStep( stepNum ) && mpFNFOutputFile != NULL )
    {
        ( *mpFNFOutputFile ) << setw( 10 ) << stepNum << "\t" << calcFirstNeighboursFract() << endl;
    }
    if ( analysisStep( stepNum ) && mpClusterStream != NULL )
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
        ( *mpStatusStream ) << "\r" << stepNum; //print status message
        ( *mpStatusStream ).flush();

    }
    if ( mpLatt->getExchanger()->hasMsd() && mpMsdOutputFile )
    {
        ( *mpMsdOutputFile ) << setw( 10 ) << stepNum << "\t" << mpLatt->getExchanger()->calcStat() << endl;
    }
}

void Metropolis::finishAnalysis( int steps )
{
    if ( mIsSetNeighOutputFile )
        for ( int i = 0; i < 7; i++ )
        {
            double freq = static_cast< double > ( mNeighHist[i] ) / ( mpLatt->getLatticeSize() * ( steps - mEquilibSteps ) );
            ( *mpNeighOutputFile ) << i 
                << " " 
                <<  freq 
                << "\t" << ( mNeighHist[ i ] ) 
                << endl;
        }
}
double Metropolis::prob( double dG )
{
    return exp( -dG / ( R * mT ));
}

double Metropolis::calcEnergyDiff( int pos1, int pos2 )
{
    int s1Diff = 6 - mpLatt->simNeighbCount( pos1 );
    int s2Diff = 6 - mpLatt->simNeighbCount( pos2 );
    
    int diff1 = s1Diff + s2Diff;
    int diff2 = 14 - ( s1Diff + s2Diff ) ;
    return ( diff2 - diff1 ) * mOmegaAB;
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

