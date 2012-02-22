 /** Host code.
 *
 */

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#endif

// std
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <math.h>
#include <assert.h>
using namespace std;

// project related
#include "../TriangularLattice.h"
#include "../Metropolis.h"

// CUDA related
#include <cutil_inline.h>

// kernels
#include "Metropolis_kernel.cu"

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
const float R = 1.986;
const float T = 325;
int* d_latt;
float energies[36]; // FIXME:36 is too big

void runTest(int argc, char** argv);
void printLatt( int *latt, int rowSize, int rowsCount );
long calcSum(int*, int);
void printDiff(float*, float*, int, int);
float probability( float dG );

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
//int main(int argc, char** argv)
//{
//    runTest(argc, argv);
//
//    cutilExit(argc, argv);
//}

////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////

long calcSum(int* data, int size)
{
    long sum( 0 );
    for (int i = 0; i < size; ++i)
        sum += data[i]; 
    return sum;
}


void printLatt( int *latt, int rowSize, int rowsCount )
{
    for ( int i = 0 ; i < rowsCount ; ++i )
    {
        cout << setw( 2 * ( i + 1 ) ) << "  " ;
        for( int j = 0 ; j < rowSize; ++j )
        {   
            cout << setw(2) << latt[ i * rowSize + j ] << "  ";
        }
        cout << endl;
    }
}
float probability( float dG )
{
    return exp( -dG / ( R * T ) );
}

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
Metropolis::Metropolis( TriangularLattice* latt, double omegaAB, int T ) 
    : mOmegaAB( omegaAB )
      , mOutputFreq( 100 )
      , mpNeighOutputFile( NULL )
      , mpFNFOutputFile( NULL )
      , mpFrameStream( NULL )
      , mpStatusStream( NULL )
      , mIsSetFrameStream( false )
      , mIsSetNeighOutputFile( false )
      , mIsSetStatusStream( false )
      , mpLatt( latt )
      , mpHistArr( NULL )
      , T( T )

{

    mStepSize = 7;
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

void Metropolis::setOmegaAB( double omegaAB )
{
    mOmegaAB = omegaAB;
}

unsigned int getBlockSize( TriangularLattice* pLatt )
{
    if ( pLatt->getLatticeSize() % 7 != 0 )
    {
        throw std::exception(); // must be divisible by 7
        // TODO: create exception chierarchy
    }
    unsigned int expectedSize = pLatt->getLatticeSize() / 7;
    unsigned int maxBlocks = 512;
    while ( expectedSize % maxBlocks != 0 )
    {
        maxBlocks /= 2;
    }
    return maxBlocks;
}

void Metropolis::run( int steps )
{
    cudaSetDevice( cutGetMaxGflopsDeviceId() );
    // set seed for rand()

    const int BLOCK_SIZE = getBlockSize( mpLatt ); //64;//512; //256; //1;
    const int BLOCKS = ( mpLatt->getLatticeSize() / 7 ) / BLOCK_SIZE; // 21;//895; //447; //5;
    assert( BLOCKS * BLOCK_SIZE == ( mpLatt->getLatticeSize() / 7 ) );
    // allocate device memory for result
    unsigned int mem_size_latt = sizeof( int ) * mpLatt->getLatticeSize();
    cutilSafeCall(cudaMalloc((void**) &d_latt, mem_size_latt));

    // allocate host memory for the result
    int* h_latt = mpLatt->getLattice(); //(int*) malloc(mem_size_latt);
    cout << calcSum(h_latt, mpLatt->getLatticeSize() ) << " Sum" << endl;
    cutilSafeCall(cudaMemcpy(d_latt, h_latt, mem_size_latt,
                              cudaMemcpyHostToDevice) );
    for ( int i = 0 ; i < 36 ; ++i )
    {
        energies[ i ] = probability( ( i - 10 ) * mOmegaAB );
    }

    float* d_energies;
    int mem_size_energies = sizeof( float ) * 36;
    cutilSafeCall(cudaMalloc((void**) &d_energies, mem_size_energies));
    cutilSafeCall(cudaMemcpy(d_energies, energies, mem_size_energies,
                              cudaMemcpyHostToDevice) );
    // create and start timer
    unsigned int timer = 0;
    cutilCheckError(cutCreateTimer(&timer));
    cutilCheckError(cutStartTimer(timer));


    // setup execution parameters
    dim3 threads( BLOCK_SIZE );
    dim3 grid( BLOCKS );

    // Allocate memory for RNG states
    curandState *d_rngStates = 0;
    cutilSafeCall( cudaMalloc((void **)&d_rngStates, grid.x * threads.x *
        sizeof(curandState)) );
    // Initialise RNG
    srand( (unsigned)time(0) );

    initRNG<<<grid, threads>>>(d_rngStates, rand() );

    long long neighHist[7] = { 0, 0, 0, 0, 0, 0, 0 };
    for ( int i = 0; i < steps; i++ )
    {
        for ( int j = 0; j < mStepSize; j++ )
        {
            int random = rand() % 7;

            metropolisStep_kernel<<< grid, threads >>>( d_latt,
            mpLatt->getLatticeSize(), mpLatt->getRowSize(),
                random, d_rngStates, d_energies );
            cudaDeviceSynchronize();
        }

        if ( i % mOutputFreq == 0 && mIsSetFrameStream )
        {
            cutilSafeCall(cudaMemcpy(h_latt, d_latt, mem_size_latt,
                              cudaMemcpyDeviceToHost) );
            ( *mpFrameStream ) << ( *mpLatt ); //print frame to output
        }
        if ( mIsSetNeighOutputFile )
        {
            createNeighHist( neighHist );
        }
        if ( mpFNFOutputFile != NULL )
        {
            ( *mpFNFOutputFile ) << setw( 10 ) << i << "\t" << calcFirstNeighboursFract() << endl;
        }

        if ( mIsSetStatusStream )
        {
            ( *mpStatusStream ) << "\r" << i; //print status message
            ( *mpStatusStream ).flush();
        }

    }
    if ( mIsSetNeighOutputFile )
        for ( int i = 0; i < 7; i++ )
        {
            double freq = static_cast< double > ( neighHist[i] ) / ( mpLatt->getLatticeSize() * ( steps ) );
            ( *mpNeighOutputFile ) << i 
                << " " 
                <<  freq 
                << "\t" << ( neighHist[ i ] ) 
                << endl;
        }
    // execute the kernel

    
    // check if kernel execution generated and error
    cutilCheckMsg("Kernel execution failed");
    cutilCheckMsg("cudaThreadSynchronize execution failed");
    // copy result from device to host
    cutilSafeCall(cudaMemcpy(h_latt, d_latt, mem_size_latt,
                              cudaMemcpyDeviceToHost) );

    // stop and destroy timer
    cutilCheckError(cutStopTimer(timer));
    printf("Processing time: %f (ms) \n", cutGetTimerValue(timer));
    cutilCheckError(cutDeleteTimer(timer));
    // compute reference solution
    //computeGold(reference, h_A, h_B, HA, WA, WB);
    cout << "CUDA res " << calcSum(h_latt, mpLatt->getLatticeSize() ) << endl;


    // clean up memory
    cutilSafeCall(cudaFree(d_latt));
    cutilSafeCall(cudaFree(d_rngStates));
    cutilSafeCall(cudaFree(d_energies));

    cudaThreadExit();

}

void Metropolis::setSampler( Sampler s )
{
    ;
}

//private functions

double Metropolis::prob( double dG )
{
    return exp( -dG / ( R * T ));
}

double Metropolis::calcEnergyDiff( int pos1, int pos2 )
{
    return 0;
}

void Metropolis::metropolisStep()
{
    ;
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
