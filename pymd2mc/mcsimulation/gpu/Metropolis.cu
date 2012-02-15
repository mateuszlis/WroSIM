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
void
runTest(int argc, char** argv)
{
    if( cutCheckCmdLineFlag(argc, (const char**)argv, "device") )
        cutilDeviceInit(argc, argv);
    else
        cudaSetDevice( cutGetMaxGflopsDeviceId() );
    // set seed for rand()
    srand(2005);
    const int LATT_ROW_SIZE = 5; // 96;//1790; //894; //5
    const int LATT_ROWS_COUNT = 7; //98;//1792; // 896; //7;
    const long long LATT_SIZE = LATT_ROW_SIZE * LATT_ROWS_COUNT;
    cout << LATT_SIZE << " size of the lattice " << endl;
    TriangularLattice *latt = new TriangularLattice( LATT_SIZE,
        LATT_ROW_SIZE, ( LATT_ROW_SIZE * LATT_ROWS_COUNT ) / 3 ); //400512);

    const int BLOCK_SIZE = 1; //64;//512; //256; //1;
    const int BLOCKS = 5; // 21;//895; //447; //5;
    assert( BLOCKS * BLOCK_SIZE == ( LATT_SIZE / 7 ) );
    const int OUTPUT_FREQ = 300;
    const int STEPS = 1000;
    // allocate device memory for result
    unsigned int size_latt = LATT_SIZE;
    unsigned int mem_size_latt = sizeof(int) * size_latt;
    int* d_latt;
    cutilSafeCall(cudaMalloc((void**) &d_latt, mem_size_latt));

    // allocate host memory for the result
    int* h_latt = latt->getLattice(); //(int*) malloc(mem_size_latt);
    cout << calcSum(h_latt, size_latt ) << " Sum" << endl;
    //printLatt( h_latt, LATT_ROW_SIZE, LATT_ROWS_COUNT );
    cutilSafeCall(cudaMemcpy(d_latt, h_latt, mem_size_latt,
                              cudaMemcpyHostToDevice) );
    float omega = 500;
    float energies[36]; // FIXME:36 is too big
    for ( int i = 0 ; i < 36 ; ++i )
    {
        energies[ i ] = probability( ( i - 10 ) * omega );
        cout << energies[ i ] << endl;
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

    //initRNG<<<grid, threads>>>(d_rngStates, rand() );

    // execute the kernel
    int start = 0;
    for ( ; start < STEPS * 7   ; ++start )
    {
        int random = rand() % 7;

        metropolisStep<<< grid, threads >>>( d_latt, LATT_SIZE, LATT_ROW_SIZE,
            random, d_rngStates, d_energies );
        cudaDeviceSynchronize();

        if ( start % ( OUTPUT_FREQ * 7 )  == 0 )
        {
            
            cutilSafeCall(cudaMemcpy(h_latt, d_latt, mem_size_latt,
                              cudaMemcpyDeviceToHost) );
            //printLatt( h_latt, LATT_ROW_SIZE, LATT_ROWS_COUNT );
            cerr << ( *latt );
            cout << start / 7 << endl;
        }
    }
    
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
    cout << "CUDA res " << calcSum(h_latt, size_latt ) << endl;
    //printLatt( h_latt, LATT_ROW_SIZE, LATT_ROWS_COUNT );

    // check result

    // clean up memory
    //free(h_latt);
    delete latt;
    cutilSafeCall(cudaFree(d_latt));
    cutilSafeCall(cudaFree(d_rngStates));
    cutilSafeCall(cudaFree(d_energies));

    cudaThreadExit();
}

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
            createNeighHist( neighHist );
        }
        if ( mpFNFOutputFile != NULL && i > EQUIB_STEPS )
        {
            ( *mpFNFOutputFile ) << setw( 10 ) << i << "\t" << calcFirstNeighboursFract() << endl;
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
    else
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
