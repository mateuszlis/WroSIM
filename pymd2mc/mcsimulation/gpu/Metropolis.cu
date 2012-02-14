 /** Host code.
 *
 */

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#endif

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "../TriangularLattice.h"
#include <vector>
#include <algorithm>
#include <ctime>
#include <math.h>
#include <assert.h>
using namespace std;

// includes, project
#include <cutil_inline.h>

// includes, kernels
#include "Metropolis_kernel.cu"

class Metropolis
{
}
;
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
int main(int argc, char** argv)
{
    runTest(argc, argv);

    cutilExit(argc, argv);
}

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
