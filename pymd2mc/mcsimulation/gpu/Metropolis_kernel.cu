#ifndef _MATRIXMUL_KERNEL_H_
#define _MATRIXMUL_KERNEL_H_

#include <stdio.h>
#include <curand_kernel.h>

//1790 
//894

// RNG init kernel
__global__ void initRNG(curandState * const rngStates,
                        const unsigned int seed)
{
    // Determine thread ID
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // Initialise the RNG
    curand_init(seed, tid, 0, &rngStates[tid]);
}

__device__ inline void exchangeSites(int* latt, int s1, int s2 )
 {
    latt[s1] ^= latt[s2];
    latt[s2] ^= latt[s1];
    latt[s1] ^= latt[s2];
 }

__device__ inline void getPosition( int & candidatePos, int lattSize )
{ 
    if ( candidatePos >= lattSize )
        candidatePos -= lattSize;
    else if ( candidatePos < 0 )
        candidatePos += lattSize;
}
        
__device__ inline int calcDiffNeighbors(int* latt
    , int lattSize
    , int lattRowSize
    , int site
    , int neighbors[] )
{
    int result = 0;
    for ( int i = 0; i < 6 ; ++i )
    {
        int neighbor = site + neighbors[ i ];
        getPosition( neighbor, lattSize ); 
        result += (int) ( latt[ neighbor ] != latt[ site ] );
    }
    return result;

}

__global__ void
metropolisStep_kernel( int * latt, int lattSize, int lattRowSize, int start,
        curandState * const rngStates, float * energies )
{
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Initialise RNG for this thread
    curandState & localState = ( rngStates[ tid ] );

    int pos1 = ( tid ) * 7 + start;
    int neighbors[] = { 1, -1, lattRowSize, -lattRowSize, lattRowSize - 1, -lattRowSize + 1 };

    unsigned int randomNeighbor = neighbors[ ( curand( &localState ) ) % 6];
    int pos2 = ( pos1 + randomNeighbor ); //neighbors[ randomNeighbor ] ); //3207680; //801024 ;
    getPosition( pos2, lattSize );

    if ( latt[ pos1 ] != latt[ pos2 ] )
    {
        int diffNeighb1( calcDiffNeighbors( latt, lattSize, lattRowSize, pos1,
        neighbors ) )
            , diffNeigh2( calcDiffNeighbors( latt, lattSize, lattRowSize, pos2,
            neighbors ) );
        int diffBefore = diffNeighb1 + diffNeigh2;
        int diffAfter = 14 - ( diffNeighb1 + diffNeigh2 );
        float prob = energies[ diffAfter - diffBefore + 10 ];

        float acceptance = curand_uniform( &localState );

        if ( prob > acceptance )
        {
            exchangeSites( latt, pos1, pos2 );
        }
    }
      //rngStates[ tid ] = localState;
}

#endif // #ifndef _MATRIXMUL_KERNEL_H_
