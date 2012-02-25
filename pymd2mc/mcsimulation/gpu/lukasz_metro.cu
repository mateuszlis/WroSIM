/** Host code.
*
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "triangularLattice.h"
#include <vector>
#include <algorithm>
#include <ctime>
#include <math.h>
#include <assert.h>
#include <sstream>

// includes, kernels
#include "curand_kernel.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <nppi.h>
#include <npp.h>
#include <FreeImage.h>
#include <ImageIO.h>

using namespace std;
using namespace npp;

 void mySaveImage(const std::string & rFileName, const ImageCPU_8u_C1 & rImage)
{
			// create the result image storage using FreeImage so we can easily 
			// save
	FIBITMAP * pResultBitmap = FreeImage_Allocate(rImage.width(), rImage.height(), 8 /* bits per pixel */);
	RGBQUAD * palettePtr = FreeImage_GetPalette(pResultBitmap);

	for(int i = 0; i < 256; i++)
	{
		(*palettePtr).rgbReserved = i;
		(*palettePtr).rgbRed = i;
		(*palettePtr).rgbGreen = i;
		(*palettePtr).rgbBlue = i;
		palettePtr++;
	}

	NPP_ASSERT_NOT_NULL(pResultBitmap);
	unsigned int nDstPitch   = FreeImage_GetPitch(pResultBitmap);
	Npp8u * pDstLine = FreeImage_GetBits(pResultBitmap) + nDstPitch * (rImage.height()-1);
	const Npp8u * pSrcLine = rImage.data();
	unsigned int nSrcPitch = rImage.pitch();


	for (size_t iLine = 0; iLine < rImage.height(); ++iLine)
	{
		memcpy(pDstLine, pSrcLine, rImage.width() * sizeof(Npp8u));
		pSrcLine += nSrcPitch;
		pDstLine -= nDstPitch;
	}

		
	// now save the result image
	bool bSuccess;
	FREE_IMAGE_FORMAT eFormat = FIF_PNG;//FreeImage_GetFileType(rFileName.c_str());
	int flags = eFormat==FIF_PNG ? PNG_DEFAULT : eFormat==FIF_BMP ? BMP_DEFAULT : 0;
	bSuccess = FreeImage_Save(eFormat, pResultBitmap, rFileName.c_str(), flags) == TRUE;
	NPP_ASSERT_MSG(bSuccess, "Failed to save result image.");
}

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

__device__ inline void exchangeSites(unsigned char * latt, int s1, int s2 )
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

__device__ inline int calcDiffNeighbors(unsigned char* latt
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

__global__ void metropolisStep( unsigned char * latt, int lattSize, int lattRowSize, int start, curandState * const rngStates, float * energies )
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
		int diffNeighb1( calcDiffNeighbors( latt, lattSize, lattRowSize, pos1, neighbors )),
			diffNeigh2( calcDiffNeighbors( latt, lattSize, lattRowSize, pos2, neighbors ));
		int diffBefore = diffNeighb1 + diffNeigh2;
		int diffAfter = 14 - ( diffNeighb1 + diffNeigh2 );
		float prob = energies[ diffAfter - diffBefore + 10 ];
		
		float acceptance = curand_uniform( &localState );
		
		if ( prob > acceptance )
		{
			exchangeSites( latt, pos1, pos2 );
		}
	}
}



class Metropolis
{
}
;
////////////////////////////////////////////////////////////////////////////////
// declaration, forward
const float R = 1.986;
const float T = 325;
void runTest(int argc, char** argv);
void printLatt( unsigned char *latt, int rowSize, int rowsCount );
long calcSum(unsigned char*, int);
void printDiff(float*, float*, int, int);
float probability( float dG );

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
	runTest(argc, argv);
}

////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
void runTest(int argc, char** argv)
{
		cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
	}
	
	// set seed for rand()
	srand(2005);

	const int LATT_ROW_SIZE = 960;//1790; //894; //5
	const int LATT_ROWS_COUNT = 980;//1792; // 896; //7;
	const long long LATT_SIZE = (long)LATT_ROW_SIZE * (long)LATT_ROWS_COUNT;
	
	cout << LATT_SIZE << " size of the lattice " << endl;
	
	TriangularLattice *latt = new TriangularLattice( LATT_SIZE, LATT_ROW_SIZE, ( LATT_ROW_SIZE * LATT_ROWS_COUNT ) / 2 ); //400512);
	
	const int BLOCK_SIZE = 384; //64;//512; //256; //1;
	const int BLOCKS = 350; // 21;//895; //447; //5;

	assert( BLOCKS * BLOCK_SIZE == ( LATT_SIZE / 7 ) );

	const int OUTPUT_FREQ = 300;
	const int STEPS = 1000;
	
	// allocate device memory for result
	unsigned int size_latt = LATT_SIZE;
	unsigned int mem_size_latt = sizeof(unsigned char) * size_latt;
	unsigned char * d_latt;

	cudaStatus = cudaMalloc((void**)&d_latt, mem_size_latt);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
	}
	
	// allocate host memory for the result
	Npp8u * h_latt = latt->getLattice(); //(int*) malloc(mem_size_latt);

	cout << calcSum(h_latt, size_latt ) << " Sum" << endl;

	//printLatt( h_latt, LATT_ROW_SIZE, LATT_ROWS_COUNT );

	cudaMemcpy(d_latt, h_latt, mem_size_latt, cudaMemcpyHostToDevice);
	
	float omega = 400;
	float energies[36]; // FIXME:36 is too big
	
	for ( int i = 0 ; i < 36 ; ++i )
	{
		energies[ i ] = probability( ( i - 10 ) * omega );
		cout << energies[ i ] << endl;
	}
	
	float* d_energies;
	int mem_size_energies = sizeof( float ) * 36;
	cudaMalloc((void**) &d_energies, mem_size_energies);
	cudaMemcpy(d_energies, energies, mem_size_energies, cudaMemcpyHostToDevice);
	

	/*
	// create and start timer
	unsigned int timer = 0;

	cutCreateTimer(&timer);
	cutStartTimer(timer);
	*/
	
	// setup execution parameters
	dim3 threads( BLOCK_SIZE );
	dim3 grid( BLOCKS );
	
	// Allocate memory for RNG states
	curandState *d_rngStates = 0;
	cudaMalloc((void **)&d_rngStates, grid.x * threads.x * sizeof(curandState));
	
	// Initialise RNG
	srand( (unsigned)time(0) );
	
	//initRNG<<<grid, threads>>>(d_rngStates, rand() );
	
	// execute the kernel
	int start = 0;
	int counter = 0;
	for ( ; start < STEPS * 7 ; ++start )
	{
		int random = rand() % 7;
		
		metropolisStep<<< grid, threads >>>( d_latt, LATT_SIZE, LATT_ROW_SIZE, 	random, d_rngStates, d_energies );
		
		cudaDeviceSynchronize();
		
		if ( start % ( OUTPUT_FREQ * 7 ) == 0 )
		{
		
	
			NppiSize roiSize = {400, 400};
			NppiRect roiRect = {0,0,400,400};
			NppiSize roiSizeSrc = {LATT_ROW_SIZE, LATT_ROWS_COUNT};
			NppiRect roiRectSrc = {0,0,LATT_ROW_SIZE, LATT_ROWS_COUNT};
			ImageNPP_8u_C1 oDeviceDst(roiSize.width, roiSize.height);
			ImageNPP_8u_C1 oDeviceSrc(roiSizeSrc.width, roiSizeSrc.height);

			

			nppiMulC_8u_C1RSfs(d_latt, LATT_ROW_SIZE, 255, oDeviceSrc.data(), oDeviceSrc.pitch(), roiSizeSrc, 0);
			nppiResize_8u_C1R(oDeviceSrc.data(), roiSizeSrc, oDeviceSrc.pitch(), roiRectSrc, oDeviceDst.data(), oDeviceDst.pitch(), roiSize, 400 / (float)LATT_ROW_SIZE, 400 / (float)LATT_ROWS_COUNT, NPPI_INTER_SUPER);

			const double coefs[2][3] = {{1.0, 0.33,0.0},{0.0,1.0,0.0}};

			ImageNPP_8u_C1 oDeviceDst2(oDeviceDst.width() + 0.33*oDeviceDst.height(), oDeviceDst.height());
			NppiSize roiSizeDst2 = {oDeviceDst2.size().nWidth, oDeviceDst2.size().nHeight};
			NppiRect roiWarpedRect = {0,0,oDeviceDst2.size().nWidth, oDeviceDst2.size().nHeight};

			nppiSet_8u_C1R(0,oDeviceDst2.data(), oDeviceDst2.pitch(), roiSizeDst2);

			nppiWarpAffine_8u_C1R(oDeviceDst.data(), roiSize,oDeviceDst.pitch(), roiRect, oDeviceDst2.data(),oDeviceDst2.pitch(), roiWarpedRect,coefs,NPPI_INTER_NN);

			ImageCPU_8u_C1 oHostDst(oDeviceDst2.size());
					// and copy the device result data into it
			oDeviceDst2.copyTo(oHostDst.data(), oHostDst.pitch());
			
			stringstream s;
			s << "d:\\image_"<< counter++ <<".png";
			
			mySaveImage(s.str(), oHostDst);

			/*
			cudaMemcpy(h_latt, d_latt, mem_size_latt, cudaMemcpyDeviceToHost);
			cerr << ( *latt );
			cout << start / 7 << endl;
			*/
		}
	}
	
	// copy result from device to host
	cudaMemcpy(h_latt, d_latt, mem_size_latt, cudaMemcpyDeviceToHost);
	
	// compute reference solution
	//computeGold(reference, h_A, h_B, HA, WA, WB);
	cout << "CUDA res " << calcSum(h_latt, size_latt ) << endl;
	
	delete latt;
	cudaFree(d_latt);
	cudaFree(d_rngStates);
	cudaFree(d_energies);
	
	cudaThreadExit();
}

long calcSum(unsigned char* data, int size)
{
	long sum( 0 );
	for (int i = 0; i < size; ++i)
		sum += data[i];
	return sum;
}


void printLatt( unsigned char *latt, int rowSize, int rowsCount )
{
	for ( int i = 0 ; i < rowsCount ; ++i )
	{
		cout << setw( 2 * ( i + 1 ) ) << " " ;
		for( int j = 0 ; j < rowSize; ++j )
		{
			cout << setw(2) << latt[ i * rowSize + j ] << " ";
		}
		cout << endl;
	}
}
float probability( float dG )
{
	return exp( -dG / ( R * T ) );
}

