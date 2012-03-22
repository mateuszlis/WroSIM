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
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <math.h>
#include <assert.h>
using namespace std;


// project related
#include "../types.h"
#include "../TriangularLattice.h"
#include "../Metropolis.h"

// CUDA related
#include <cutil_inline.h>
#include <nppi.h>
#include <npp.h>
#include <FreeImage.h>
#include <ImageIO.h>

using namespace npp;
// kernels
#include "Metropolis_kernel.cu"

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
const float R = 1.986;
const float T = 325;
lattMember* d_latt;
float energies[36]; // FIXME:36 is too big

void runTest(int argc, char** argv);
void printLatt( lattMember *latt, int rowSize, int rowsCount );
long calcSum(lattMember*, int);
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

long calcSum(lattMember* data, int size)
{
    long sum( 0 );
    for (int i = 0; i < size; ++i)
        sum += data[i]; 
    return sum;
}


void printLatt( lattMember *latt, int rowSize, int rowsCount )
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
Metropolis::Metropolis( TriangularLattice* latt, double omegaAB, int T, int equilibSteps) 
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
      , mEquilibSteps( equilibSteps )
      , mT( T )

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

void Metropolis::setClusterStream( ostream &clusterStream )
{
    mpClusterStream = &clusterStream;
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
    unsigned int mem_size_latt = sizeof( TriangularLattice::lattMember ) * mpLatt->getLatticeSize();
    cutilSafeCall(cudaMalloc((void**) &d_latt, mem_size_latt));

    // allocate host memory for the result
    lattMember* h_latt = mpLatt->getLattice(); //(int*) malloc(mem_size_latt);
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
		if ( i % mOutputFreq == 0 ) //FIXME: export this to a function
        // creates interpolated images
		{
            int width = 400;
            int height = 400;
			NppiSize roiSize = { width, height };
			NppiRect roiRect = {0,0, width, height};
			NppiSize roiSizeSrc = { mpLatt->getRowSize(), mpLatt->getLatticeSize() / mpLatt->getRowSize() };
			NppiRect roiRectSrc = { 0,0, roiSizeSrc.width, roiSizeSrc.height };
			ImageNPP_8u_C1 oDeviceDst(roiSize.width, roiSize.height);
			ImageNPP_8u_C1 oDeviceSrc(roiSizeSrc.width, roiSizeSrc.height);

			

			//nppiMulC_8u_C1RSfs(d_latt, mpLatt->getRowSize(), 1, oDeviceSrc.data(), oDeviceSrc.pitch(), roiSizeSrc, 0);
            nppiCopy_8u_C1R( d_latt, mpLatt->getRowSize(), oDeviceSrc.data(),
                oDeviceSrc.pitch(), roiSizeSrc );
			nppiResize_8u_C1R(oDeviceSrc.data(), roiSizeSrc, oDeviceSrc.pitch(),
            roiRectSrc, oDeviceDst.data(), oDeviceDst.pitch(), roiSize, width /
            (float)roiSizeSrc.width, height / (float)roiSizeSrc.height, NPPI_INTER_SUPER);

			const double coefs[2][3] = {{1.0, 0.33,0.0},{0.0,1.0,0.0}};

			ImageNPP_8u_C1 oDeviceDst2(oDeviceDst.width() + 0.33*oDeviceDst.height(), oDeviceDst.height());
			NppiSize roiSizeDst2 = {oDeviceDst2.size().nWidth, oDeviceDst2.size().nHeight};
			NppiRect roiWarpedRect = {0,0,oDeviceDst2.size().nWidth, oDeviceDst2.size().nHeight};

			nppiSet_8u_C1R(0,oDeviceDst2.data(), oDeviceDst2.pitch(), roiSizeDst2);

			nppiWarpAffine_8u_C1R(oDeviceDst.data(), roiSize,oDeviceDst.pitch(), roiRect, oDeviceDst2.data(),oDeviceDst2.pitch(), roiWarpedRect,coefs,NPPI_INTER_NN );

			ImageCPU_8u_C1 oHostDst(oDeviceDst2.size());
					// and copy the device result data into it
			oDeviceDst2.copyTo(oHostDst.data(), oHostDst.pitch());
			
			stringstream s;
			s << "image_"<< i <<".png";
			
			mySaveImage(s.str(), oHostDst);
		}
        //FIXME: currently disabled
        //if ( i % mOutputFreq == 0 && mIsSetFrameStream )
        //{
        //    cutilSafeCall(cudaMemcpy(h_latt, d_latt, mem_size_latt,
        //                      cudaMemcpyDeviceToHost) );
        //    ( *mpFrameStream ) << ( *mpLatt ); //print frame to output
        //}
        if ( analysisStep( i ) )
        {
            cutilSafeCall(cudaMemcpy(h_latt, d_latt, mem_size_latt,
                              cudaMemcpyDeviceToHost) );
        }
        if ( analysisStep( i ) && mIsSetNeighOutputFile )
        {
            createNeighHist( neighHist );
        }
        if ( analysisStep( i ) && mpFNFOutputFile != NULL )
        {
            ( *mpFNFOutputFile ) << setw( 10 ) << i << "\t" << calcFirstNeighboursFract() << endl;
        }

        if ( mIsSetStatusStream )
        {
            ( *mpStatusStream ) << "\r" << i; //print status message
            ( *mpStatusStream ).flush();
        }
        if ( analysisStep( i ) && mpClusterStream != NULL )
        {
            TriangularLattice::clustersMap map;
            mpLatt->calculateClusters( map );
            int sum = 0;
            for( TriangularLattice::clustersMap::const_iterator it = map.begin() ; it != map.end() ; ++it )
            {
                ( *mpClusterStream ) << ( *it ).first << "\t" << ( *it ).second << std::endl;
            }
            cout << sum << " " << map.size()<<  std::endl;
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

