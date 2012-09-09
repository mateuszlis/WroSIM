 /** Host code.
 *
 */

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#endif

// project related
#include "MPKK.h"

// CUDA related
#include <cutil_inline.h>
#include <FreeImage.h>
#include <ImageIO.h>
#include <nppi.h>
#include <npp.h>
using namespace npp;
// kernels
#include "Metropolis_kernel.cu"

/**
 * @brief non-class function for saving the image generated from simulation
 * Our policy is that all CUDA related code belongs here, thus this function is
 * not assigned to any class (and is not declared in any header that might be
 * included elsewhere)
 **/
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

MPKK::MPKK( TriangularLattice* latt, double omegaAB, int T, int equilibSteps, unsigned int imageRes ) 
    : Metropolis( latt, omegaAB, T, equilibSteps )
    , IMAGE_RESOLUTION( imageRes )
{
    precalculateEnergies();
}

void MPKK::run( int steps )
{
    TriangularLattice::lattMember *d_latt( NULL ); // device variables are prefixed with d_
    TriangularLattice::lattMember *h_latt( NULL ); // important host variables are prefixed with h_
    float *d_energies( NULL );
    unsigned int blockSize;
    unsigned int blocksCount;
    unsigned int memSizeLatt;

    initDevice( d_latt
              , h_latt
              , memSizeLatt
              , d_energies
              , blockSize
              , blocksCount ); // by reference sets up all variables

    // setup execution parameters
    dim3 threads( blockSize );
    dim3 grid( blocksCount );

    // Allocate memory for RNG states
    curandState *d_rngStates = 0;
    cutilSafeCall( cudaMalloc( ( void ** )&d_rngStates, grid.x * threads.x * sizeof( curandState ) ) );
    // Initialise RNG
    srand( ( unsigned )time( 0 ) );

    initRNG<<<grid, threads>>>( d_rngStates, rand() );

    long long neighHist[7] = { 0, 0, 0, 0, 0, 0, 0 };
    for ( int i = 0; i < steps; i++ )
    {
        for ( int j = 0; j < mStepSize; j++ )
        {
            int random = rand() % 7;

            metropolisStep_kernel<<< grid, threads >>>( d_latt, mpLatt->getLatticeSize(), mpLatt->getRowSize(),
                random, d_rngStates, d_energies );
            cudaDeviceSynchronize();
        }
		if ( i % mOutputFreq == 0 ) //FIXME: export this to a function
        // creates interpolated images
		{
            generateImage( d_latt, i );
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
            cutilSafeCall( cudaMemcpy( h_latt, d_latt, memSizeLatt,
                              cudaMemcpyDeviceToHost) );
        }
        //FIXME: This functionalities go up
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
            for( TriangularLattice::clustersMap::const_iterator it = map.begin() ; it != map.end() ; ++it )
            {
                ( *mpClusterStream ) << i << "\t" << ( *it ).first << "\t" << ( *it ).second << std::endl;
            }
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
    
    // check if kernel execution generated and error
    cutilCheckMsg( "Kernel execution failed" );
    cutilCheckMsg( "cudaThreadSynchronize execution failed" );

    // copy result from device to host
    cutilSafeCall( cudaMemcpy( h_latt, d_latt, memSizeLatt,
                              cudaMemcpyDeviceToHost ) );

    // clean up memory
    cutilSafeCall( cudaFree( d_latt ) );
    cutilSafeCall( cudaFree( d_rngStates ) );
    cutilSafeCall( cudaFree( d_energies ) );

    cudaThreadExit();
}


MPKK::~MPKK()
{
}

/// protected members

long MPKK::calcSum( TriangularLattice::lattMember *data, int size )
{
    long sum( 0 );
    for (int i = 0; i < size; ++i)
        sum += data[i]; 
    return sum;
}

void MPKK::printLatt( TriangularLattice::lattMember *data, int rowSize, int rowsCount )
{
    for ( int i = 0 ; i < rowsCount ; ++i )
    {
        cout << setw( 2 * ( i + 1 ) ) << "  " ;
        for( int j = 0 ; j < rowSize; ++j )
        {   
            cout << setw(2) << data[ i * rowSize + j ] << "  ";
        }
        cout << endl;
    }
}

unsigned int MPKK::calcBlockSize()
{
    cudaDeviceProp deviceProps;
    cudaGetDeviceProperties( &deviceProps, cutGetMaxGflopsDeviceId() );

    if ( mpLatt->getLatticeSize() % 7 != 0 )
    {
        throw std::exception(); // must be divisible by 7
        // TODO: create exception chierarchy
    }
    unsigned int expectedSize = mpLatt->getLatticeSize() / 7;
    unsigned int maxBlocks = deviceProps.maxThreadsPerBlock;

    while ( expectedSize % maxBlocks != 0 )
    {
        maxBlocks /= 2;
    }
    return maxBlocks;
}

void MPKK::precalculateEnergies()
{
    for ( int i = 0 ; i < 36 ; ++i )
    {
        energies[ i ] = prob( ( i - 10 ) * mOmegaAB );
    }
}

/**
 * Most of code here does not obey any code conventions.
 * Dirty hack by Lukasz P. :)
 */
void MPKK::generateImage( TriangularLattice::lattMember *d_latt, int stepNum )
{
    NppiSize roiSize = { IMAGE_RESOLUTION, IMAGE_RESOLUTION };
    NppiRect roiRect = { 0,0, IMAGE_RESOLUTION, IMAGE_RESOLUTION };
    NppiSize roiSizeSrc = { mpLatt->getRowSize(), mpLatt->getLatticeSize() / mpLatt->getRowSize() };
    NppiRect roiRectSrc = { 0,0, roiSizeSrc.width, roiSizeSrc.height };
    ImageNPP_8u_C1 oDeviceDst( roiSize.width, roiSize.height );
    ImageNPP_8u_C1 oDeviceSrc( roiSizeSrc.width, roiSizeSrc.height );

    //nppiMulC_8u_C1RSfs(d_latt, mpLatt->getRowSize(), 1, oDeviceSrc.data(), oDeviceSrc.pitch(), roiSizeSrc, 0);
    nppiCopy_8u_C1R( d_latt, mpLatt->getRowSize(), oDeviceSrc.data(),
        oDeviceSrc.pitch(), roiSizeSrc );
    nppiResize_8u_C1R( oDeviceSrc.data(), roiSizeSrc, oDeviceSrc.pitch(),
        roiRectSrc, oDeviceDst.data(), oDeviceDst.pitch(), roiSize, IMAGE_RESOLUTION /
        ( float )roiSizeSrc.width, IMAGE_RESOLUTION / ( float )roiSizeSrc.height, NPPI_INTER_SUPER );

    const double coefs[2][3] = { { 1.0, 0.33,0.0 },{ 0.0,1.0,0.0 } };

    ImageNPP_8u_C1 oDeviceDst2( oDeviceDst.width() + 0.33*oDeviceDst.height(), oDeviceDst.height() );
    NppiSize roiSizeDst2 = { oDeviceDst2.size().nWidth, oDeviceDst2.size().nHeight };
    NppiRect roiWarpedRect = { 0,0,oDeviceDst2.size().nWidth, oDeviceDst2.size().nHeight };

    nppiSet_8u_C1R( 0,oDeviceDst2.data(), oDeviceDst2.pitch(), roiSizeDst2 );

    nppiWarpAffine_8u_C1R( oDeviceDst.data(), roiSize,oDeviceDst.pitch(), roiRect, oDeviceDst2.data(),oDeviceDst2.pitch(), roiWarpedRect,coefs,NPPI_INTER_NN );

    ImageCPU_8u_C1 oHostDst( oDeviceDst2.size() );
            // and copy the device result data into it
    oDeviceDst2.copyTo( oHostDst.data(), oHostDst.pitch() );
    
    stringstream s;
    s << "image_"<< stepNum <<".png";
    
    mySaveImage( s.str(), oHostDst );
}

void MPKK::initDevice( TriangularLattice::lattMember * &d_latt
               , TriangularLattice::lattMember * &h_latt
               , unsigned int &memSizeLatt
               , float * &d_energies
               , unsigned int &blockSize
               , unsigned int &blocksCount )
{
    cudaSetDevice( cutGetMaxGflopsDeviceId() );
    
    blockSize = calcBlockSize(); 
    blocksCount = ( mpLatt->getLatticeSize() / mStepSize ) / blockSize; 

    // allocate device memory for result
    memSizeLatt = sizeof( TriangularLattice::lattMember ) * mpLatt->getLatticeSize();
    cutilSafeCall( cudaMalloc( ( void** ) &d_latt, memSizeLatt ) );

    h_latt = mpLatt->getLattice();

    cutilSafeCall( cudaMemcpy( d_latt, h_latt, memSizeLatt, cudaMemcpyHostToDevice) );

    int mem_size_energies = sizeof( float ) * PRECALCULATED_ENERGIES_COUNT;
    cutilSafeCall( cudaMalloc( ( void** ) &d_energies, mem_size_energies ) );
    cutilSafeCall( cudaMemcpy( d_energies, energies, mem_size_energies, cudaMemcpyHostToDevice) );

    // Initialise RNG
    srand( ( unsigned )time( 0 ) );
}

