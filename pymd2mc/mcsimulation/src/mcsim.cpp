/*
 * mcsim.cpp
 *
 *  Created on: 23-04-2011
 *      Author: lisu
 */

// std
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iostream>

// project related
#include "TriangularLattice.h"
#include "LattExchanger.h"
#include "MsdEnabledLattEx.h"
#include "KawasakiProteins.h"
#include "KawasakiSimulation.h"
#include "gpu/MPKK.h"
#include "ProteinTriangularLattice.h"
#include "TriangularLattice.h"

// project constants
#include "mcsim.h"

using namespace std;

#include "mcsimCommandLine.hxx"

const int LATT_SIZE = 10000;

void usage()
{
    cerr << "usage: mcsim [options]" << endl << "options:" << endl;
    options::print_usage( cerr);
}

int main( int argc, char* argv[] )
{
    try
    {
        srand( ( unsigned )time( 0));
        int end; // End of options.
        options opt( argc, argv, end);

        if ( opt.help() )
        {
            usage();
            return 0;
        }

        double omegaAB = opt.omegaAB();
        double omegaAC = opt.omegaAC();
        double omegaBC = opt.omegaBC();
        int proteinsNum = opt.proteins_count();
        int aLipidsNum = opt.first_particles_count();
        int lattSize = opt.latt_row_size() * opt.latt_row_count();
        int lattRowSize = opt.latt_row_size();
        int steps = opt.steps();
        int outputFreq = opt.output_freq();
        int outputTemperature = opt.T();
        int eqSteps = opt.eq_steps();
        int proteinStepFreq = opt.protein_step_freq();
        int proteinStepSize = opt.protein_step_size();
        bool chooseStartRandomly = !opt.no_random_start();
        bool enableMsd = opt.enable_calc_msd();

        string sampler = opt.sampling();
        string outputFilename = opt.o();
        ofstream outputFile, neighHistFile, fnfFile, clusterFile, msdFile;
        outputFile.open( outputFilename.c_str());
        neighHistFile.open( "neigh_hist_omegaAB.dat");
        fnfFile.open( "fraction_of_first_neighbors.dat" );
        clusterFile.open( "clusters.dat" );
        std::cout << eqSteps << std::endl;

        //FIXME: Currently KawasakiProteins is the only option - add command line support
        ProteinTriangularLattice *lattice = new ProteinTriangularLattice( lattSize, lattRowSize, aLipidsNum, proteinsNum, chooseStartRandomly );
        LattExchanger* exchanger = NULL;
        if ( enableMsd && !exchanger )
        {
            std::cout << "Enabled Mean Square Displacement calculation " << std::endl;
            msdFile.open( "msd.dat" );
            exchanger = new MsdEnabledLattEx( lattice, lattice->getProteins(), lattice->getProteinCnt() );
            exchanger->setProteins( lattice->getProteins(), lattice->getProteinCnt() );
            lattice->setExchanger( exchanger );
        }

#ifdef BUILD_CUDA
        std::cout << "CUDA simulation " << omegaAC << " " << omegaBC << "ProteinStepSize " << proteinStepSize << proteinStepFreq << std::endl;
        MPKK *simulation = new MPKK( lattice, omegaAB, outputTemperature, eqSteps );
#else
        KawasakiProteins *simulation = new KawasakiProteins( lattice
                                                           , omegaAB
                                                           , omegaAC
                                                           , omegaBC
                                                           , outputTemperature
                                                           , eqSteps 
                                                           , proteinStepFreq
                                                           , proteinStepSize );
        if ( sampler == "Kawasaki" )
            simulation->setSampler( Kawasaki );
        if ( sampler == "Almeida" )
            simulation->setSampler( Almeida );
        if ( sampler == "MassiveParallelKawasaki" )
        {
            simulation->setSampler( MassiveParallelKawasaki );
        }
#endif
        simulation->setNeighOutput( neighHistFile );

        simulation->setOutput( outputFile );
        simulation->setFNFStream( fnfFile );
        simulation->setMsdOutput( msdFile );
        simulation->setOutputFreq( outputFreq );
        simulation->setStatus( cout );
        simulation->setClusterStream( clusterFile );


        simulation->run( steps );

        outputFile.close();
        clusterFile.close();
        fnfFile.close();
        delete lattice;
        if ( exchanger && enableMsd )
        {
            msdFile.close();
            delete exchanger;
        }
        delete simulation;

    }
    catch ( const cli::exception& e )
    {
        cerr << e << endl;
        usage();
        return 1;
    }

}
