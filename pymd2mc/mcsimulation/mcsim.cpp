/*
 * mcsim.cpp
 *
 *  Created on: 23-04-2011
 *      Author: lisu
 */

#include <iostream>

#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>

#include "TriangularLattice.h"
#include "Metropolis.h"

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

        float omega = opt.omega();
        int aLipidsNum = opt.first_particles_count();
        int lattSize = opt.latt_row_size() * opt.latt_row_count();
        int lattRowSize = opt.latt_row_size();
        int steps = opt.steps();
        int outputFreq = opt.output_freq();
        int outputTemperature = opt.T();
        string sampler = opt.sampling();
        string outputFilename = opt.o();
        ofstream outputFile, neighHistFile, fnfFile;
        outputFile.open( outputFilename.c_str());
        neighHistFile.open( "neigh_hist_omega.dat");
        fnfFile.open( "fraction_of_first_neighbors.dat" );

        TriangularLattice *lattice = new TriangularLattice( lattSize, lattRowSize, aLipidsNum);

        Metropolis *simulation = new Metropolis( lattice, omega, outputTemperature );

        simulation->setNeighOutput( neighHistFile );

        simulation->setOutput( outputFile);
        simulation->setFNFStream( fnfFile );
        simulation->setOutputFreq( outputFreq);
        simulation->setStatus( cout);

        if ( sampler == "Kawasaki" )
            simulation->setSampler( Kawasaki );
        if ( sampler == "Almeida" )
            simulation->setSampler( Almeida );
        if ( sampler == "MassiveParallelKawasaki" )
        {
            simulation->setSampler( MassiveParallelKawasaki );
        }

        simulation->run( steps);

        outputFile.close();
        delete lattice;
        delete simulation;

    }
    catch ( const cli::exception& e )
    {
        cerr << e << endl;
        usage();
        return 1;
    }

}
