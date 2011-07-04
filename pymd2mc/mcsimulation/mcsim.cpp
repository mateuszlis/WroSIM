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

#include "commandLine.hxx"

const int LATT_SIZE = 10000;

void usage()
{
    cerr << "usage: driver [options] <names>" << endl << "options:" << endl;
    options::print_usage( cerr);
}

int main( int argc, char* argv[] )
{
    try
    {
        int end; // End of options.
        options opt( argc, argv, end);

        if ( opt.help() )
        {
            usage();
            return 0;
        }

        float omega = opt.omega();
        int aLipidsNum = opt.first_particles_count();
        int lattSize = opt.latt_row_size() * opt.latt_row_size();
        int lattRowSize = opt.latt_row_size();
        int steps = opt.steps();
        int outputFreq = opt.output_freq();
        string sampler = opt.sampling();
        string outputFilename = opt.o();
        ofstream outputFile, neighHistFile;
        outputFile.open( outputFilename.c_str());
        neighHistFile.open( "neigh_hist_omega.dat");

        TriangularLattice *lattice = new TriangularLattice( lattSize, lattRowSize, aLipidsNum);

        Metropolis *simulation = new Metropolis( lattice, omega);

        simulation->setOutput( outputFile);
        simulation->setOutputFreq( outputFreq);
        simulation->setStatus( cout);

        if ( sampler == "Kawasaki" )
            simulation->setSampler( Kawasaki);
        if ( sampler == "Almeida" )
            simulation->setSampler( Almeida);

        simulation->run( steps);

        outputFile.close();

    }
    catch ( const cli::exception& e )
    {
        cerr << e << endl;
        usage();
        return 1;
    }
    srand( ( unsigned )time( 0));

}
