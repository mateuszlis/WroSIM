/*
 * Metropolis.h
 *
 *  Created on: 20-06-2011
 *      Author: lisu
 */

#pragma once
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cstdlib>

#include "TriangularLattice.h"

using namespace std;

enum Sampler
{
    Kawasaki, Almeida, MassiveParallelKawasaki
};

class Metropolis
{
    private:
        double mOmegaAB;
        int mOutputFreq;

        ostream* mpNeighOutputFile;
        ostream* mpFNFOutputFile; // first neighbours fraction
        ostream* mpFrameStream;
        ostream* mpStatusStream;

        bool mIsSetFrameStream;
        bool mIsSetNeighOutputFile;
        bool mIsSetStatusStream;

        TriangularLattice *mpLatt;
        long long *mpHistArr;
        int mStepSize;

        // we use raw pointers instead of std::tr1::function because 
        // we are extremely concerned about the performance
        void (*mpSampler)( TriangularLattice *latt, int &, int & ); 

        double prob( double dG );
        double calcEnergyDiff( int pos1, int pos2 );
        void metropolisStep();
        void createNeighHist( long long *histArr );
        double calcFirstNeighboursFract();

    public:
        static const double R = 1.986;
        static const int EQUIB_STEPS = 1000;
        const int T;

        Metropolis( TriangularLattice* latt, double omegaAB = 0.0, int T = 325 );
        void setOutputFreq( int freq );
        void setOmegaAB( double omegaAB );
        void setNeighOutput( ostream &neighOutputFile );
        void setOutput( ostream &frameStream );
        void setStatus( ostream &statusStream );
        void setFNFStream( ostream &fnfStream );
        void setSampler( Sampler s );
        void run( int steps );
        virtual ~Metropolis();
};
