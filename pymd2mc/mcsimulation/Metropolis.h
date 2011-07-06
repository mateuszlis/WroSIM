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
    Kawasaki, Almeida
};

class Metropolis
{
    private:
        int mOutputFreq;
        double mOmegaAB;

        ostream* mpNeighOutputFile;
        ostream* mpFrameStream;
        ostream* mpStatusStream;

        bool mIsSetFrameStream;
        bool mIsSetNeighOutputFile;
        bool mIsSetStatusStream;

        TriangularLattice *mpLatt;
        long long *mpHistArr;
        int mStepSize;

        void (*mpSampler)( TriangularLattice *latt, int &, int & );

        double prob( double dG );
        double calcEnergyDiff( int pos1, int pos2 );
        void metropolisStep();
        void createNeighHist( long long *histArr );

    public:
        static const double R = 1.986;
        static const int T = 310;
        static const int EQUIB_STEPS = 1000;

        Metropolis( TriangularLattice* latt, double omegaAB = 0.0 );
        void setOutputFreq( int freq );
        void setOmegaAB( double omegaAB );
        void setNeighOutput( ostream &neighOutputFile );
        void setOutput( ostream &frameStream );
        void setStatus( ostream &statusStream );
        void setSampler( Sampler s );
        void run( int steps );
        virtual ~Metropolis();
};