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
    public:
        static const double R = 1.986; //TODO: find place for constants

        const int mEquilibSteps; 
        const int mT;

        Metropolis( TriangularLattice* latt, double omegaAB = 0.0, int T = 325, int equilibSteps = 0 );
        virtual void setOutputFreq( int freq );
        virtual void setOmegaAB( double omegaAB );
        virtual void setNeighOutput( ostream &neighOutputFile );
        virtual void setOutput( ostream &frameStream );
        virtual void setStatus( ostream &statusStream );
        virtual void setFNFStream( ostream &fnfStream );
        virtual void setClusterStream( ostream &clusterStream );
        virtual void setSampler( Sampler s );
        virtual void run( int steps );
        virtual ~Metropolis();

    protected:
        int mOutputFreq;

        ostream* mpNeighOutputFile; //FIXME: should not use ostream*
        ostream* mpFNFOutputFile; // first neighbours fraction
        ostream* mpFrameStream;
        ostream* mpStatusStream;
        ostream* mpClusterStream;

        bool mIsSetFrameStream; //FIXME: this variables are unnecesary since we can check if output streams are present
        bool mIsSetNeighOutputFile;
        bool mIsSetStatusStream;
        bool mIsSetClusterStream;

        TriangularLattice *mpLatt;
        long long *mpHistArr;
        int mStepSize;

        // we use raw pointers instead of std::tr1::function because 
        // we are extremely concerned about the performance
        void (*mpSampler)( TriangularLattice *latt, int &, int & ); 

        bool analysisStep( int stepNum ) { return stepNum % mOutputFreq == 0 && stepNum > mEquilibSteps; }
        double prob( double dG );
        double calcEnergyDiff( int pos1, int pos2 );
        void metropolisStep();
        void createNeighHist( long long *histArr );
        double calcFirstNeighboursFract();
        double mOmegaAB;

};
