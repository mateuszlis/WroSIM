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

/**
 * @brief Abstract class for performing Metrpopolis based simulations.
 * It has tools for monitoring the simulation (output to the stdout or somewhere else)
 * and basic analysis methods (those which are independent on methodology of simulation.
 **/
class Metropolis
{
    public: // fields
        static const double R = 1.986; //TODO: find place for constants

        const int mEquilibSteps;  ///< number of equilibration steps
        const int mT; ///< temperature of simulation

    public: // member functions

        /**
         * @brief Initializes most of the fields in the class
         *
         * @param latt - lattice to perform the simulation
         * @param omegaAB - energy difference between lipids in lattice
         * @param T - temperature of simulation
         * @param equilibSteps - number of steps for which we don't perform analysis
         **/
        Metropolis( TriangularLattice* latt, double omegaAB = 0.0, int T = 325, int equilibSteps = 0 );
        
        /**
         * @brief Setter for frequency of creating analysis
         *
         * @param freq - frequency (10 - means that every 10 steps there will be analysis step)
         **/
        virtual void setOutputFreq( int freq );

        /**
         * @brief Setter for the energy of interactions between lipids
         **/
        virtual void setOmegaAB( double omegaAB );

        /**
         * @brief Setter for neighbor sites calculator output stream
         **/
        virtual void setNeighOutput( ostream &neighOutputFile );

        /**
         * @brief Setter for output xyz trajectory
         **/
        virtual void setOutput( ostream &frameStream );
        
        /**
         * @brief Setter for status stream (it outputs number of steps performed and so on)
         **/
        virtual void setStatus( ostream &statusStream );

        /**
         * @brief Setter for stream for fraction of first neighbors analysis method
         **/
        virtual void setFNFStream( ostream &fnfStream );

        /**
         * @brief Setter for stream for cluster analysis 
         **/
        virtual void setClusterStream( ostream &clusterStream );

        /**
         * @brief Setter for stream for Mean Square Displacement analysis
         **/
        virtual void setMsdOutput( ostream &MsdOutput );

        /** 
         * @brief Perform simulation for a number of iterations
         * 
         * param steps - number of iterations
         **/
        virtual void run( int steps ) = 0;

        /**
         * @brief default dtor
         **/
        virtual ~Metropolis();

    protected:
        // we use raw pointers instead of std::tr1::function because 
        // we are extremely concerned about the performance
        void (*mpSampler)( TriangularLattice *latt, int &, int & ); 

        int mOutputFreq;

        ostream* mpNeighOutputFile; //FIXME: should not use ostream*
        ostream* mpFNFOutputFile; // first neighbours fraction
        ostream* mpFrameStream;
        ostream* mpStatusStream;
        ostream* mpClusterStream;
        ostream* mpMsdOutputFile;

        bool mIsSetFrameStream; //FIXME: this variables are unnecesary since we can check if output streams are present
        bool mIsSetNeighOutputFile;
        bool mIsSetStatusStream;
        bool mIsSetClusterStream;

        TriangularLattice *mpLatt;
        long long *mpHistArr;
        int mStepSize;
        static const unsigned int LATTICE_NEIGHB_COUNT = 7;
        long long mNeighHist[ LATTICE_NEIGHB_COUNT ];

    protected: // functions

        bool analysisStep( int stepNum ) { return stepNum % mOutputFreq == 0 && stepNum > mEquilibSteps; }
        void performAnalysis( int stepNum );
        void finishAnalysis( int stepNum );
        double prob( double dG );
        double calcEnergyDiff( int pos1, int pos2 );
        void createNeighHist( long long *histArr );
        double calcFirstNeighboursFract();
        double mOmegaAB;

};
