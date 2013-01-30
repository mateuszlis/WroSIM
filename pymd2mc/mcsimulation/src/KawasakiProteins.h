/*
 * KawasakiProteins.h
 *
 *  Created on: 21-10-2012
 *      Author: lisu
 */

#pragma once

// std
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cstdlib>

// project related
#include "KawasakiSimulation.h"
#include "Metropolis.h"
#include "ProteinTriangularLattice.h"

using namespace std;


/**
 * @brief Class performs Kawasak simulations on triangular lattices with two types of lipids and Protein
 **/
class KawasakiProteins : public KawasakiSimulation 
{
    public: // member functions
        /**
         * @brief Creates instance of Kawasaki simulator.
         * 
         * @param latt - lattice to perform simulation (must support proteins)
         * @param omegaAB - energy of interaction between lipids
         * @param T - temperature
         **/
        KawasakiProteins( ProteinTriangularLattice* latt
                        , double omegaAB = 0.0
                        , double omegaAC = 0.0
                        , double omegaBC = 0.0
                        , int T = 325
                        , int equilibSteps = 0 
                        , int proteinStepFreq = 1
                        , int proteinStepSize = 1 );

        /**
         * @brief Starts a simulation for defined number of iterations
         *
         * @param steps - number of iterations
         **/
        virtual void run( int steps );

        /**
         * @brief Default dtor
         **/
        virtual ~KawasakiProteins();

    protected: // fields
        ProteinTriangularLattice* mpLatt;
        const double mOmegaAC;
        const double mOmegaBC;
        const int mProteinStepFreq;
        const int mProteinStepSize;
    protected: // member functions
        /**
         * @brief Performs single step of for proteins
         */
        void proteinStep(); 

        /**
         * @brief checks if given lattMember is a lipid
         **/
        bool isLipid( lattMember site );

        /**
         * @brief performs metropolis step of simulation. Slightly different than in KawasakiSimulation
         **/
        void metropolisStep();

        double calcEnergyDiff( int pos1, int pos2 );
        double calcEnergy();
};
